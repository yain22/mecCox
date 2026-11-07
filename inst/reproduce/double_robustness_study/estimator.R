# Original MEC-Cox point-estimation workflow with controlled learner inputs.
# The only nuisance-model adaptation is allowing different PS and OR covariate
# sets. All learner, landmark, and calibration helpers are original functions.
# Independent nuisance fits are cached across landmark counts; this is equivalent
# to repeating the original fits because the training folds are unchanged.
suppressPackageStartupMessages(library(survival))
.mec_source <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.mec_folder <- if (is.null(.mec_source)) "." else dirname(.mec_source)
.mec_original <- new.env(parent = globalenv())
sys.source(file.path(.mec_folder, "original_project_functions.R"), .mec_original)

mec_capture <- function(fun) {
  messages <- character()
  value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = value, warning = paste(unique(messages), collapse = " | "),
       error = if (inherits(value, "error")) conditionMessage(value) else "")
}

mec_prepare_nuisances <- function(data, landmark_counts = c(1L, 2L, 5L, 10L, 20L),
                                  fold_seed = 20260427L, Kfold = 10L) {
  vars <- list(Correct = paste0("X", 1:10), Misspecified = "X1")
  grids <- lapply(landmark_counts, function(L)
    .mec_original$choose_landmark_times(data, n_landmarks = L,
                                       landmark_range = c(0.10, 0.90)))
  names(grids) <- paste0("MEC", landmark_counts)
  stopifnot(vapply(grids, length, integer(1)) == landmark_counts)
  times <- sort(unique(unlist(grids)))
  fold_id <- .mec_original$make_stratified_folds(data$A, K = Kfold, seed = fold_seed)
  ehat <- lapply(vars, function(x) rep(NA_real_, nrow(data)))
  surv <- lapply(vars, function(x) matrix(NA_real_, nrow(data), length(times)))
  warnings <- character()
  for (k in seq_len(Kfold)) {
    tr <- data[fold_id != k, , drop = FALSE]
    va <- data[fold_id == k, , drop = FALSE]
    for (kind in names(vars)) {
      ps <- mec_capture(function() .mec_original$estimate_ps_ml(
        data = tr, newdata = va, xvars = vars[[kind]], ps_method = "glm",
        trim = c(0.01, 0.99), auto_tune = FALSE))
      if (nzchar(ps$error)) stop(ps$error)
      or <- mec_capture(function() .mec_original$predict_survival_cox_train(
        train_data = tr, newdata = va, xvars = vars[[kind]], landmark_times = times))
      if (nzchar(or$error)) stop(or$error)
      ehat[[kind]][fold_id == k] <- ps$value$ehat
      surv[[kind]][fold_id == k, ] <- or$value$S
      warnings <- c(warnings, ps$warning, or$warning)
    }
  }
  stopifnot(all(vapply(ehat, function(x) all(is.finite(x)), logical(1))),
            all(vapply(surv, function(x) all(is.finite(x)), logical(1))))
  list(ehat = ehat, survival = surv, grids = grids, times = times,
       fold_id = fold_id, warnings = unique(warnings[nzchar(warnings)]))
}

mec_feature_matrix <- function(survival) {
  H <- cbind(Intercept = 1, survival)
  if (ncol(H) > 1) {
    keep <- c(TRUE, apply(H[, -1, drop = FALSE], 2, sd) > 1e-10)
    H <- H[, keep, drop = FALSE]
  }
  H
}

mec_weighted_cox <- function(data, W) {
  if (length(W) != nrow(data) || any(!is.finite(W)) || any(W <= 0))
    return(list(theta = NA_real_, HR = NA_real_, converged = FALSE,
                reason = "Invalid weights", warning = ""))
  result <- mec_capture(function() coxph(Surv(time, delta) ~ A, data = data,
      weights = W, robust = TRUE, ties = "breslow"))
  if (nzchar(result$error)) return(list(theta = NA_real_, HR = NA_real_,
    converged = FALSE, reason = result$error, warning = result$warning))
  theta <- unname(coef(result$value)["A"])
  bad_warning <- grepl("converg|infinite|singular|iterations|overflow|NaN",
                       result$warning, ignore.case = TRUE)
  ok <- length(theta) == 1L && is.finite(theta) && !bad_warning
  list(theta = theta, HR = exp(theta), converged = ok,
       reason = if (ok) "Converged" else "Cox convergence issue",
       warning = result$warning)
}

mec_original_calibration <- function(H, idx0, idx1, d0) {
  n1 <- length(idx1)
  H0 <- H[idx0, , drop = FALSE]; H1 <- H[idx1, , drop = FALSE]
  # Original coordinates and numerical settings are attempted first.
  result <- mec_capture(function()
    .mec_original$calibrate_weights_bregman_from_design_return(
      Xs = H0, Xp = H1, d = d0, divergence = "kl", alpha = 1 / 2))
  coordinates <- "Original"
  raw_error <- result$error
  if (nzchar(result$error) || !isTRUE(result$value$converged)) {
    # Full square invertible coordinate transformation: no column is discarded,
    # and neither the KL objective nor any original constraint is changed.
    sv <- svd(H)
    rotation <- sweep(sv$v, 2, pmax(sv$d, 1e-5), "/") * sqrt(nrow(H))
    H_numeric <- H %*% rotation
    result <- mec_capture(function()
      .mec_original$calibrate_weights_bregman_from_design_return(
        Xs = H_numeric[idx0, , drop = FALSE],
        Xp = H_numeric[idx1, , drop = FALSE], d = d0,
        divergence = "kl", alpha = 1 / 2, maxit = 1000L, tol = 1e-8))
    coordinates <- "Full invertible SVD"
  }
  result$coordinates <- coordinates
  result$initial_error <- raw_error
  result$original_mean_imbalance <- NA_real_
  result$kkt_residual <- NA_real_
  if (!nzchar(result$error) && all(is.finite(result$value$w)) && all(result$value$w > 0)) {
    w <- result$value$w
    result$original_mean_imbalance <- max(abs(drop(crossprod(H0, w)) - colSums(H1))) / n1
    U <- svd(H0)$u
    log_ratio <- log(w / d0)
    result$kkt_residual <- max(abs(log_ratio - U %*% crossprod(U, log_ratio)))
  }
  result
}

mec_fit_point <- function(data, ehat, H) {
  idx0 <- which(data$A == 0); idx1 <- which(data$A == 1)
  n1 <- length(idx1)
  qhat <- ehat / (1 - ehat)
  d0 <- n1 * qhat[idx0] / sum(qhat[idx0])
  H0 <- H[idx0, , drop = FALSE]; H1 <- H[idx1, , drop = FALSE]
  result <- mec_original_calibration(H, idx0, idx1, d0)
  if (nzchar(result$error)) return(list(theta = NA_real_, HR = NA_real_,
      converged = FALSE, reason = result$error, warning = result$warning,
      calibration_converged = FALSE, calibration_grad_norm = NA_real_,
      original_mean_imbalance = NA_real_, calibration_iterations = NA_integer_))
  cal <- result$value
  W <- numeric(nrow(data)); W[idx1] <- 1; W[idx0] <- cal$w
  fit <- mec_weighted_cox(data, W)
  fit$calibration_converged <- isTRUE(cal$converged)
  fit$calibration_coordinates <- result$coordinates
  fit$calibration_initial_error <- result$initial_error
  fit$calibration_kkt_residual <- result$kkt_residual
  fit$calibration_grad_norm <- cal$grad_norm
  fit$original_mean_imbalance <- max(abs(drop(crossprod(H0, cal$w)) - colSums(H1))) / n1
  fit$calibration_iterations <- cal$iterations
  fit$W <- W
  fit$H <- H
  fit$ehat <- ehat
  fit$warning <- paste(c(result$warning, fit$warning)[nzchar(c(result$warning, fit$warning))], collapse = " | ")
  if (!fit$calibration_converged || !is.finite(fit$original_mean_imbalance) ||
      fit$original_mean_imbalance > 1e-10 || !is.finite(fit$calibration_kkt_residual) ||
      fit$calibration_kkt_residual > 1e-6) {
    fit$converged <- FALSE
    fit$reason <- "Calibration failed convergence, original-moment, or KKT check"
  }
  fit
}

mec_ipw_point <- function(data, xvars) {
  result <- mec_capture(function() .mec_original$make_att_ipw_data_glm(
    data = data, xvars = xvars, trim = c(0.01, 0.99), normalize_controls = TRUE))
  if (nzchar(result$error)) return(list(theta = NA_real_, HR = NA_real_,
      converged = FALSE, reason = result$error, warning = result$warning))
  fit <- mec_weighted_cox(data, result$value$W)
  fit$W <- result$value$W
  fit$warning <- paste(c(result$warning, fit$warning)[nzchar(c(result$warning, fit$warning))], collapse = " | ")
  fit
}

mec_fit_replication <- function(data, config, fold_seed = 20260427L, Kfold = 10L) {
  nuis <- mec_prepare_nuisances(data, fold_seed = fold_seed, Kfold = Kfold)
  vars <- list(Correct = paste0("X", 1:10), Misspecified = "X1")
  # Same original unweighted treatment-only Cox coefficient; no SE is reported.
  naive <- mec_weighted_cox(data, rep(1, nrow(data)))
  ipw <- lapply(vars, function(x) mec_ipw_point(data, x))
  oracle_times <- nuis$grids$MEC5
  oracle <- mec_fit_point(data, dgp_true_ps(data, config),
    mec_feature_matrix(dgp_true_survival(data, oracle_times, config)))
  rows <- list()
  for (PS in names(vars)) for (OR in names(vars)) {
    fits <- list(Naive = naive, IPW = ipw[[PS]])
    for (method in names(nuis$grids)) {
      pos <- match(nuis$grids[[method]], nuis$times)
      H <- mec_feature_matrix(nuis$survival[[OR]][, pos, drop = FALSE])
      fits[[method]] <- mec_fit_point(data, nuis$ehat[[PS]], H)
    }
    fits$Oracle <- oracle
    for (method in names(fits)) {
      fit <- fits[[method]]
      value <- function(key, fallback = NA_real_) if (is.null(fit[[key]])) fallback else fit[[key]]
      rows[[length(rows) + 1L]] <- data.frame(
        PS = PS, OR = OR, method = method, theta = fit$theta, HR = fit$HR,
        true_HR = config$hr, n1 = sum(data$A == 1), n0 = sum(data$A == 0),
        converged = fit$converged, reason = fit$reason, warning = fit$warning,
        calibration_converged = value("calibration_converged", NA),
        calibration_grad_norm = value("calibration_grad_norm"),
        original_mean_imbalance = value("original_mean_imbalance"),
        calibration_iterations = value("calibration_iterations", NA_integer_),
        calibration_coordinates = value("calibration_coordinates", ""),
        calibration_initial_error = value("calibration_initial_error", ""),
        calibration_kkt_residual = value("calibration_kkt_residual"),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

mec_validate_adapter <- function(data, counts = c(2L, 20L), fold_seed = 20260427L, Kfold = 10L) {
  nuis <- mec_prepare_nuisances(data, landmark_counts = counts, fold_seed = fold_seed, Kfold = Kfold)
  ans <- list()
  for (kind in c("Correct", "Misspecified")) for (L in counts) {
    xvars <- if (kind == "Correct") paste0("X", 1:10) else "X1"
    original <- .mec_original$fit_mec_cox(data, xvars = xvars,
      ps_method = "glm", or_method = "cox", use_crossfit = TRUE, Kfold = Kfold,
      n_landmarks = L, landmark_range = c(0.10, 0.90), divergence = "kl",
      trim = c(0.01, 0.99), normalize_baseline = TRUE, auto_tune = FALSE)
    method <- paste0("MEC", L)
    H <- mec_feature_matrix(nuis$survival[[kind]][,
       match(nuis$grids[[method]], nuis$times), drop = FALSE])
    adapter <- mec_fit_point(data, nuis$ehat[[kind]], H)
    ans[[length(ans) + 1L]] <- data.frame(kind = kind, landmarks = L,
      max_ps_difference = max(abs(original$ehat - adapter$ehat)),
      max_feature_difference = max(abs(unname(original$H) - unname(adapter$H))),
      max_weight_difference = max(abs(original$W - adapter$W)),
      absolute_theta_difference = abs(original$theta - adapter$theta),
      original_converged = original$calibration_converged,
      adapter_converged = adapter$calibration_converged)
  }
  ans <- do.call(rbind, ans)
  stopifnot(all(ans$max_ps_difference < 1e-12),
            all(ans$max_feature_difference < 1e-12),
            all(ans$max_weight_difference < 1e-10),
            all(ans$absolute_theta_difference < 1e-10))
  ans
}
