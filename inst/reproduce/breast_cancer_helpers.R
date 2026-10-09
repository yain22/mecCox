# Case-specific MEC-Cox calculations for the public GBSG/Rotterdam example.
# Adapted from the original breast-cancer analysis. These helpers preserve its
# factor coding, including its treatment of grade 1, which occurs only in GBSG.
# They do not change the stricter factor-support checks in the package API.
# All results stay in memory; this file has no top-level analysis or exports.

.breast_with_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(expression)
}

.breast_solve <- function(matrix, rhs = NULL, ridge = 1e-8) {
  matrix <- as.matrix(matrix)
  tryCatch(if (is.null(rhs)) solve(matrix) else solve(matrix, rhs),
    error = function(error) {
      regularized <- matrix + diag(ridge, nrow(matrix))
      tryCatch(if (is.null(rhs)) solve(regularized) else solve(regularized, rhs),
        error = function(error) {
          inverse <- MASS::ginv(matrix)
          if (is.null(rhs)) inverse else inverse %*% rhs
        })
    })
}

.breast_jacobian <- function(function_, point, eps = 1e-5) {
  result <- matrix(NA_real_, length(function_(point)), length(point))
  for (j in seq_along(point)) {
    increment <- rep(0, length(point))
    increment[j] <- eps * max(1, abs(point[j]))
    result[, j] <- (function_(point + increment) -
                     function_(point - increment)) / (2 * increment[j])
  }
  result
}

.breast_folds <- function(source, n_folds, seed) {
  .breast_with_seed(seed, {
    fold <- integer(length(source))
    for (a in c(0L, 1L)) {
      rows <- which(source == a)
      fold[rows] <- sample(rep(seq_len(n_folds), length.out = length(rows)))
    }
    fold
  })
}

.breast_propensity <- function(training, validation, covariates, learner,
                              seed, mlp_epochs) {
  if (learner == "glm") {
    fit <- stats::glm(stats::reformulate(covariates, response = "A"),
                      data = training, family = stats::binomial())
    prediction <- as.numeric(stats::predict(fit, validation, type = "response"))
  } else {
    for (package in c("brulee", "torch")) {
      if (!requireNamespace(package, quietly = TRUE))
        stop(sprintf("Install {%s} for the DL/RSF breast-cancer analysis.",
                     package), call. = FALSE)
    }
    formula <- stats::reformulate(covariates)
    x_train <- stats::model.matrix(formula, training)[, -1L, drop = FALSE]
    x_valid <- stats::model.matrix(formula, validation)[, -1L, drop = FALSE]
    stopifnot(identical(colnames(x_train), colnames(x_valid)))
    center <- colMeans(x_train)
    scale <- apply(x_train, 2L, stats::sd)
    scale[scale == 0] <- 1
    train_frame <- as.data.frame(base::scale(x_train, center, scale))
    valid_frame <- as.data.frame(base::scale(x_valid, center, scale))
    names(train_frame) <- names(valid_frame) <- make.names(
      names(train_frame), unique = TRUE)
    train_frame$A <- factor(training$A, levels = c(0L, 1L))
    previous_threads <- torch::torch_get_num_threads()
    previous_torch_state <- torch::torch_get_rng_state()
    on.exit({
      torch::torch_set_rng_state(previous_torch_state)
      torch::torch_set_num_threads(previous_threads)
    }, add = TRUE)
    torch::torch_set_num_threads(2L)
    torch::torch_manual_seed(as.integer(seed))
    fit <- .breast_with_seed(seed, brulee::brulee_mlp(
      A ~ ., data = train_frame, hidden_units = c(32L, 16L),
      activation = "relu", dropout = 0.1, penalty = 1e-3,
      epochs = as.integer(mlp_epochs), learn_rate = 0.01,
      batch_size = min(256L, nrow(train_frame)), stop_iter = 10L,
      verbose = FALSE))
    prediction <- as.numeric(stats::predict(
      fit, new_data = valid_frame, type = "prob")[[".pred_1"]])
  }
  if (any(!is.finite(prediction)))
    stop("The breast-cancer propensity predictions are not finite.",
         call. = FALSE)
  pmin(pmax(prediction, 0.01), 0.99)
}

.breast_cox_survival <- function(training, validation, covariates, landmarks) {
  controls <- training[training$A == 0L, , drop = FALSE]
  formula <- stats::reformulate(covariates,
                                response = "survival::Surv(time, delta)")
  # Preserve the original factor levels and singular.ok = TRUE convention.
  # Rotterdam has no grade-1 patients. With grade 1 as the factor reference,
  # its grade-2 and grade-3 columns are aliased in the control Cox fit. coxph
  # omits one coefficient and predict.coxph treats that coefficient as zero.
  # Predictions for grade-1 GBSG patients therefore depend on this modeling
  # convention; they are extrapolations, not evidence of covariate overlap.
  fit <- survival::coxph(formula, data = controls, ties = "breslow",
                         x = TRUE, singular.ok = TRUE)
  baseline <- survival::basehaz(fit, centered = FALSE)
  hazard <- stats::approx(baseline$time, baseline$hazard, xout = landmarks,
                          method = "constant", f = 0, rule = 2)$y
  linear_predictor <- as.numeric(stats::predict(
    fit, validation, type = "lp", reference = "zero"))
  prediction <- exp(-outer(exp(linear_predictor), hazard))
  list(survival = pmin(pmax(prediction, 1e-8), 1 - 1e-8),
       aliased_coefficients = names(stats::coef(fit))[is.na(stats::coef(fit))])
}

.breast_tune_rsf <- function(controls, covariates, landmarks, seed) {
  grid <- expand.grid(
    mtry = sort(unique(pmax(1, c(floor(sqrt(length(covariates))),
                                ceiling(length(covariates) / 2))))),
    min.node.size = c(15L, 30L))
  train_rows <- .breast_with_seed(seed, {
    event_rows <- which(controls$delta == 1L)
    censor_rows <- which(controls$delta == 0L)
    sort(c(sample(event_rows, floor(0.6 * length(event_rows))),
           sample(censor_rows, floor(0.6 * length(censor_rows)))))
  })
  validation <- controls[-train_rows, , drop = FALSE]
  if (nrow(validation) < 10L || sum(validation$delta) < 3L)
    return(list(mtry = max(1, floor(sqrt(length(covariates)))),
                min.node.size = 15L))
  training <- controls[train_rows, c("time", "delta", covariates), drop = FALSE]
  concordance <- rep(-Inf, nrow(grid))
  for (j in seq_len(nrow(grid))) {
    concordance[j] <- tryCatch({
      fit <- ranger::ranger(
        survival::Surv(time, delta) ~ ., data = training,
        num.trees = 200L, mtry = grid$mtry[j],
        min.node.size = grid$min.node.size[j], splitrule = "logrank",
        write.forest = TRUE, num.threads = 1L, seed = as.integer(seed))
      prediction <- stats::predict(fit, data = validation[, covariates,
                                                           drop = FALSE])
      times <- prediction$unique.death.times
      if (is.null(times)) times <- fit$unique.death.times
      position <- min(max(findInterval(stats::median(landmarks), times), 1L),
                      length(times))
      risk <- 1 - as.numeric(prediction$survival[, position])
      value <- survival::concordance(
        survival::Surv(time, delta) ~ risk, data = validation)$concordance
      max(value, 1 - value)
    }, error = function(error) -Inf)
  }
  if (!any(is.finite(concordance)))
    stop("No RSF tuning candidate succeeded in the breast-cancer analysis.",
         call. = FALSE)
  selected <- which.max(concordance)
  list(mtry = grid$mtry[selected], min.node.size = grid$min.node.size[selected])
}

.breast_rsf_survival <- function(training, validation, covariates, landmarks,
                                seed, tuning_seed, auto_tune) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("Install {ranger} for the RSF breast-cancer analysis.", call. = FALSE)
  controls <- training[training$A == 0L, , drop = FALSE]
  tuning <- if (auto_tune) {
    .breast_tune_rsf(controls, covariates, landmarks, tuning_seed)
  } else list(mtry = max(1, floor(sqrt(length(covariates)))),
              min.node.size = 15L)
  # Keep all factor levels as in the original analysis. Grade 1 has no
  # external-control support; its forest prediction also extrapolates.
  fit <- ranger::ranger(
    survival::Surv(time, delta) ~ .,
    data = controls[, c("time", "delta", covariates), drop = FALSE],
    num.trees = if (auto_tune) 500L else 200L,
    mtry = tuning$mtry, min.node.size = tuning$min.node.size,
    splitrule = "logrank", write.forest = TRUE, num.threads = 1L,
    seed = as.integer(seed))
  prediction <- stats::predict(fit, data = validation[, covariates, drop = FALSE])
  times <- prediction$unique.death.times
  if (is.null(times)) times <- fit$unique.death.times
  positions <- pmin(pmax(findInterval(landmarks, times), 1L), length(times))
  list(survival = pmin(pmax(prediction$survival[, positions, drop = FALSE],
                            1e-8), 1 - 1e-8), tuning = tuning)
}

.breast_cox_contributions <- function(theta, data, weights) {
  source <- data$A
  relative_risk <- exp(theta * source)
  order <- order(data$time, decreasing = TRUE)
  distinct_times <- rle(data$time[order])
  group_end <- cumsum(distinct_times$lengths)
  s0 <- cumsum((weights * relative_risk)[order])[group_end]
  s1 <- cumsum((weights * relative_risk * source)[order])[group_end]
  index <- match(data$time, distinct_times$values)
  source_mean <- s1[index] / s0[index]
  direct <- weights * data$delta * (source - source_mean)
  event_times <- sort(unique(data$time[data$delta == 1L]))
  weighted_events <- as.numeric(tapply(weights * data$delta,
    factor(data$time, levels = event_times), sum))
  event_index <- match(event_times, distinct_times$values)
  q1 <- cumsum(weighted_events / s0[event_index])
  qbar <- cumsum(weighted_events * s1[event_index] / s0[event_index]^2)
  position <- findInterval(data$time, event_times)
  correction1 <- correction_bar <- numeric(nrow(data))
  present <- position > 0L
  correction1[present] <- q1[position[present]]
  correction_bar[present] <- qbar[position[present]]
  list(score = sum(direct), contribution = direct -
         weights * relative_risk * (source * correction1 - correction_bar))
}

.breast_calibrate <- function(control_basis, target_basis, baseline) {
  target <- colSums(target_basis)
  lambda <- numeric(ncol(control_basis))
  converged <- FALSE
  # Newton updates with backtracking minimize the same KL dual objective as
  # the original script, with a safeguard against an overflowing full step.
  objective <- function(candidate) {
    weights <- baseline * exp(as.numeric(control_basis %*% candidate))
    if (any(!is.finite(weights)) || any(weights <= 0)) return(Inf)
    sum(weights) - sum(target * candidate)
  }
  for (iteration in seq_len(200L)) {
    weights <- baseline * exp(as.numeric(control_basis %*% lambda))
    gradient <- as.numeric(crossprod(control_basis, weights) - target)
    if (max(abs(gradient)) / nrow(target_basis) <= 1e-8) {
      converged <- TRUE
      break
    }
    hessian <- crossprod(control_basis * weights, control_basis)
    direction <- .breast_solve(hessian + diag(1e-8, ncol(hessian)), gradient)
    old_value <- objective(lambda)
    step <- 1
    accepted <- FALSE
    for (attempt in seq_len(40L)) {
      candidate <- lambda - step * as.numeric(direction)
      value <- objective(candidate)
      if (is.finite(value) && value <= old_value -
          1e-4 * step * sum(gradient * direction)) {
        lambda <- candidate
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) break
  }
  weights <- baseline * exp(as.numeric(control_basis %*% lambda))
  residual <- as.numeric(crossprod(control_basis, weights) - target)
  maximum_gap <- max(abs(residual)) / nrow(target_basis)
  if (!is.finite(maximum_gap) || maximum_gap > 1e-8)
    stop(sprintf("Breast-cancer calibration failed (mean residual %.4g).",
                 maximum_gap), call. = FALSE)
  list(weights = weights, lambda = lambda,
       hessian = crossprod(control_basis * weights, control_basis),
       max_mean_balance_residual = maximum_gap, iterations = iteration,
       converged = TRUE)
}

fit_breast_mec <- function(data, covariates, ps_learner = c("glm", "mlp"),
                           survival_learner = c("cox", "rsf"),
                           seed = 20260427L, n_folds = 10L,
                           n_landmarks = 20L, mlp_epochs = 200L,
                           rsf_auto_tune = TRUE) {
  ps_learner <- match.arg(ps_learner)
  survival_learner <- match.arg(survival_learner)
  data <- as.data.frame(data)
  required <- c("time", "delta", "A", covariates)
  if (!all(required %in% names(data)) || anyNA(data[, required, drop = FALSE]))
    stop("Breast-cancer analysis columns must be present and complete.",
         call. = FALSE)
  if (!all(data$A %in% c(0, 1)) || !all(data$delta %in% c(0, 1)) ||
      any(!is.finite(data$time)) || any(data$time <= 0))
    stop("Invalid source, event, or follow-up values.", call. = FALSE)
  controls <- which(data$A == 0L)
  treated <- which(data$A == 1L)
  if (n_folds < 2L || n_folds != as.integer(n_folds) ||
      min(length(controls), length(treated)) < n_folds)
    stop("Each cohort must have at least n_folds observations.", call. = FALSE)
  landmarks <- sort(unique(as.numeric(stats::quantile(
    data$time[data$A == 0L & data$delta == 1L],
    seq(0.10, 0.90, length.out = n_landmarks), names = FALSE, type = 8))))
  fold <- .breast_folds(data$A, n_folds, seed)
  propensity <- rep(NA_real_, nrow(data))
  survival <- matrix(NA_real_, nrow(data), length(landmarks))
  learner_diagnostics <- vector("list", n_folds)
  for (k in seq_len(n_folds)) {
    message(sprintf("  Cross-fitting fold %d of %d (%s/%s)",
                    k, n_folds, toupper(ps_learner), toupper(survival_learner)))
    training <- data[fold != k, , drop = FALSE]
    validation <- data[fold == k, , drop = FALSE]
    propensity[fold == k] <- .breast_propensity(
      training, validation, covariates, ps_learner, seed + k, mlp_epochs)
    prediction <- if (survival_learner == "cox") {
      .breast_cox_survival(training, validation, covariates, landmarks)
    } else {
      .breast_rsf_survival(training, validation, covariates, landmarks,
                           seed, seed + k, rsf_auto_tune)
    }
    survival[fold == k, ] <- prediction$survival
    learner_diagnostics[[k]] <- prediction[setdiff(names(prediction), "survival")]
  }
  if (any(!is.finite(survival)))
    stop("Control-survival predictions are not finite.", call. = FALSE)
  basis <- cbind(intercept = 1, survival)
  basis <- basis[, c(TRUE, apply(survival, 2L, stats::sd) > 1e-10), drop = FALSE]
  baseline <- propensity[controls] / (1 - propensity[controls])
  baseline <- baseline * length(treated) / sum(baseline)
  h0 <- basis[controls, , drop = FALSE]
  h1 <- basis[treated, , drop = FALSE]
  calibration <- .breast_calibrate(h0, h1, baseline)
  weights <- rep(1, nrow(data))
  weights[controls] <- calibration$weights
  fit <- survival::coxph(survival::Surv(time, delta) ~ A, data = data,
                         weights = weights, robust = TRUE, ties = "breslow")
  theta <- unname(stats::coef(fit)["A"])
  contribution <- .breast_cox_contributions(theta, data, weights)$contribution
  derivative_theta <- as.numeric(.breast_jacobian(function(value)
    .breast_cox_contributions(value[1L], data, weights)$score, theta))
  derivative_lambda <- .breast_jacobian(function(value) {
    perturbed <- weights
    perturbed[controls] <- baseline * exp(as.numeric(h0 %*% value))
    .breast_cox_contributions(theta, data, perturbed)$score
  }, calibration$lambda)
  rho <- matrix(0, nrow(data), ncol(basis))
  rho[controls, ] <- h0 * weights[controls]
  rho[treated, ] <- -h1
  derivative <- rbind(c(derivative_theta, as.numeric(derivative_lambda)),
    cbind(rep(0, ncol(basis)), calibration$hessian))
  influence <- -cbind(contribution, rho) %*% t(.breast_solve(derivative))
  se <- sqrt(sum(influence[, 1L]^2))
  if (!is.finite(theta) || !is.finite(se) || se <= 0)
    stop("The breast-cancer estimate or stacked-sandwich SE is invalid.",
         call. = FALSE)
  list(theta = theta, se = se, hr = exp(theta),
       conf_int = exp(theta + c(-1, 1) * stats::qnorm(0.975) * se),
       weights = weights, baseline_control_weights = baseline,
       propensity = propensity, landmarks = landmarks, fold_id = fold,
       basis = basis, calibration = calibration,
       learner_diagnostics = learner_diagnostics,
       diagnostics = list(ess = sum(weights[controls])^2 /
                            sum(weights[controls]^2),
                          cv = stats::sd(weights[controls]) /
                            mean(weights[controls])))
}
