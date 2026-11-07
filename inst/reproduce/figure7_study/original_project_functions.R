# Exact original MEC-Cox function definitions; no study script is executed.
# Supported adapter uses GLM PS and Cox OR with tuning disabled.

expit <- function(x) {
  1 / (1 + exp(-x))
}

clip01 <- function(x, lower = 0.01, upper = 0.99) {
  pmin(pmax(as.numeric(x), lower), upper)
}

solve_safe <- function(A, b = NULL, ridge = 1e-8) {
  A <- as.matrix(A)
  
  ans <- tryCatch(
    {
      if (is.null(b)) {
        solve(A)
      } else {
        solve(A, b)
      }
    },
    error = function(e) {
      A2 <- A + diag(ridge, nrow(A))
      tryCatch(
        {
          if (is.null(b)) {
            solve(A2)
          } else {
            solve(A2, b)
          }
        },
        error = function(e2) {
          if (is.null(b)) {
            MASS::ginv(A)
          } else {
            MASS::ginv(A) %*% b
          }
        }
      )
    }
  )
  
  ans
}

jacobian_fd <- function(f, x, eps = 1e-5) {
  p <- length(x)
  f0 <- f(x)
  J <- matrix(NA_real_, nrow = length(f0), ncol = p)
  
  for (j in seq_len(p)) {
    step <- rep(0, p)
    step[j] <- eps * max(1, abs(x[j]))
    J[, j] <- (f(x + step) - f(x - step)) / (2 * step[j])
  }
  
  J
}

make_stratified_folds <- function(A, K = 5, seed = 20260427) {
  set.seed(seed)
  
  A <- as.integer(A)
  fold_id <- integer(length(A))
  
  for (a in sort(unique(A))) {
    idx <- which(A == a)
    fold_id[idx] <- sample(rep(seq_len(K), length.out = length(idx)))
  }
  
  fold_id
}

make_att_ipw_data_glm <- function(data,
                                  xvars,
                                  trim = c(0.01, 0.99),
                                  normalize_controls = TRUE) {
  d <- as.data.frame(data)
  
  ps_fit <- glm(
    reformulate(xvars, response = "A"),
    data = d,
    family = binomial()
  )
  
  ehat <- clip01(predict(ps_fit, type = "response"), trim[1], trim[2])
  qhat <- ehat / (1 - ehat)
  
  if (normalize_controls) {
    n1 <- sum(d$A == 1)
    qhat[d$A == 0] <- qhat[d$A == 0] * n1 / sum(qhat[d$A == 0])
  }
  
  W <- with(d, A + (1 - A) * qhat)
  
  d$ehat_glm <- ehat
  d$qhat_glm <- qhat
  d$W_att <- W
  
  list(
    data = d,
    ps_fit = ps_fit,
    ehat = ehat,
    qhat = qhat,
    W = W
  )
}

att_weights_from_gamma <- function(gamma,
                                   data,
                                   xvars,
                                   trim = c(0.01, 0.99),
                                   normalize_controls = TRUE) {
  d <- as.data.frame(data)
  Z <- model.matrix(reformulate(xvars, response = "A"), data = d)
  
  ehat <- clip01(expit(as.vector(Z %*% gamma)), trim[1], trim[2])
  qhat <- ehat / (1 - ehat)
  
  if (normalize_controls) {
    n1 <- sum(d$A == 1)
    qhat[d$A == 0] <- qhat[d$A == 0] * n1 / sum(qhat[d$A == 0])
  }
  
  with(d, A + (1 - A) * qhat)
}

cox_score_eta_fixed_weight <- function(theta, data, W) {
  d <- as.data.frame(data)
  
  A <- d$A
  time <- d$time
  delta <- d$delta
  n <- nrow(d)
  
  r <- exp(theta * A)
  
  ord_desc <- order(time, decreasing = TRUE)
  time_desc <- time[ord_desc]
  
  wr_desc <- W[ord_desc] * r[ord_desc]
  wrA_desc <- W[ord_desc] * r[ord_desc] * A[ord_desc]
  
  cum_S0 <- cumsum(wr_desc)
  cum_S1 <- cumsum(wrA_desc)
  
  rle_time <- rle(time_desc)
  times_desc <- rle_time$values
  group_end <- cumsum(rle_time$lengths)
  
  S0_desc <- cum_S0[group_end]
  S1_desc <- cum_S1[group_end]
  
  time_index <- match(time, times_desc)
  
  S0_at <- S0_desc[time_index]
  S1_at <- S1_desc[time_index]
  Abar_at <- S1_at / S0_at
  
  direct <- W * delta * (A - Abar_at)
  Ucox <- sum(direct)
  
  event_times_asc <- sort(unique(time[delta == 1]))
  
  D0_event <- as.numeric(
    tapply(
      W * delta,
      factor(time, levels = event_times_asc),
      sum
    )
  )
  
  event_index_desc <- match(event_times_asc, times_desc)
  
  S0_event <- S0_desc[event_index_desc]
  S1_event <- S1_desc[event_index_desc]
  Abar_event <- S1_event / S0_event
  
  q1_event <- D0_event / S0_event
  qbar_event <- D0_event * Abar_event / S0_event
  
  cum_q1 <- cumsum(q1_event)
  cum_qbar <- cumsum(qbar_event)
  
  pos_event <- findInterval(time, event_times_asc)
  
  Q1_i <- numeric(n)
  Qbar_i <- numeric(n)
  
  has_event_before <- pos_event > 0
  
  Q1_i[has_event_before] <- cum_q1[pos_event[has_event_before]]
  Qbar_i[has_event_before] <- cum_qbar[pos_event[has_event_before]]
  
  risk_correction <- W * r * (A * Q1_i - Qbar_i)
  
  eta <- direct - risk_correction
  
  list(
    U = Ucox,
    eta = eta
  )
}

manual_linwei_se <- function(theta_hat, data, W) {
  ce <- cox_score_eta_fixed_weight(theta_hat, data, W)
  
  D_theta <- as.numeric(
    jacobian_fd(
      function(th) {
        c(cox_score_eta_fixed_weight(th[1], data, W)$U)
      },
      x = c(theta_hat)
    )
  )
  
  phi <- -as.numeric(ce$eta / D_theta)
  
  sqrt(sum(phi^2))
}

fit_unweighted_cox <- function(data) {
  d <- as.data.frame(data)
  
  fit <- coxph(
    Surv(time, delta) ~ A,
    data = d,
    robust = FALSE,
    ties = "breslow"
  )
  
  theta <- unname(coef(fit)["A"])
  
  ## Usual model-based Cox standard error for the unweighted Cox fit.
  se <- sqrt(as.numeric(vcov(fit)["A", "A"]))
  
  list(theta = theta, se = se)
}

fit_att_ipw_naive_linwei_shu <- function(data,
                                         xvars,
                                         trim = c(0.01, 0.99),
                                         normalize_controls = TRUE) {
  obj <- make_att_ipw_data_glm(
    data = data,
    xvars = xvars,
    trim = trim,
    normalize_controls = normalize_controls
  )
  
  d <- obj$data
  W <- obj$W
  
  ## ATT-IPW Cox point estimator.
  ## Naive model-based variance: usual Cox partial-likelihood information
  ## formula, treating the estimated IPW weights as fixed constants.
  fit_naive <- coxph(
    Surv(time, delta) ~ A,
    data = d,
    weights = W,
    robust = FALSE,
    ties = "breslow"
  )
  
  theta_hat <- unname(coef(fit_naive)["A"])
  se_naive <- sqrt(as.numeric(vcov(fit_naive)["A", "A"]))
  
  ## Fixed-weight Lin-Wei/Binder robust sandwich SE.
  se_lw <- manual_linwei_se(theta_hat, d, W)
  
  ## Shu corrected sandwich SE for logistic-GLM PS.
  ps_fit <- obj$ps_fit
  gamma_hat <- coef(ps_fit)
  
  Z <- model.matrix(reformulate(xvars, response = "A"), data = d)
  A <- d$A
  ehat <- obj$ehat
  
  ce <- cox_score_eta_fixed_weight(theta_hat, d, W)
  eta_theta <- ce$eta
  
  score_gamma <- sweep(Z, 1, A - ehat, "*")
  
  D_theta_theta <- as.numeric(
    jacobian_fd(
      function(th) {
        c(cox_score_eta_fixed_weight(th[1], d, W)$U)
      },
      x = c(theta_hat)
    )
  )
  
  D_theta_gamma <- jacobian_fd(
    function(gam) {
      W_gam <- att_weights_from_gamma(
        gamma = gam,
        data = d,
        xvars = xvars,
        trim = trim,
        normalize_controls = normalize_controls
      )
      
      c(cox_score_eta_fixed_weight(theta_hat, d, W_gam)$U)
    },
    x = gamma_hat
  )
  
  D_gamma_gamma <- -crossprod(Z * (ehat * (1 - ehat)), Z)
  
  K <- length(gamma_hat)
  
  D <- rbind(
    c(D_theta_theta, as.numeric(D_theta_gamma)),
    cbind(rep(0, K), D_gamma_gamma)
  )
  
  U_i <- cbind(eta_theta, score_gamma)
  
  IF <- -U_i %*% t(solve_safe(D))
  phi_theta <- IF[, 1]
  
  se_shu <- sqrt(sum(phi_theta^2))
  
  list(
    theta = theta_hat,
    se_naive = se_naive,
    se_lw = se_lw,
    se_shu = se_shu,
    W = W,
    ehat = ehat
  )
}

.to_matrix <- function(X) {
  if (is.null(dim(X))) matrix(as.numeric(X), ncol = 1) else as.matrix(X)
}

.bregman_family <- function(name, alpha = 1/2) {
  name <- tolower(name)
  
  switch(
    name,
    "quadratic" = list(
      g = function(w) w,
      ginv = function(nu) nu,
      gprime = function(w) rep(1, length(w)),
      ok = function(nu) TRUE
    ),
    "kl" = list(
      g = function(w) log(pmax(w, 1e-16)) + 1,
      ginv = function(nu) exp(nu - 1),
      gprime = function(w) 1 / pmax(w, 1e-16),
      ok = function(nu) TRUE
    ),
    "el" = list(
      g = function(w) -1 / pmax(w, 1e-16),
      ginv = function(nu) -1 / nu,
      gprime = function(w) 1 / pmax(w, 1e-16)^2,
      ok = function(nu) all(nu < 0)
    ),
    "hellinger" = list(
      g = function(w) 1 - 1 / sqrt(pmax(w, 1e-16)),
      ginv = function(nu) 1 / (1 - nu)^2,
      gprime = function(w) 0.5 / pmax(w, 1e-16)^(3 / 2),
      ok = function(nu) all(nu < 1)
    ),
    "renyi" = list(
      ## G(w) = w^(alpha + 1) / (alpha + 1), so g(w) = G'(w) = w^alpha.
      ## The default alpha = 1/2 gives g(w) = sqrt(w).
      g = function(w) pmax(w, 1e-16)^alpha,
      ginv = function(nu) pmax(nu, 1e-16)^(1 / alpha),
      gprime = function(w) alpha * pmax(w, 1e-16)^(alpha - 1),
      ok = function(nu) all(nu > 0)
    ),
    stop("Unknown divergence.")
  )
}

.weight_diag <- function(w) {
  n <- length(w)
  s1 <- sum(w)
  s2 <- sum(w^2)
  
  ess <- (s1^2) / s2
  rel_ess <- ess / n
  cvw <- sd(w) / mean(w)
  
  list(
    ESS = ess,
    Rel_ESS = rel_ess,
    W_CV = cvw,
    W_Min = min(w),
    W_Max = max(w),
    SumW = s1
  )
}

calibrate_weights_bregman_from_design_return <- function(Xs,
                                                         Xp,
                                                         d,
                                                         divergence = c("kl", "quadratic", "el", "hellinger", "renyi"),
                                                         alpha = 1/2,
                                                         maxit = 100,
                                                         tol = 1e-10,
                                                         ridge = 1e-8) {
  divergence <- match.arg(divergence)
  
  fam <- .bregman_family(divergence, alpha)
  
  Xs <- .to_matrix(Xs)
  Xp <- .to_matrix(Xp)
  
  Tx <- colSums(Xp)
  
  if (length(d) == 1) d <- rep(d, nrow(Xs))
  
  gd <- fam$g(d)
  lam <- rep(0, ncol(Xs))
  
  step_shrink <- 0.5
  converged <- FALSE
  
  for (it in seq_len(maxit)) {
    nu <- as.numeric(gd + Xs %*% lam)
    
    bt <- 0
    while (!fam$ok(nu) && bt < 30) {
      lam <- lam * step_shrink
      nu <- as.numeric(gd + Xs %*% lam)
      bt <- bt + 1
    }
    
    w <- fam$ginv(nu)
    grad <- as.numeric(crossprod(Xs, w) - Tx)
    grad_norm <- sqrt(sum(grad^2))
    
    if (grad_norm < tol) {
      converged <- TRUE
      break
    }
    
    dwdnu <- 1 / fam$gprime(w)
    Hmat <- crossprod(Xs * dwdnu, Xs) + diag(ridge, ncol(Xs))
    
    step <- tryCatch(
      solve(Hmat, grad),
      error = function(e) MASS::ginv(Hmat) %*% grad
    )
    
    t <- 1
    
    repeat {
      lam_new <- lam - t * as.numeric(step)
      nu_new <- as.numeric(gd + Xs %*% lam_new)
      
      if (fam$ok(nu_new)) {
        lam <- lam_new
        break
      }
      
      t <- t * step_shrink
      
      if (t < 1e-8) {
        lam <- lam_new
        break
      }
    }
  }
  
  nu <- as.numeric(gd + Xs %*% lam)
  w <- fam$ginv(nu)
  
  grad <- as.numeric(crossprod(Xs, w) - Tx)
  grad_norm <- sqrt(sum(grad^2))
  
  dwdnu <- 1 / fam$gprime(w)
  D_lambda_lambda <- crossprod(Xs * dwdnu, Xs)
  
  list(
    w = as.numeric(w),
    lambda = as.numeric(lam),
    grad = grad,
    grad_norm = grad_norm,
    D_lambda_lambda = D_lambda_lambda,
    converged = converged,
    iterations = it,
    diagnostics = .weight_diag(w)
  )
}

estimate_ps_ml <- function(data,
                           xvars,
                           ps_method = c("glm", "nnet", "rf", "bart"),
                           newdata = NULL,
                           trim = c(0.01, 0.99),
                           nnet_size = 4,
                           nnet_decay = 1e-3,
                           nnet_maxit = 200,
                           rf_num_trees = 300,
                           rf_mtry = NULL,
                           rf_min_node_size = 5,
                           bart_ntree = 50,
                           bart_ndpost = 200,
                           bart_nskip = 100,
                           auto_tune = FALSE,
                           tune_n = 100,
                           tune_seed = 20260427,
                           tune_verbose = FALSE) {
  ps_method <- match.arg(ps_method)
  
  d_train <- as.data.frame(data)
  d_pred <- if (is.null(newdata)) d_train else as.data.frame(newdata)
  
  ## Lightweight automatic tuning using a small balanced subset.
  ## We tune only true ML methods; glm is left unchanged.
  tuned_params <- list()
  if (isTRUE(auto_tune) && ps_method %in% c("nnet", "rf", "bart")) {
    tuned_params <- tune_ps_ml(
      data = d_train,
      xvars = xvars,
      ps_method = ps_method,
      trim = trim,
      tune_n = tune_n,
      tune_seed = tune_seed,
      tune_verbose = tune_verbose
    )
    
    if (!is.null(tuned_params$nnet_size)) nnet_size <- tuned_params$nnet_size
    if (!is.null(tuned_params$nnet_decay)) nnet_decay <- tuned_params$nnet_decay
    if (!is.null(tuned_params$rf_num_trees)) rf_num_trees <- tuned_params$rf_num_trees
    if (!is.null(tuned_params$rf_mtry)) rf_mtry <- tuned_params$rf_mtry
    if (!is.null(tuned_params$rf_min_node_size)) rf_min_node_size <- tuned_params$rf_min_node_size
    if (!is.null(tuned_params$bart_ntree)) bart_ntree <- tuned_params$bart_ntree
    if (!is.null(tuned_params$bart_ndpost)) bart_ndpost <- tuned_params$bart_ndpost
    if (!is.null(tuned_params$bart_nskip)) bart_nskip <- tuned_params$bart_nskip
  }
  
  X_train <- as.matrix(d_train[, xvars, drop = FALSE])
  X_pred <- as.matrix(d_pred[, xvars, drop = FALSE])
  A_train <- as.integer(d_train$A)
  
  if (ps_method == "glm") {
    fit <- glm(
      reformulate(xvars, response = "A"),
      data = d_train,
      family = binomial()
    )
    
    ehat <- predict(fit, newdata = d_pred, type = "response")
  }
  
  if (ps_method == "nnet") {
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("Package {nnet} is required for ps_method='nnet'.")
    }
    
    mu <- colMeans(X_train)
    sdv <- apply(X_train, 2, sd)
    sdv[sdv == 0] <- 1
    
    Xs_train <- scale(X_train, center = mu, scale = sdv)
    Xs_pred <- scale(X_pred, center = mu, scale = sdv)
    
    fit <- nnet::nnet(
      x = Xs_train,
      y = A_train,
      size = nnet_size,
      decay = nnet_decay,
      maxit = nnet_maxit,
      entropy = TRUE,
      trace = FALSE,
      MaxNWts = 10000
    )
    
    ehat <- as.numeric(predict(fit, Xs_pred, type = "raw"))
  }
  
  if (ps_method == "rf") {
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("Package {ranger} is required for ps_method='rf'.")
    }
    
    df_rf <- d_train[, c("A", xvars), drop = FALSE]
    df_rf$A <- factor(df_rf$A, levels = c(0, 1))
    
    rf_mtry_use <- if (is.null(rf_mtry)) {
      max(1, floor(sqrt(length(xvars))))
    } else {
      max(1, min(length(xvars), as.integer(rf_mtry)))
    }
    
    fit <- ranger::ranger(
      A ~ .,
      data = df_rf,
      probability = TRUE,
      num.trees = rf_num_trees,
      mtry = rf_mtry_use,
      min.node.size = rf_min_node_size,
      num.threads = 1,
      seed = 20260427
    )
    
    pred <- predict(fit, data = d_pred[, xvars, drop = FALSE])$predictions
    ehat <- pred[, "1"]
  }
  
  if (ps_method == "bart") {
    if (!requireNamespace("dbarts", quietly = TRUE)) {
      stop("Package {dbarts} is required for ps_method='bart'.")
    }
    
    pA <- clip01(mean(A_train), 0.01, 0.99)
    binary_offset <- qnorm(pA)
    
    fit <- tryCatch(
      dbarts::bart(
        x.train = X_train,
        y.train = A_train,
        x.test = X_pred,
        ntree = bart_ntree,
        ndpost = bart_ndpost,
        nskip = bart_nskip,
        binaryOffset = binary_offset,
        verbose = FALSE,
        keeptrees = FALSE,
        nthread = 1
      ),
      error = function(e) {
        dbarts::bart(
          x.train = X_train,
          y.train = A_train,
          x.test = X_pred,
          ntree = bart_ntree,
          ndpost = bart_ndpost,
          nskip = bart_nskip,
          verbose = FALSE,
          keeptrees = FALSE
        )
      }
    )
    
    if (!is.null(fit$prob.test.mean)) {
      ehat <- as.numeric(fit$prob.test.mean)
    } else if (!is.null(fit$prob.test)) {
      ehat <- as.numeric(colMeans(fit$prob.test))
    } else if (!is.null(fit[["yhat.test"]])) {
      # dbarts latent draws already include binaryOffset; transform before averaging.
      ehat <- as.numeric(colMeans(pnorm(fit[["yhat.test"]])))
    } else if (!is.null(fit$yhat.test.mean)) {
      # Mean-only fallback is a plug-in approximation, not a posterior probability mean.
      ehat <- pnorm(as.numeric(fit$yhat.test.mean))
    } else {
      stop("Could not extract BART predicted probabilities.")
    }
  }
  
  ehat <- clip01(ehat, trim[1], trim[2])
  
  list(
    ehat = ehat,
    ps_method = ps_method,
    tuned_params = tuned_params
  )
}

choose_landmark_times <- function(data,
                                  n_landmarks = 4,
                                  landmark_range = c(0.25, 0.70),
                                  landmark_probs = NULL,
                                  landmark_times = NULL) {
  if (!is.null(landmark_times)) {
    return(sort(unique(as.numeric(landmark_times))))
  }
  
  event_times_control <- data$time[data$A == 0 & data$delta == 1]
  
  if (length(event_times_control) < 5) {
    stop("Too few control events to define landmark times.")
  }
  
  if (is.null(landmark_probs)) {
    landmark_probs <- seq(
      landmark_range[1],
      landmark_range[2],
      length.out = n_landmarks
    )
  }
  
  out <- as.numeric(
    quantile(
      event_times_control,
      probs = landmark_probs,
      names = FALSE,
      type = 8
    )
  )
  
  sort(unique(out))
}

predict_survival_cox_train <- function(train_data,
                                       newdata,
                                       xvars,
                                       landmark_times) {
  d_train <- as.data.frame(train_data)
  d_pred <- as.data.frame(newdata)
  d0 <- d_train[d_train$A == 0, , drop = FALSE]
  
  if (sum(d0$delta == 1) < 5) {
    stop("Too few control events in training fold for control-only Cox.")
  }
  
  form <- as.formula(
    paste0("Surv(time, delta) ~ ", paste(xvars, collapse = " + "))
  )
  
  fit <- coxph(
    form,
    data = d0,
    ties = "breslow",
    x = TRUE
  )
  
  bh <- basehaz(fit, centered = FALSE)
  
  H0_t <- approx(
    x = bh$time,
    y = bh$hazard,
    xout = landmark_times,
    method = "constant",
    f = 0,
    rule = 2
  )$y
  
  lp <- as.numeric(
    predict(fit, newdata = d_pred, type = "lp", reference = "zero")
  )
  
  S <- exp(-outer(exp(lp), H0_t))
  S <- pmin(pmax(S, 1e-8), 1 - 1e-8)
  
  colnames(S) <- paste0("S0_t", seq_along(landmark_times))
  
  list(
    S = S,
    fit = fit
  )
}

predict_survival_rsf_train <- function(train_data,
                                       newdata,
                                       xvars,
                                       landmark_times,
                                       num.trees = 300,
                                       mtry = NULL,
                                       min.node.size = 15,
                                       auto_tune = FALSE,
                                       tune_n = 100,
                                       tune_seed = 20260427,
                                       tune_verbose = FALSE) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package {ranger} is required for or_method='rsf'.")
  }
  
  d_train <- as.data.frame(train_data)
  d_pred <- as.data.frame(newdata)
  
  d0 <- d_train[d_train$A == 0, , drop = FALSE]
  
  if (sum(d0$delta == 1) < 5) {
    stop("Too few control events in training fold for RSF.")
  }
  
  tuned_params <- list()
  if (isTRUE(auto_tune)) {
    tuned_params <- tune_rsf_survival(
      train_data = d_train,
      xvars = xvars,
      landmark_times = landmark_times,
      tune_n = tune_n,
      tune_seed = tune_seed,
      tune_verbose = tune_verbose,
      fast_rsf = TRUE
    )
    
    if (!is.null(tuned_params$num.trees)) num.trees <- tuned_params$num.trees
    if (!is.null(tuned_params$mtry)) mtry <- tuned_params$mtry
    if (!is.null(tuned_params$min.node.size)) min.node.size <- tuned_params$min.node.size
  }
  
  df0 <- d0[, c("time", "delta", xvars), drop = FALSE]
  df_pred <- d_pred[, xvars, drop = FALSE]
  
  mtry_use <- if (is.null(mtry)) {
    max(1, floor(sqrt(length(xvars))))
  } else {
    max(1, min(length(xvars), as.integer(mtry)))
  }
  
  fit <- ranger::ranger(
    formula = Surv(time, delta) ~ .,
    data = df0,
    num.trees = num.trees,
    mtry = mtry_use,
    min.node.size = min.node.size,
    #_____________________________
    # Current (slow)
    #splitrule = "logrank",
    # Fast
    splitrule = "extratrees",
    num.random.splits = 1,
    replace = FALSE,
    sample.fraction = 0.632,
    #_____________________________
    write.forest = TRUE,
    num.threads = 1,
    seed = 20260427
  )
  
  pred <- predict(fit, data = df_pred)
  
  surv_mat <- pred$survival
  times <- pred$unique.death.times
  
  if (is.null(times)) {
    times <- fit$unique.death.times
  }
  
  if (is.null(surv_mat) || is.null(times)) {
    stop("RSF prediction failed: survival matrix or death times are NULL.")
  }
  
  ## Vectorized extraction: findInterval gives the last time <= each landmark.
  ## Equivalent to approx(method = "constant", f = 0, rule = 2) but avoids
  ## row-wise R-level looping.
  idx_lm <- findInterval(landmark_times, times)
  idx_lm <- pmin(pmax(idx_lm, 1L), length(times))
  S <- surv_mat[, idx_lm, drop = FALSE]
  
  S <- pmin(pmax(S, 1e-8), 1 - 1e-8)
  colnames(S) <- paste0("S0_t", seq_along(landmark_times))
  
  list(
    S = S,
    fit = fit,
    tuned_params = tuned_params
  )
}

estimate_mec_nuisance <- function(data,
                                  xvars,
                                  ps_method = c("glm", "nnet", "rf", "bart"),
                                  or_method = c("cox", "rsf"),
                                  use_crossfit = TRUE,
                                  Kfold = 5,
                                  n_landmarks = 4,
                                  landmark_range = c(0.25, 0.70),
                                  landmark_probs = NULL,
                                  landmark_times = NULL,
                                  trim = c(0.01, 0.99),
                                  fold_seed = 20260427,
                                  auto_tune = FALSE,
                                  tune_n = 100,
                                  tune_seed = 20260427,
                                  tune_verbose = FALSE) {
  ps_method <- match.arg(ps_method)
  or_method <- match.arg(or_method)
  
  d <- as.data.frame(data)
  n <- nrow(d)
  
  t_landmark <- choose_landmark_times(
    data = d,
    n_landmarks = n_landmarks,
    landmark_range = landmark_range,
    landmark_probs = landmark_probs,
    landmark_times = landmark_times
  )
  
  L <- length(t_landmark)
  
  if (!use_crossfit) {
    ps_obj <- estimate_ps_ml(
      data = d,
      xvars = xvars,
      ps_method = ps_method,
      newdata = d,
      trim = trim,
      auto_tune = auto_tune,
      tune_n = tune_n,
      tune_seed = tune_seed,
      tune_verbose = tune_verbose
    )
    
    if (or_method == "cox") {
      pred <- predict_survival_cox_train(
        train_data = d,
        newdata = d,
        xvars = xvars,
        landmark_times = t_landmark
      )
    } else {
      pred <- predict_survival_rsf_train(
        train_data = d,
        newdata = d,
        xvars = xvars,
        landmark_times = t_landmark,
        auto_tune = auto_tune,
        tune_n = tune_n,
        tune_seed = tune_seed,
        tune_verbose = tune_verbose
      )
    }
    
    H <- cbind(Intercept = 1, pred$S)
    
    if (ncol(H) > 1) {
      keep <- c(TRUE, apply(H[, -1, drop = FALSE], 2, sd) > 1e-10)
      H <- H[, keep, drop = FALSE]
    }
    
    return(
      list(
        ehat = ps_obj$ehat,
        H = H,
        landmark_times = t_landmark,
        use_crossfit = FALSE,
        fold_id = NULL
      )
    )
  }
  
  fold_id <- make_stratified_folds(d$A, K = Kfold, seed = fold_seed)
  
  ehat_cf <- rep(NA_real_, n)
  S_cf <- matrix(NA_real_, nrow = n, ncol = L)
  colnames(S_cf) <- paste0("S0_t", seq_len(L))
  
  for (k in seq_len(Kfold)) {
    idx_valid <- which(fold_id == k)
    idx_train <- which(fold_id != k)
    
    d_train <- d[idx_train, , drop = FALSE]
    d_valid <- d[idx_valid, , drop = FALSE]
    
    ps_obj_k <- estimate_ps_ml(
      data = d_train,
      xvars = xvars,
      ps_method = ps_method,
      newdata = d_valid,
      trim = trim,
      auto_tune = auto_tune,
      tune_n = tune_n,
      tune_seed = tune_seed + k,
      tune_verbose = tune_verbose
    )
    
    ehat_cf[idx_valid] <- ps_obj_k$ehat
    
    if (or_method == "cox") {
      pred_k <- predict_survival_cox_train(
        train_data = d_train,
        newdata = d_valid,
        xvars = xvars,
        landmark_times = t_landmark
      )
    } else {
      pred_k <- predict_survival_rsf_train(
        train_data = d_train,
        newdata = d_valid,
        xvars = xvars,
        landmark_times = t_landmark,
        auto_tune = auto_tune,
        tune_n = tune_n,
        tune_seed = tune_seed + k,
        tune_verbose = tune_verbose
      )
    }
    
    S_cf[idx_valid, ] <- pred_k$S
  }
  
  if (any(!is.finite(ehat_cf))) {
    stop("Cross-fitted PS contains non-finite values.")
  }
  
  if (any(!is.finite(S_cf))) {
    stop("Cross-fitted survival features contain non-finite values.")
  }
  
  H <- cbind(Intercept = 1, S_cf)
  
  if (ncol(H) > 1) {
    keep <- c(TRUE, apply(H[, -1, drop = FALSE], 2, sd) > 1e-10)
    H <- H[, keep, drop = FALSE]
  }
  
  list(
    ehat = ehat_cf,
    H = H,
    landmark_times = t_landmark,
    use_crossfit = TRUE,
    fold_id = fold_id
  )
}

fit_mec_cox <- function(data,
                        xvars,
                        ps_method = c("glm", "nnet", "rf", "bart"),
                        or_method = c("cox", "rsf"),
                        use_crossfit = TRUE,
                        Kfold = 5,
                        n_landmarks = 4,
                        landmark_range = c(0.25, 0.70),
                        landmark_probs = NULL,
                        landmark_times = NULL,
                        divergence = c("kl", "quadratic", "el", "hellinger", "renyi"),
                        trim = c(0.01, 0.99),
                        normalize_baseline = TRUE,
                        auto_tune = FALSE,
                        tune_n = 100,
                        tune_seed = 20260427,
                        tune_verbose = FALSE) {
  ps_method <- match.arg(ps_method)
  or_method <- match.arg(or_method)
  divergence <- match.arg(divergence)
  
  d <- as.data.frame(data)
  n <- nrow(d)
  
  nuis <- estimate_mec_nuisance(
    data = d,
    xvars = xvars,
    ps_method = ps_method,
    or_method = or_method,
    use_crossfit = use_crossfit,
    Kfold = Kfold,
    n_landmarks = n_landmarks,
    landmark_range = landmark_range,
    landmark_probs = landmark_probs,
    landmark_times = landmark_times,
    trim = trim,
    auto_tune = auto_tune,
    tune_n = tune_n,
    tune_seed = tune_seed,
    tune_verbose = tune_verbose
  )
  
  ehat <- nuis$ehat
  H <- nuis$H
  
  qhat <- ehat / (1 - ehat)
  
  idx1 <- which(d$A == 1)
  idx0 <- which(d$A == 0)
  
  n1 <- length(idx1)
  
  d0 <- qhat[idx0]
  
  if (normalize_baseline) {
    d0 <- n1 * d0 / sum(d0)
  }
  
  H1 <- H[idx1, , drop = FALSE]
  H0 <- H[idx0, , drop = FALSE]
  
  cal <- calibrate_weights_bregman_from_design_return(
    Xs = H0,
    Xp = H1,
    d = d0,
    divergence = divergence,
    alpha = 1 / 2
  )
  
  w0_mec <- cal$w
  
  W_mec <- numeric(n)
  W_mec[idx1] <- 1
  W_mec[idx0] <- w0_mec
  
  fit <- coxph(
    Surv(time, delta) ~ A,
    data = d,
    weights = W_mec,
    robust = TRUE,
    ties = "breslow"
  )
  
  theta_hat <- unname(coef(fit)["A"])
  
  ce <- cox_score_eta_fixed_weight(
    theta = theta_hat,
    data = d,
    W = W_mec
  )
  
  eta_theta <- ce$eta
  
  D_theta_theta <- as.numeric(
    jacobian_fd(
      function(th) {
        c(cox_score_eta_fixed_weight(th[1], d, W_mec)$U)
      },
      x = c(theta_hat)
    )
  )
  
  fam <- .bregman_family(divergence, alpha = 1 / 2)
  gd0 <- fam$g(d0)
  
  weights_from_lambda <- function(lambda) {
    nu0 <- as.numeric(gd0 + H0 %*% lambda)
    w0 <- fam$ginv(nu0)
    
    W <- numeric(n)
    W[idx1] <- 1
    W[idx0] <- w0
    W
  }
  
  D_theta_lambda <- jacobian_fd(
    function(lambda) {
      W_lambda <- weights_from_lambda(lambda)
      c(cox_score_eta_fixed_weight(theta_hat, d, W_lambda)$U)
    },
    x = cal$lambda
  )
  
  rho <- matrix(0, nrow = n, ncol = ncol(H))
  rho[idx0, ] <- sweep(H0, 1, w0_mec, "*")
  rho[idx1, ] <- -H1
  
  D_lambda_lambda <- cal$D_lambda_lambda
  K <- ncol(H)
  
  D <- rbind(
    c(D_theta_theta, as.numeric(D_theta_lambda)),
    cbind(rep(0, K), D_lambda_lambda)
  )
  
  U_i <- cbind(eta_theta, rho)
  
  IF <- -U_i %*% t(solve_safe(D))
  phi_theta <- IF[, 1]
  
  se_mec <- sqrt(sum(phi_theta^2))
  
  if (!is.finite(se_mec)) {
    se_mec <- sqrt(as.numeric(vcov(fit)["A", "A"]))
  }
  
  wd <- .weight_diag(w0_mec)
  
  list(
    theta = theta_hat,
    se = se_mec,
    W = W_mec,
    ehat = ehat,
    qhat = qhat,
    H = H,
    lambda = cal$lambda,
    calibration_grad_norm = cal$grad_norm,
    calibration_converged = cal$converged,
    landmark_times = nuis$landmark_times,
    use_crossfit = use_crossfit,
    Kfold = Kfold,
    ESS = wd$ESS,
    Rel_ESS = wd$Rel_ESS,
    W_CV = wd$W_CV,
    W_Min = wd$W_Min,
    W_Max = wd$W_Max
  )
}

