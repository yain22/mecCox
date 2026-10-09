# Manuscript simulation functions: study implementation.
# Function definitions from Sim_MEC_Cox_landmark_feature.R.
# Sourcing this file defines functions only; see the adjacent adapters.
# BART uses the study conversion of latent means with its range heuristic
# and binaryOffset calculation. The calculation is described in README.md.
# See README.md. These functions do not replace mecCox::fit_mec_cox.

expit <- function(x) {
  1 / (1 + exp(-x))
}

clip01 <- function(x, lower = 0.01, upper = 0.99) {
  pmin(pmax(as.numeric(x), lower), upper)
}

round_numeric_df <- function(x, digits = 4) {
  x <- as.data.frame(x)
  num_cols <- vapply(x, is.numeric, logical(1))
  x[num_cols] <- lapply(x[num_cols], round, digits = digits)
  x
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

source_score <- function(X,
                         alpha0 = -0.2,
                         ps_truth = c("linear", "nonlinear"),
                         ps_nonlinearity = 1) {
  ps_truth <- match.arg(ps_truth)
  X <- as.matrix(X)

  ## The first five covariates drive source selection.
  ## This keeps the PS model interpretable but allows a scalar nonlinearity
  ## parameter to control the difficulty of estimating ATT weights.
  x1 <- X[, 1]
  x2 <- X[, 2]
  x3 <- X[, 3]
  x4 <- X[, 4]
  x5 <- X[, 5]

  ## Linear component: correctly specified by logistic-GLM PS.
  #lin_eta <- 0.35 * rowSums(X[, 1:5, drop = FALSE])
  lin_eta <- 0.75 * x1 + 0.75 * x2 + 0.65 * x3 + 0.65 * x4 + 0.55 * x5

  ## Nonlinear component: intentionally difficult for logistic-GLM PS,
  ## but learnable by flexible ML learners such as BART/RF/nnet.
  ## Terms are roughly centered to avoid changing the marginal treated-control
  ## ratio too aggressively when ps_nonlinearity changes.
  nonlin_eta <-
    0.70 * sin(1.25 * x1) +
    0.45 * (x2^2 - 1) -
    0.55 * (as.numeric(x3 > 0) - 0.5) +
    0.35 * x4 * x5 +
    0.25 * (cos(x1 + x2) - exp(-1))

  if (ps_truth == "linear") {
    eta <- alpha0 + lin_eta
  } else {
    eta <- alpha0 + lin_eta + ps_nonlinearity * nonlin_eta
  }

  clip01(expit(eta), 0.02, 0.98)
}

sample_source_covariates <- function(n1,
                                     n0,
                                     p = 10,
                                     alpha0 = -0.2,
                                     ps_truth = "linear",
                                     ps_nonlinearity = 1) {
  X1 <- matrix(NA_real_, nrow = 0, ncol = p)
  X0 <- matrix(NA_real_, nrow = 0, ncol = p)

  while (nrow(X1) < n1 || nrow(X0) < n0) {
    m <- max(4000, 4 * (n1 + n0))

    Xcand <- matrix(rnorm(m * p), nrow = m, ncol = p)
    pi <- source_score(
      Xcand,
      alpha0 = alpha0,
      ps_truth = ps_truth,
      ps_nonlinearity = ps_nonlinearity
    )
    A <- rbinom(m, size = 1, prob = pi)

    add1 <- Xcand[A == 1, , drop = FALSE]
    add0 <- Xcand[A == 0, , drop = FALSE]

    if (nrow(X1) < n1) X1 <- rbind(X1, add1)
    if (nrow(X0) < n0) X0 <- rbind(X0, add0)
  }

  X1 <- X1[seq_len(n1), , drop = FALSE]
  X0 <- X0[seq_len(n0), , drop = FALSE]

  colnames(X1) <- paste0("X", seq_len(p))
  colnames(X0) <- paste0("X", seq_len(p))

  list(X1 = X1, X0 = X0)
}

domain_shift <- function(X, type = c("constant", "nonlinear")) {
  type <- match.arg(type)
  X <- as.matrix(X)

  if (type == "constant") {
    return(rep(1, nrow(X)))
  }

  1 + rowSums(sin(2 * X[, 1:5, drop = FALSE])) / sqrt(5)
}

prognostic_lp <- function(X,
                          outcome_truth = c("linear_ph", "nonlinear_ph"),
                          or_nonlinearity = 1,
                          p_signal_outcome = 10) {
  outcome_truth <- match.arg(outcome_truth)
  X <- as.matrix(X)
  p <- ncol(X)

  if (p < 5) {
    stop("p must be at least 5 because the nonlinear prognostic component uses X1,...,X5.")
  }

  ## True nonzero outcome coefficients.
  ## Only the first p_signal_outcome variables can affect the true outcome.
  b_base <- c(
    log(1.75),
    log(1.75),
    log(1.60),
    log(1.60),
    log(1.50),
    log(1.25),
    log(1.25),
    log(1.25),
    log(1.25),
    log(1.25)
  )

  ## If p > 10, extra variables are pure noise in the true outcome model.
  ## If p < 10, use only the available coefficients.
  b <- rep(0, p)
  J <- min(length(b_base), p_signal_outcome, p)
  b[seq_len(J)] <- b_base[seq_len(J)]

  lin_lp <- drop(X %*% b)

  x1 <- X[, 1]
  x2 <- X[, 2]
  x3 <- X[, 3]
  x4 <- X[, 4]
  x5 <- X[, 5]

  nonlin_lp <-
    0.45 * sin(x2) +
    0.35 * (x3^2 - 1) +
    0.30 * (as.numeric(x4 > 0) - 0.5) +
    0.25 * x1 * x5 +
    0.20 * (cos(x2 + x5) - exp(-1))

  if (outcome_truth == "linear_ph") {
    lp <- lin_lp
  } else {
    lp <- lin_lp + or_nonlinearity * nonlin_lp
  }

  lp
}

weibull_time <- function(lp, lambda0 = 0.00008, eta_shape = 2.0) {
  U <- runif(length(lp))
  (-log(U) / (lambda0 * exp(lp)))^(1 / eta_shape)
}

generate_observed_data <- function(n1 = 250,
                                   n0 = 500,
                                   p = 10,
                                   alpha0 = -0.2,
                                   ps_truth = "linear",
                                   outcome_truth = "linear_ph",
                                   ps_nonlinearity = 1,
                                   or_nonlinearity = 1,
                                   beta_cond = log(0.70),
                                   xi = 0,
                                   h_type = "constant",
                                   lambda0 = 0.00008,
                                   eta_shape = 2.0,
                                   censor_rate = 0.0008) {
  SX <- sample_source_covariates(
    n1 = n1,
    n0 = n0,
    p = p,
    alpha0 = alpha0,
    ps_truth = ps_truth,
    ps_nonlinearity = ps_nonlinearity
  )

  X1 <- SX$X1
  X0 <- SX$X0

  lp0_trial_ext <- prognostic_lp(
    X1,
    outcome_truth = outcome_truth,
    or_nonlinearity = or_nonlinearity
  )
  lp0_trial_true <- lp0_trial_ext + xi * domain_shift(X1, h_type)
  lp1_trial_true <- lp0_trial_true + beta_cond

  lp0_external <- prognostic_lp(
    X0,
    outcome_truth = outcome_truth,
    or_nonlinearity = or_nonlinearity
  )

  T1_trial <- weibull_time(lp1_trial_true, lambda0, eta_shape)
  T0_external <- weibull_time(lp0_external, lambda0, eta_shape)

  C1 <- rexp(n1, rate = censor_rate)
  C0 <- rexp(n0, rate = censor_rate)

  dat1 <- data.frame(
    A = 1L,
    time = pmin(T1_trial, C1),
    delta = as.integer(T1_trial <= C1),
    X1
  )

  dat0 <- data.frame(
    A = 0L,
    time = pmin(T0_external, C0),
    delta = as.integer(T0_external <= C0),
    X0
  )

  dat <- rbind(dat1, dat0)
  dat$id <- seq_len(nrow(dat))

  dat
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

make_balanced_tuning_indices <- function(A,
                                         tune_n = 100) {
  A <- as.integer(A)
  n <- length(A)

  if (n <= tune_n) {
    return(seq_len(n))
  }

  idx1 <- which(A == 1)
  idx0 <- which(A == 0)

  if (length(idx1) == 0 || length(idx0) == 0) {
    return(seq_len(min(tune_n, n)))
  }

  n1_target <- min(length(idx1), floor(tune_n / 2))
  n0_target <- min(length(idx0), tune_n - n1_target)

  ## Fill any leftover capacity from the other group.
  left <- tune_n - n1_target - n0_target
  if (left > 0 && length(idx1) > n1_target) {
    add <- min(left, length(idx1) - n1_target)
    n1_target <- n1_target + add
    left <- left - add
  }
  if (left > 0 && length(idx0) > n0_target) {
    add <- min(left, length(idx0) - n0_target)
    n0_target <- n0_target + add
  }

  sort(c(head(idx1, n1_target), head(idx0, n0_target)))
}

make_stratified_train_valid_split <- function(A,
                                              train_frac = 0.70,
                                              seed = 20260427) {
  set.seed(seed)
  A <- as.integer(A)
  n <- length(A)

  idx_train <- integer(0)
  idx_valid <- integer(0)

  for (a in sort(unique(A))) {
    idx_a <- which(A == a)
    if (length(idx_a) <= 2) {
      idx_train <- c(idx_train, idx_a)
    } else {
      n_train_a <- max(1, floor(train_frac * length(idx_a)))
      train_a <- sample(idx_a, size = n_train_a, replace = FALSE)
      valid_a <- setdiff(idx_a, train_a)
      idx_train <- c(idx_train, train_a)
      idx_valid <- c(idx_valid, valid_a)
    }
  }

  if (length(idx_valid) == 0) {
    idx_all <- seq_len(n)
    n_train <- max(1, floor(train_frac * n))
    idx_train <- sample(idx_all, size = n_train, replace = FALSE)
    idx_valid <- setdiff(idx_all, idx_train)
  }

  list(
    train = sort(idx_train),
    valid = sort(idx_valid)
  )
}

binary_logloss <- function(y,
                           p,
                           eps = 1e-6) {
  p <- pmin(pmax(as.numeric(p), eps), 1 - eps)
  y <- as.integer(y)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

safe_metric_min <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) {
    return(NA_integer_)
  }
  which.min(ifelse(is.finite(x), x, Inf))
}

safe_metric_max <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) {
    return(NA_integer_)
  }
  which.max(ifelse(is.finite(x), x, -Inf))
}

tune_ps_ml <- function(data,
                       xvars,
                       ps_method = c("glm", "nnet", "rf", "bart"),
                       trim = c(0.01, 0.99),
                       tune_n = 100,
                       tune_seed = 20260427,
                       tune_verbose = FALSE) {
  ps_method <- match.arg(ps_method)

  if (ps_method == "glm") {
    return(list())
  }

  d <- as.data.frame(data)
  idx_tune <- make_balanced_tuning_indices(d$A, tune_n = tune_n)
  d_tune <- d[idx_tune, , drop = FALSE]

  if (length(unique(d_tune$A)) < 2 || nrow(d_tune) < 20) {
    if (isTRUE(tune_verbose)) {
      message("PS auto-tuning skipped: too few observations or only one class.")
    }
    return(list())
  }

  sp <- make_stratified_train_valid_split(
    A = d_tune$A,
    train_frac = 0.70,
    seed = tune_seed
  )

  d_train <- d_tune[sp$train, , drop = FALSE]
  d_valid <- d_tune[sp$valid, , drop = FALSE]

  if (length(unique(d_train$A)) < 2 || length(unique(d_valid$A)) < 2) {
    if (isTRUE(tune_verbose)) {
      message("PS auto-tuning skipped: train/validation split has one class.")
    }
    return(list())
  }

  p_dim <- length(xvars)

  if (ps_method == "nnet") {
    grid <- expand.grid(
      nnet_size = c(2, 4, 6),
      nnet_decay = c(1e-4, 1e-3, 1e-2),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (ps_method == "rf") {
    grid <- expand.grid(
      rf_mtry = sort(unique(pmax(1, c(floor(sqrt(p_dim)), ceiling(p_dim / 2), p_dim)))),
      rf_min_node_size = c(5, 10, 20),
      rf_num_trees = c(200),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (ps_method == "bart") {
    grid <- expand.grid(
      bart_ntree = c(25, 50, 100),
      bart_ndpost = c(100),
      bart_nskip = c(50),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }

  loss <- rep(Inf, nrow(grid))

  for (gg in seq_len(nrow(grid))) {
    pars <- as.list(grid[gg, , drop = FALSE])

    ans <- tryCatch(
      {
        pred <- estimate_ps_ml(
          data = d_train,
          xvars = xvars,
          ps_method = ps_method,
          newdata = d_valid,
          trim = trim,
          nnet_size = if (!is.null(pars$nnet_size)) pars$nnet_size else 4,
          nnet_decay = if (!is.null(pars$nnet_decay)) pars$nnet_decay else 1e-3,
          rf_num_trees = if (!is.null(pars$rf_num_trees)) pars$rf_num_trees else 300,
          rf_mtry = if (!is.null(pars$rf_mtry)) pars$rf_mtry else NULL,
          rf_min_node_size = if (!is.null(pars$rf_min_node_size)) pars$rf_min_node_size else 5,
          bart_ntree = if (!is.null(pars$bart_ntree)) pars$bart_ntree else 50,
          bart_ndpost = if (!is.null(pars$bart_ndpost)) pars$bart_ndpost else 200,
          bart_nskip = if (!is.null(pars$bart_nskip)) pars$bart_nskip else 100,
          auto_tune = FALSE
        )
        binary_logloss(d_valid$A, pred$ehat)
      },
      error = function(e) Inf
    )

    loss[gg] <- ans
  }

  best_id <- safe_metric_min(loss)
  if (is.na(best_id)) {
    if (isTRUE(tune_verbose)) {
      message("PS auto-tuning failed for all candidates; using defaults.")
    }
    return(list())
  }

  best <- as.list(grid[best_id, , drop = FALSE])

  if (isTRUE(tune_verbose)) {
    message(
      "PS auto-tuning selected for ", ps_method, ": ",
      paste(names(best), unlist(best), sep = "=", collapse = ", "),
      "; validation log-loss = ", signif(loss[best_id], 4)
    )
  }

  best
}

tune_rsf_survival <- function(train_data,
                              xvars,
                              landmark_times,
                              tune_n = 100,
                              tune_seed = 20260427,
                              tune_verbose = FALSE,
                              fast_rsf = TRUE) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package {ranger} is required for RSF tuning.")
  }

  d <- as.data.frame(train_data)
  d0 <- d[d$A == 0, , drop = FALSE]

  if (nrow(d0) < 30 || sum(d0$delta == 1) < 8) {
    if (isTRUE(tune_verbose)) {
      message("RSF auto-tuning skipped: too few control rows or events.")
    }
    return(list())
  }

  d_tune <- head(d0, min(tune_n, nrow(d0)))

  if (nrow(d_tune) < 30 || sum(d_tune$delta == 1) < 8) {
    if (isTRUE(tune_verbose)) {
      message("RSF auto-tuning skipped: too few events in tuning subset.")
    }
    return(list())
  }

  set.seed(tune_seed)
  idx_all <- seq_len(nrow(d_tune))
  idx_event <- which(d_tune$delta == 1)
  idx_cens <- which(d_tune$delta == 0)

  train_event <- sample(
    idx_event,
    size = max(1, floor(0.50 * length(idx_event))),
    replace = FALSE
  )

  train_cens <- if (length(idx_cens) > 0) {
    sample(
      idx_cens,
      size = max(1, floor(0.50 * length(idx_cens))),
      replace = FALSE
    )
  } else {
    integer(0)
  }

  idx_train <- sort(c(train_event, train_cens))
  idx_valid <- setdiff(idx_all, idx_train)

  if (length(idx_valid) < 10 || sum(d_tune$delta[idx_valid] == 1) < 3) {
    if (isTRUE(tune_verbose)) {
      message("RSF auto-tuning skipped: too few validation events.")
    }
    return(list())
  }

  d_train <- d_tune[idx_train, , drop = FALSE]
  d_valid <- d_tune[idx_valid, , drop = FALSE]

  p_dim <- length(xvars)

  if (isTRUE(fast_rsf)) {
    ## RSF tuning grid.
    ## This is recommended for large Monte Carlo simulations.
    grid <- expand.grid(
      rsf_mtry = sort(unique(pmax(
        1,
        pmin(p_dim, c(floor(sqrt(p_dim)), ceiling(p_dim / 3), ceiling(p_dim / 2)))
      ))),
      rsf_min_node_size = c(15, 30, 50),
      rsf_num_trees = c(100),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    ## Original, more expensive tuning grid.
    grid <- expand.grid(
      rsf_mtry = sort(unique(pmax(
        1,
        pmin(p_dim, c(floor(sqrt(p_dim)), ceiling(p_dim / 2), p_dim))
      ))),
      rsf_min_node_size = c(5, 15, 30),
      rsf_num_trees = c(150, 300),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }

  score <- rep(-Inf, nrow(grid))
  t0 <- median(landmark_times)

  for (gg in seq_len(nrow(grid))) {
    pars <- as.list(grid[gg, , drop = FALSE])

    score[gg] <- tryCatch(
      {
        if (isTRUE(fast_rsf)) {
          fit <- ranger::ranger(
            formula = Surv(time, delta) ~ .,
            data = d_train[, c("time", "delta", xvars), drop = FALSE],
            num.trees = pars$rsf_num_trees,
            mtry = pars$rsf_mtry,
            min.node.size = pars$rsf_min_node_size,

            ## Extremely randomized survival splits.
            splitrule = "extratrees",
            num.random.splits = 1,

            ## Sampling without replacement.
            replace = FALSE,
            sample.fraction = 0.632,

            ## Retain trees for validation predictions.
            write.forest = TRUE,

            ## Avoid unnecessary OOB-error calculation during tuning.
            oob.error = FALSE,

            ## Keep 1 because the outer Monte Carlo is already parallelized.
            num.threads = 1,
            seed = tune_seed
          )
        } else {
          fit <- ranger::ranger(
            formula = Surv(time, delta) ~ .,
            data = d_train[, c("time", "delta", xvars), drop = FALSE],
            num.trees = pars$rsf_num_trees,
            mtry = pars$rsf_mtry,
            min.node.size = pars$rsf_min_node_size,
            #_____________________________
            # Alternative split rule:
            #splitrule = "logrank",
            # Extremely randomized splits:
            splitrule = "extratrees",
            num.random.splits = 1,
            replace = FALSE,
            sample.fraction = 0.632,
            #_____________________________
            write.forest = TRUE,
            num.threads = 1,
            seed = tune_seed
          )
        }

        pred <- predict(fit, data = d_valid[, xvars, drop = FALSE])
        surv_mat <- pred$survival
        times <- pred$unique.death.times

        if (is.null(times)) {
          times <- fit$unique.death.times
        }

        if (is.null(surv_mat) || is.null(times)) {
          stop("RSF prediction failed.")
        }

        ## Extremely randomized splits:er extraction of S(t0 | X) than row-wise approx().
        idx_t0 <- findInterval(t0, times)
        idx_t0 <- pmin(pmax(idx_t0, 1), length(times))

        S_t0 <- surv_mat[, idx_t0]

        risk <- 1 - as.numeric(S_t0)

        cc <- survival::concordance(
          survival::Surv(time, delta) ~ risk,
          data = d_valid
        )$concordance

        ## Orientation-free C-index.
        max(cc, 1 - cc)
      },
      error = function(e) -Inf
    )
  }

  if (all(!is.finite(score))) {
    if (isTRUE(tune_verbose)) {
      message("RSF auto-tuning failed for all candidates; using defaults.")
    }
    return(list())
  }

  best_id <- which.max(ifelse(is.finite(score), score, -Inf))
  best_raw <- as.list(grid[best_id, , drop = FALSE])

  best <- list(
    num.trees = best_raw$rsf_num_trees,
    mtry = best_raw$rsf_mtry,
    min.node.size = best_raw$rsf_min_node_size
  )

  if (isTRUE(tune_verbose)) {
    message(
      "RSF auto-tuning selected: ",
      paste(names(best), unlist(best), sep = "=", collapse = ", "),
      "; validation C-index = ", signif(score[best_id], 4),
      "; fast_rsf = ", fast_rsf
    )
  }

  best
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
  ## Tune machine-learning methods; GLM uses maximum likelihood.
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
    } else if (!is.null(fit$yhat.test.mean)) {
      raw <- as.numeric(fit$yhat.test.mean)

      if (all(raw >= -0.05 & raw <= 1.05, na.rm = TRUE)) {
        ehat <- raw
      } else {
        ehat <- pnorm(raw + binary_offset)
      }
    } else if (!is.null(fit$yhat.test)) {
      raw <- as.numeric(colMeans(fit$yhat.test))

      if (all(raw >= -0.05 & raw <= 1.05, na.rm = TRUE)) {
        ehat <- raw
      } else {
        ehat <- pnorm(raw + binary_offset)
      }
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

fit_att_ipw_linwei_ml <- function(data,
                                  xvars,
                                  ps_method = c("glm", "nnet", "rf", "bart"),
                                  trim = c(0.01, 0.99),
                                  normalize_controls = TRUE,
                                  auto_tune = FALSE,
                                  tune_n = 100,
                                  tune_seed = 20260427,
                                  tune_verbose = FALSE) {
  ps_method <- match.arg(ps_method)

  d <- as.data.frame(data)

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

  ehat <- ps_obj$ehat
  qhat <- ehat / (1 - ehat)

  if (normalize_controls) {
    n1 <- sum(d$A == 1)
    qhat[d$A == 0] <- n1 * qhat[d$A == 0] / sum(qhat[d$A == 0])
  }

  W <- with(d, A + (1 - A) * qhat)

  fit <- coxph(
    Surv(time, delta) ~ A,
    data = d,
    weights = W,
    robust = TRUE,
    ties = "breslow"
  )

  theta_hat <- unname(coef(fit)["A"])

  ## Fixed-weight Lin-Wei/Binder variance.
  se_lw <- manual_linwei_se(
    theta_hat = theta_hat,
    data = d,
    W = W
  )

  idx0 <- which(d$A == 0)
  wd <- .weight_diag(W[idx0])

  list(
    theta = theta_hat,
    se_lw = se_lw,
    W = W,
    ehat = ehat,
    qhat = qhat,
    ps_method = ps_method,
    ESS = wd$ESS,
    Rel_ESS = wd$Rel_ESS,
    W_CV = wd$W_CV,
    W_Min = wd$W_Min,
    W_Max = wd$W_Max
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
    # Alternative split rule:
    #splitrule = "logrank",
    # Extremely randomized splits:
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

compute_true_att_target <- function(n1_super = 30000,
                                    n0_super = 60000,
                                    p = 10,
                                    ps_truth = "linear",
                                    outcome_truth = "linear_ph",
                                    ps_nonlinearity = 1,
                                    or_nonlinearity = 1,
                                    beta_cond = log(0.70),
                                    xi = 0,
                                    h_type = "constant",
                                    lambda0 = 0.00008,
                                    eta_shape = 2.0,
                                    censor_rate = 0.0008,
                                    seed = 12345) {
  set.seed(seed)

  dat <- generate_observed_data(
    n1 = n1_super,
    n0 = n0_super,
    p = p,
    ps_truth = ps_truth,
    outcome_truth = outcome_truth,
    ps_nonlinearity = ps_nonlinearity,
    or_nonlinearity = or_nonlinearity,
    beta_cond = beta_cond,
    xi = xi,
    h_type = h_type,
    lambda0 = lambda0,
    eta_shape = eta_shape,
    censor_rate = censor_rate
  )

  xvars <- paste0("X", seq_len(p))
  X <- as.matrix(dat[, xvars, drop = FALSE])

  e_true <- source_score(
    X,
    ps_truth = ps_truth,
    ps_nonlinearity = ps_nonlinearity
  )
  q_true <- e_true / (1 - e_true)

  idx0 <- which(dat$A == 0)
  idx1 <- which(dat$A == 1)

  q_true[idx0] <- length(idx1) * q_true[idx0] / sum(q_true[idx0])

  W_true <- with(dat, A + (1 - A) * q_true)

  fit <- coxph(
    Surv(time, delta) ~ A,
    data = dat,
    weights = W_true,
    ties = "breslow"
  )

  unname(coef(fit)["A"])
}

empty_row <- function(rep_id,
                      n1,
                      n0,
                      method,
                      theta_true,
                      error_msg = NA_character_) {
  data.frame(
    rep = rep_id,
    n1 = n1,
    n0 = n0,
    Method = method,
    theta_true = theta_true,
    Estimate = NA_real_,
    SE = NA_real_,
    CI_L = NA_real_,
    CI_U = NA_real_,
    ESS = NA_real_,
    Rel_ESS = NA_real_,
    W_CV = NA_real_,
    W_Min = NA_real_,
    W_Max = NA_real_,
    Cal_Grad = NA_real_,
    Cal_Converged = NA,
    Error = error_msg,
    stringsAsFactors = FALSE
  )
}

fit_one_rep <- function(rep_id,
                        n1,
                        n0,
                        theta_true,
                        p = 10,
                        ps_truth = "linear",
                        outcome_truth = "linear_ph",
                        ps_nonlinearity = 1,
                        or_nonlinearity = 1,
                        beta_cond = log(0.70),
                        xi = 0,
                        h_type = "constant",
                        lambda0 = 0.00008,
                        eta_shape = 2.0,
                        censor_rate = 0.0008,
                        mec_ps_method = "glm",
                        mec_or_method = "cox",
                        Kfold_mec = 5,
                        n_landmarks = 4,
                        landmark_range = c(0.25, 0.70),
                        landmark_probs = NULL,
                        landmark_times = NULL,
                        mec_divergence = "kl",
                        trim = c(0.01, 0.99),
                        include_mec_no_cf = FALSE,
                        include_unweighted_cox = TRUE,
                        include_linwei_ml = TRUE,
                        Lin_Wei_ps_method = "glm",
                        auto_tune_ml = FALSE,
                        auto_tune_n = 100,
                        auto_tune_seed = 20260427,
                        auto_tune_verbose = FALSE) {
  dat <- generate_observed_data(
    n1 = n1,
    n0 = n0,
    p = p,
    ps_truth = ps_truth,
    outcome_truth = outcome_truth,
    ps_nonlinearity = ps_nonlinearity,
    or_nonlinearity = or_nonlinearity,
    beta_cond = beta_cond,
    xi = xi,
    h_type = h_type,
    lambda0 = lambda0,
    eta_shape = eta_shape,
    censor_rate = censor_rate
  )

  xvars <- paste0("X", seq_len(p))
  rows <- list()

  ## 1. Optional Unweighted Cox benchmark
  if (isTRUE(include_unweighted_cox)) {
    ans_unweighted <- tryCatch(
      fit_unweighted_cox(dat),
      error = function(e) e
    )

    if (inherits(ans_unweighted, "error")) {
      rows[[length(rows) + 1]] <- empty_row(
        rep_id, n1, n0, "Unweighted Cox", theta_true, conditionMessage(ans_unweighted)
      )
    } else {
      rows[[length(rows) + 1]] <- data.frame(
        rep = rep_id,
        n1 = n1,
        n0 = n0,
        Method = "Unweighted Cox",
        theta_true = theta_true,
        Estimate = ans_unweighted$theta,
        SE = ans_unweighted$se,
        CI_L = ans_unweighted$theta - 1.96 * ans_unweighted$se,
        CI_U = ans_unweighted$theta + 1.96 * ans_unweighted$se,
        ESS = NA_real_,
        Rel_ESS = NA_real_,
        W_CV = NA_real_,
        W_Min = NA_real_,
        W_Max = NA_real_,
        Cal_Grad = NA_real_,
        Cal_Converged = NA,
        Error = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  ## 2. Standard ATT-IPW Cox using logistic GLM PS:
  ##    naive model-based SE, Lin-Wei/Binder SE, and Shu corrected SE.
  ans_ipw <- tryCatch(
    fit_att_ipw_naive_linwei_shu(
      data = dat,
      xvars = xvars,
      trim = trim,
      normalize_controls = TRUE
    ),
    error = function(e) e
  )

  if (inherits(ans_ipw, "error")) {
    rows[[length(rows) + 1]] <- empty_row(
      rep_id, n1, n0, "ATT-IPW Cox: Naive model-based", theta_true, conditionMessage(ans_ipw)
    )
    rows[[length(rows) + 1]] <- empty_row(
      rep_id, n1, n0, "ATT-IPW Cox: Lin-Wei", theta_true, conditionMessage(ans_ipw)
    )
    rows[[length(rows) + 1]] <- empty_row(
      rep_id, n1, n0, "ATT-IPW Cox: Shu", theta_true, conditionMessage(ans_ipw)
    )
  } else {
    rows[[length(rows) + 1]] <- data.frame(
      rep = rep_id,
      n1 = n1,
      n0 = n0,
      Method = "ATT-IPW Cox: Naive model-based",
      theta_true = theta_true,
      Estimate = ans_ipw$theta,
      SE = ans_ipw$se_naive,
      CI_L = ans_ipw$theta - 1.96 * ans_ipw$se_naive,
      CI_U = ans_ipw$theta + 1.96 * ans_ipw$se_naive,
      ESS = NA_real_,
      Rel_ESS = NA_real_,
      W_CV = NA_real_,
      W_Min = NA_real_,
      W_Max = NA_real_,
      Cal_Grad = NA_real_,
      Cal_Converged = NA,
      Error = NA_character_,
      stringsAsFactors = FALSE
    )

    rows[[length(rows) + 1]] <- data.frame(
      rep = rep_id,
      n1 = n1,
      n0 = n0,
      Method = "ATT-IPW Cox: Lin-Wei",
      theta_true = theta_true,
      Estimate = ans_ipw$theta,
      SE = ans_ipw$se_lw,
      CI_L = ans_ipw$theta - 1.96 * ans_ipw$se_lw,
      CI_U = ans_ipw$theta + 1.96 * ans_ipw$se_lw,
      ESS = NA_real_,
      Rel_ESS = NA_real_,
      W_CV = NA_real_,
      W_Min = NA_real_,
      W_Max = NA_real_,
      Cal_Grad = NA_real_,
      Cal_Converged = NA,
      Error = NA_character_,
      stringsAsFactors = FALSE
    )

    rows[[length(rows) + 1]] <- data.frame(
      rep = rep_id,
      n1 = n1,
      n0 = n0,
      Method = "ATT-IPW Cox: Shu",
      theta_true = theta_true,
      Estimate = ans_ipw$theta,
      SE = ans_ipw$se_shu,
      CI_L = ans_ipw$theta - 1.96 * ans_ipw$se_shu,
      CI_U = ans_ipw$theta + 1.96 * ans_ipw$se_shu,
      ESS = NA_real_,
      Rel_ESS = NA_real_,
      W_CV = NA_real_,
      W_Min = NA_real_,
      W_Max = NA_real_,
      Cal_Grad = NA_real_,
      Cal_Converged = NA,
      Error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  ## 3. Additional Lin-Wei with ML-estimated PS weights, no CF
  if (isTRUE(include_linwei_ml)) {
    linwei_ml_label <- paste0(
      "ATT-IPW Cox: Lin-Wei (PS=", Lin_Wei_ps_method, ", no CF)"
    )

    ans_lw_ml <- tryCatch(
      fit_att_ipw_linwei_ml(
        data = dat,
        xvars = xvars,
        ps_method = Lin_Wei_ps_method,
        trim = trim,
        normalize_controls = TRUE,
        auto_tune = auto_tune_ml,
        tune_n = auto_tune_n,
        tune_seed = auto_tune_seed + rep_id,
        tune_verbose = auto_tune_verbose
      ),
      error = function(e) e
    )

    if (inherits(ans_lw_ml, "error")) {
      rows[[length(rows) + 1]] <- empty_row(
        rep_id, n1, n0, linwei_ml_label, theta_true, conditionMessage(ans_lw_ml)
      )
    } else {
      rows[[length(rows) + 1]] <- data.frame(
        rep = rep_id,
        n1 = n1,
        n0 = n0,
        Method = linwei_ml_label,
        theta_true = theta_true,
        Estimate = ans_lw_ml$theta,
        SE = ans_lw_ml$se_lw,
        CI_L = ans_lw_ml$theta - 1.96 * ans_lw_ml$se_lw,
        CI_U = ans_lw_ml$theta + 1.96 * ans_lw_ml$se_lw,
        ## Weight-stability diagnostics are reported only for MEC-Cox.
        ## Keep these as NA for non-MEC methods to avoid mixing diagnostics
        ## from ordinary IPW weights with calibrated MEC weights.
        ESS = NA_real_,
        Rel_ESS = NA_real_,
        W_CV = NA_real_,
        W_Min = NA_real_,
        W_Max = NA_real_,
        Cal_Grad = NA_real_,
        Cal_Converged = NA,
        Error = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  ## 4. MEC-Cox with cross-fitting by default
  mec_label_cf <- paste0(
    "MEC-Cox: PS=", mec_ps_method,
    ", OR=", mec_or_method,
    ", L=", n_landmarks,
    ", K=", Kfold_mec,
    ", ", mec_divergence
  )

  ans_mec_cf <- tryCatch(
    fit_mec_cox(
      data = dat,
      xvars = xvars,
      ps_method = mec_ps_method,
      or_method = mec_or_method,
      use_crossfit = TRUE,
      Kfold = Kfold_mec,
      n_landmarks = n_landmarks,
      landmark_range = landmark_range,
      landmark_probs = landmark_probs,
      landmark_times = landmark_times,
      divergence = mec_divergence,
      trim = trim,
      auto_tune = auto_tune_ml,
      tune_n = auto_tune_n,
      tune_seed = auto_tune_seed + rep_id,
      tune_verbose = auto_tune_verbose
    ),
    error = function(e) e
  )

  if (inherits(ans_mec_cf, "error")) {
    rows[[length(rows) + 1]] <- empty_row(
      rep_id, n1, n0, mec_label_cf, theta_true, conditionMessage(ans_mec_cf)
    )
  } else {
    mec_label_cf <- paste0(
      "MEC-Cox: PS=", mec_ps_method,
      ", OR=", mec_or_method,
      ", L=", ncol(ans_mec_cf$H) - 1,
      ", K=", Kfold_mec,
      ", ", mec_divergence
    )

    rows[[length(rows) + 1]] <- data.frame(
      rep = rep_id,
      n1 = n1,
      n0 = n0,
      Method = mec_label_cf,
      theta_true = theta_true,
      Estimate = ans_mec_cf$theta,
      SE = ans_mec_cf$se,
      CI_L = ans_mec_cf$theta - 1.96 * ans_mec_cf$se,
      CI_U = ans_mec_cf$theta + 1.96 * ans_mec_cf$se,
      ESS = ans_mec_cf$ESS,
      Rel_ESS = ans_mec_cf$Rel_ESS,
      W_CV = ans_mec_cf$W_CV,
      W_Min = ans_mec_cf$W_Min,
      W_Max = ans_mec_cf$W_Max,
      Cal_Grad = ans_mec_cf$calibration_grad_norm,
      Cal_Converged = ans_mec_cf$calibration_converged,
      Error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  ## 5. Optional MEC-Cox without cross-fitting
  if (isTRUE(include_mec_no_cf)) {
    mec_label_nocf <- paste0(
      "MEC-Cox (no CF): PS=", mec_ps_method,
      ", OR=", mec_or_method,
      ", L=", n_landmarks,
      ", ", mec_divergence
    )

    ans_mec_nocf <- tryCatch(
      fit_mec_cox(
        data = dat,
        xvars = xvars,
        ps_method = mec_ps_method,
        or_method = mec_or_method,
        use_crossfit = FALSE,
        Kfold = Kfold_mec,
        n_landmarks = n_landmarks,
        landmark_range = landmark_range,
        landmark_probs = landmark_probs,
        landmark_times = landmark_times,
        divergence = mec_divergence,
        trim = trim,
        auto_tune = auto_tune_ml,
        tune_n = auto_tune_n,
        tune_seed = auto_tune_seed + rep_id,
        tune_verbose = auto_tune_verbose
      ),
      error = function(e) e
    )

    if (inherits(ans_mec_nocf, "error")) {
      rows[[length(rows) + 1]] <- empty_row(
        rep_id, n1, n0, mec_label_nocf, theta_true, conditionMessage(ans_mec_nocf)
      )
    } else {
      mec_label_nocf <- paste0(
        "MEC-Cox (no CF): PS=", mec_ps_method,
        ", OR=", mec_or_method,
        ", L=", ncol(ans_mec_nocf$H) - 1,
        ", ", mec_divergence
      )

      rows[[length(rows) + 1]] <- data.frame(
        rep = rep_id,
        n1 = n1,
        n0 = n0,
        Method = mec_label_nocf,
        theta_true = theta_true,
        Estimate = ans_mec_nocf$theta,
        SE = ans_mec_nocf$se,
        CI_L = ans_mec_nocf$theta - 1.96 * ans_mec_nocf$se,
        CI_U = ans_mec_nocf$theta + 1.96 * ans_mec_nocf$se,
        ESS = ans_mec_nocf$ESS,
        Rel_ESS = ans_mec_nocf$Rel_ESS,
        W_CV = ans_mec_nocf$W_CV,
        W_Min = ans_mec_nocf$W_Min,
        W_Max = ans_mec_nocf$W_Max,
        Cal_Grad = ans_mec_nocf$calibration_grad_norm,
        Cal_Converged = ans_mec_nocf$calibration_converged,
        Error = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

required_worker_packages <- function(mec_ps_method,
                                     mec_or_method,
                                     Lin_Wei_ps_method = "glm") {
  pkgs <- c("survival", "MASS", "dplyr")

  ps_methods_needed <- unique(c(mec_ps_method, Lin_Wei_ps_method))

  if ("bart" %in% ps_methods_needed) pkgs <- c(pkgs, "dbarts")
  if ("nnet" %in% ps_methods_needed) pkgs <- c(pkgs, "nnet")
  if ("rf" %in% ps_methods_needed) pkgs <- c(pkgs, "ranger")

  if (mec_or_method == "rsf") pkgs <- c(pkgs, "ranger")

  unique(pkgs)
}
