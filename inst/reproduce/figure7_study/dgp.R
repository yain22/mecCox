# Data generation for the original MEC-Cox workflow.
# Five binary and five continuous covariates; no nuisance fitting is performed here.

dgp_gauss_legendre <- function(n = 32L) {
  stopifnot(length(n) == 1L, n >= 2L)
  j <- seq_len(n - 1L)
  off <- j / sqrt(4 * j^2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(j, j + 1L)] <- off
  jacobi[cbind(j + 1L, j)] <- off
  e <- eigen(jacobi, symmetric = TRUE)
  order_index <- order(e$values)
  list(nodes = e$values[order_index],
       weights = 2 * e$vectors[1L, order_index]^2)
}

dgp_irwin_hall5_density <- function(v) {
  # Reflection avoids cancellation in the upper tail of the polynomial density.
  x <- pmin(v, 5 - v)
  ans <- numeric(length(v))
  inside <- x > 0 & x <= 2.5
  for (j in 0:2) {
    ans[inside] <- ans[inside] + (-1)^j * choose(5, j) *
      pmax(x[inside] - j, 0)^4 / factorial(4)
  }
  ans
}

dgp_target_quadrature <- function(beta, nodes_per_interval = 32L) {
  stopifnot(length(beta) == 10L,
            max(abs(beta[6:10] - beta[6L])) < 1e-14)
  gl <- dgp_gauss_legendre(nodes_per_interval)
  v <- unlist(lapply(0:4, function(j) j + (gl$nodes + 1) / 2),
              use.names = FALSE)
  v_weights <- rep(gl$weights / 2, 5L) * dgp_irwin_hall5_density(v)
  binary <- as.matrix(expand.grid(rep(list(c(-1, 1)), 5L)))
  eta_binary <- drop(binary %*% beta[1:5])
  eta_continuous <- beta[6L] * (2 * v - 5)
  rates <- exp(rep(eta_binary, each = length(v)) +
                 rep(eta_continuous, times = nrow(binary)))
  weights <- rep(v_weights / nrow(binary), times = nrow(binary))
  stopifnot(all(weights > 0), abs(sum(weights) - 1) < 1e-12)
  list(rates = rates, weights = weights)
}

dgp_marginal_survival_direct <- function(time, quadrature) {
  stopifnot(all(is.finite(time)), all(time >= 0))
  answer <- numeric(length(time))
  # Small batches keep precomputation inexpensive in memory.
  groups <- split(seq_along(time), ceiling(seq_along(time) / 128L))
  for (index in groups) {
    answer[index] <- drop(exp(-outer(time[index], quadrature$rates)) %*%
                            quadrature$weights)
  }
  pmin(answer, 1)
}

dgp_make_transform <- function(beta, tau, nodes_per_interval = 32L,
                               grid_size = 4097L) {
  quadrature <- dgp_target_quadrature(beta, nodes_per_interval)
  time <- seq(0, tau, length.out = grid_size)
  survival <- dgp_marginal_survival_direct(time, quadrature)
  H <- -log(survival)
  H[1L] <- 0
  stopifnot(all(diff(H) > 0), all(is.finite(H)))
  list(time = time, H = H, quadrature = quadrature,
       H_function = splinefun(time, H, method = "monoH.FC"),
       inverse_function = splinefun(H, time, method = "monoH.FC"),
       nodes_per_interval = nodes_per_interval, grid_size = grid_size)
}

.dgp_transform_cache <- new.env(parent = emptyenv())

dgp_config <- function(hr = 0.70, target_cohort_probability = 1/3) {
  stopifnot(length(hr) == 1L, hr %in% c(1, 0.70),
            target_cohort_probability > 0, target_cohort_probability < 1)
  beta <- c(0.45, 0.35, 0.25, 0.20, 0.15, rep(0.20, 5L))
  tau <- log(8)
  cache_key <- paste(c(beta, tau), collapse = ":")
  if (!exists(cache_key, envir = .dgp_transform_cache, inherits = FALSE)) {
    assign(cache_key, dgp_make_transform(beta, tau, grid_size = 12289L), envir = .dgp_transform_cache)
  }
  list(hr = hr, beta = beta,
       p_control = c(0.40, 0.42, 0.44, 0.46, 0.48),
       p_target = rep(0.5, 5L),
       continuous_tilt = c(0.30, 0.25, 0.20, 0.15, 0.10),
       censor_rate = c(0.30, 0.50), tau = tau,
       target_cohort_probability = target_cohort_probability,
       covariate_names = paste0("X", 1:10),
       transform = get(cache_key, envir = .dgp_transform_cache, inherits = FALSE))
}

dgp_covariates <- function(data_or_X) {
  wanted <- paste0("X", 1:10)
  if (!is.null(colnames(data_or_X)) && all(wanted %in% colnames(data_or_X))) {
    X <- as.matrix(data_or_X[, wanted, drop = FALSE])
  } else {
    X <- as.matrix(data_or_X)
  }
  stopifnot(ncol(X) == 10L, all(is.finite(X)), all(abs(X) <= 1))
  X
}

dgp_true_log_odds <- function(data_or_X, config) {
  X <- dgp_covariates(data_or_X)
  p <- config$p_control
  k <- config$continuous_tilt
  intercept <- sum(-log(2) - 0.5 * log(p * (1 - p))) +
    sum(log(sinh(k) / k)) + qlogis(config$target_cohort_probability)
  drop(intercept + X %*% c(0.5 * log((1 - p) / p), k))
}

dgp_true_ps <- function(data_or_X, config) {
  plogis(dgp_true_log_odds(data_or_X, config))
}

dgp_true_survival <- function(data_or_X, times, config) {
  X <- dgp_covariates(data_or_X)
  stopifnot(all(is.finite(times)), all(times >= 0))
  exp(-outer(exp(drop(X %*% config$beta)), times))
}

dgp_target_survival <- function(times, config) {
  dgp_marginal_survival_direct(times, config$transform$quadrature)
}

dgp_draw_arm <- function(n, arm, config) {
  stopifnot(length(n) == 1L, n >= 2L, arm %in% c(0L, 1L))
  p <- if (arm == 0L) config$p_control else config$p_target
  binary <- vapply(p, function(prob) 2 * rbinom(n, 1, prob) - 1, numeric(n))
  continuous <- vapply(config$continuous_tilt, function(k) {
    v <- runif(n)
    if (arm == 1L) 2 * v - 1 else -log(exp(k) - 2 * v * sinh(k)) / k
  }, numeric(n))
  X <- cbind(binary, continuous)
  colnames(X) <- config$covariate_names
  T0 <- rexp(n, rate = exp(drop(X %*% config$beta)))
  if (arm == 0L || config$hr == 1) {
    event <- T0
  } else {
    # H0 is the target marginal control cumulative hazard. The full potential
    # event time is T1 = H0^{-1}{H0(T0)/hr}. Events beyond tau need not be
    # transformed, because they cannot affect administratively censored data.
    H_tau <- tail(config$transform$H, 1L)
    cutoff <- config$transform$inverse_function(config$hr * H_tau)
    before_tau <- T0 <= cutoff
    event <- rep(Inf, n)
    if (any(before_tau)) {
      transformed_H <- config$transform$H_function(T0[before_tau]) / config$hr
      transformed_H <- pmin(pmax(transformed_H, 0), H_tau)
      event[before_tau] <- config$transform$inverse_function(transformed_H)
    }
  }
  censor <- rexp(n, rate = config$censor_rate[arm + 1L])
  data.frame(A = arm, X, time = pmin(event, censor, config$tau),
             delta = as.integer(event <= censor & event <= config$tau),
             row.names = NULL)
}

dgp_draw_data <- function(n0, n1, config) {
  rbind(dgp_draw_arm(n0, 0L, config), dgp_draw_arm(n1, 1L, config))
}
