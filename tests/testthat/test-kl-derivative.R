kl_derivative_test_data <- function() {
  set.seed(6192)
  source <- c(rep(1L, 80L), rep(0L, 160L))
  x <- rnorm(length(source), mean = 0.25 * source)
  z <- rnorm(length(source), mean = 0.15 * source)
  event_time <- rexp(length(source), 0.08 * exp(-0.3 * source + 0.4 * x))
  censor_time <- rexp(length(source), 0.025)
  data.frame(
    time = pmax(0.25, round(pmin(event_time, censor_time) * 2) / 2),
    event = as.integer(event_time <= censor_time), source = source, x = x, z = z
  )
}

kl_variance_arguments <- function(fit) {
  controls <- which(fit$analysis_data$A == 0L)
  targets <- which(fit$analysis_data$A == 1L)
  control_basis <- fit$basis[controls, , drop = FALSE]
  list(
    theta = fit$theta, data = fit$analysis_data, weights = fit$weights,
    base_weights = fit$baseline_weights[controls],
    control_basis = control_basis,
    target_basis = fit$basis[targets, , drop = FALSE],
    control_rows = controls, target_rows = targets,
    lambda = fit$calibration$lambda,
    calibration_hessian = crossprod(control_basis * fit$weights[controls],
                                   control_basis)
  )
}

# Differentiate each tied event set and its complete risk set explicitly.
kl_explicit_breslow_variance <- function(arguments) {
  d <- arguments$data
  w <- arguments$weights
  theta <- arguments$theta
  n <- nrow(d)
  p <- ncol(arguments$control_basis)
  controls <- arguments$control_rows
  targets <- arguments$target_rows
  derivative_weights <- matrix(0, n, p)
  derivative_weights[controls, ] <- arguments$control_basis * w[controls]
  rho <- derivative_weights
  rho[targets, ] <- -arguments$target_basis
  event_times <- sort(unique(d$time[d$delta == 1L]))

  score <- function(coefficient) {
    total <- 0
    for (time in event_times) {
      events <- which(d$time == time & d$delta == 1L)
      risk <- which(d$time >= time)
      weighted_risk <- w[risk] * exp(coefficient * d$A[risk])
      mean_source <- sum(weighted_risk * d$A[risk]) / sum(weighted_risk)
      total <- total + sum(w[events] * (d$A[events] - mean_source))
    }
    total
  }

  eta <- numeric(n)
  b <- numeric(p)
  for (time in event_times) {
    events <- which(d$time == time & d$delta == 1L)
    risk <- which(d$time >= time)
    relative_risk <- exp(theta * d$A[risk])
    s0 <- sum(w[risk] * relative_risk)
    mean_source <- sum(w[risk] * relative_risk * d$A[risk]) / s0
    event_mass <- sum(w[events])
    ds0 <- colSums(derivative_weights[risk, , drop = FALSE] * relative_risk)
    ds1 <- colSums(derivative_weights[risk, , drop = FALSE] *
                    (relative_risk * d$A[risk]))
    b <- b + colSums(derivative_weights[events, , drop = FALSE] *
                      (d$A[events] - mean_source)) -
      event_mass * (ds1 - mean_source * ds0) / s0
    eta[events] <- eta[events] + w[events] * (d$A[events] - mean_source)
    eta[risk] <- eta[risk] - w[risk] * relative_risk * event_mass / s0 *
      (d$A[risk] - mean_source)
  }

  # Retain the existing scalar numerical derivative and SVD convention.
  step <- 1e-5 * max(1, abs(theta))
  a <- (score(theta + step) - score(theta - step)) / (2 * step)
  jacobian <- rbind(c(a, b), cbind(0, arguments$calibration_hessian))
  decomposition <- svd(jacobian)
  retained <- decomposition$d > max(decomposition$d) * 1e-8
  inverse <- decomposition$v[, retained, drop = FALSE] %*%
    (t(decomposition$u[, retained, drop = FALSE]) / decomposition$d[retained])
  influence <- -(cbind(eta, rho) %*% t(inverse))[, 1L]
  list(se = sqrt(sum(influence^2)), influence = influence,
       rank = sum(retained), used_pseudoinverse = !all(retained))
}

test_that("KL derivative choice preserves fitting and the numerical default", {
  sample_data <- kl_derivative_test_data()
  expect_true(anyDuplicated(sample_data$time[sample_data$event == 1L]) > 0L)
  arguments <- list(
    data = sample_data, time = "time", event = "event", source = "source",
    covariates = c("x", "z"), ps_learner = "glm", survival_learner = "cox",
    n_folds = 3L, n_landmarks = 1L, seed = 715L
  )
  previous_rng <- .Random.seed
  default <- do.call(fit_mec_cox, arguments)
  expect_identical(.Random.seed, previous_rng)
  numerical <- do.call(fit_mec_cox, c(arguments, list(kl_derivative = "numerical")))
  expect_identical(.Random.seed, previous_rng)
  analytic <- (function() {
    scalar_difference <- mecCox:::.mec_central_difference
    local_mocked_bindings(
      .mec_central_difference = function(fun, x, step = 1e-5) {
        if (length(x) != 1L) {
          stop("Analytic KL must not perturb the calibration coefficients.")
        }
        scalar_difference(fun, x, step)
      },
      .package = "mecCox"
    )
    do.call(fit_mec_cox, c(arguments, list(kl_derivative = "analytic")))
  })()
  expect_identical(.Random.seed, previous_rng)

  expect_identical(default$kl_derivative, "numerical")
  expect_identical(numerical$kl_derivative, "numerical")
  expect_identical(analytic$kl_derivative, "analytic")
  expect_identical(analytic$variance_method, "stacked_working")
  expect_gt(ncol(analytic$basis), 1L)
  fitted_fields <- c("theta", "hr", "weights", "baseline_weights", "ps_oof",
                     "survival_oof", "basis", "fold_id", "landmark_times",
                     "propensity_tuning", "calibration")
  for (field in fitted_fields) {
    expect_identical(default[[field]], numerical[[field]], info = field)
    expect_identical(default[[field]], analytic[[field]], info = field)
  }
  for (field in c("se", "conf_int", "variance_diagnostics")) {
    expect_identical(default[[field]], numerical[[field]], info = field)
  }

  variance_arguments <- kl_variance_arguments(default)
  default_variance <- do.call(mecCox:::.mec_stacked_variance, variance_arguments)
  numerical_variance <- do.call(mecCox:::.mec_stacked_variance,
                               c(variance_arguments, list(kl_derivative = "numerical")))
  expect_identical(default_variance, numerical_variance)
  expect_identical(default_variance$standard_error, default$se)
  expect_identical(.Random.seed, previous_rng)

  analytic_variance <- do.call(mecCox:::.mec_stacked_variance,
                              c(variance_arguments, list(kl_derivative = "analytic")))
  reference <- kl_explicit_breslow_variance(variance_arguments)
  expect_identical(analytic_variance$standard_error, analytic$se)
  expect_equal(analytic$se, reference$se, tolerance = 1e-8)
  expect_equal(analytic_variance$influence, reference$influence, tolerance = 1e-8)
  expect_equal(analytic_variance$jacobian_rank, reference$rank)
  expect_identical(analytic_variance$used_pseudoinverse, reference$used_pseudoinverse)
  expected_interval <- exp(analytic$theta + c(-1, 1) *
                             qnorm((1 + analytic$conf_level) / 2) * reference$se)
  expect_equal(unname(analytic$conf_int), expected_interval, tolerance = 1e-8)
  expect_identical(.Random.seed, previous_rng)
})

test_that("analytic KL derivatives retain the existing generalized inverse", {
  d <- data.frame(
    A = rep(c(1L, 0L), 6), time = rep(1:6, each = 2),
    delta = c(1L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L)
  )
  x <- seq(-1, 1, length.out = nrow(d))
  basis <- cbind(1, x, 2 * x)
  controls <- which(d$A == 0L)
  targets <- which(d$A == 1L)
  weights <- rep(1, nrow(d))
  weights[controls] <- c(0.7, 1.2, 0.9, 1.1, 0.8, 1.3)
  control_basis <- basis[controls, , drop = FALSE]
  arguments <- list(
    theta = -0.25, data = d, weights = weights,
    base_weights = weights[controls], control_basis = control_basis,
    target_basis = basis[targets, , drop = FALSE], control_rows = controls,
    target_rows = targets, lambda = numeric(ncol(basis)),
    calibration_hessian = crossprod(control_basis * weights[controls], control_basis)
  )
  expect_warning(
    analytic <- do.call(mecCox:::.mec_stacked_variance,
                        c(arguments, list(kl_derivative = "analytic"))),
    "rank deficient; using its generalized inverse"
  )
  expect_warning(
    numerical <- do.call(mecCox:::.mec_stacked_variance,
                         c(arguments, list(kl_derivative = "numerical"))),
    "rank deficient; using its generalized inverse"
  )
  reference <- kl_explicit_breslow_variance(arguments)
  expect_true(analytic$used_pseudoinverse)
  expect_identical(analytic$used_pseudoinverse, numerical$used_pseudoinverse)
  expect_equal(analytic$jacobian_rank, numerical$jacobian_rank)
  expect_equal(analytic$jacobian_rank, reference$rank)
  expect_lt(analytic$jacobian_rank, ncol(basis) + 1L)
  expect_equal(analytic$standard_error, reference$se, tolerance = 1e-9)
  expect_equal(analytic$influence, reference$influence, tolerance = 1e-9)
})

test_that("invalid KL derivative choices are rejected before data processing", {
  local_mocked_bindings(
    .pkg_validate_data = function(...) stop("Data validation must not run"),
    .package = "mecCox"
  )
  for (choice in c("unsupported", "ANALYTIC")) {
    expect_error(
      fit_mec_cox(NULL, "time", "event", "source", "x", kl_derivative = choice),
      "numerical.*analytic"
    )
  }
})
