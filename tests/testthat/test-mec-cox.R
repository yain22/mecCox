test_that("MEC-Cox balances cross-fitted control prognosis", {
  set.seed(20261008)
  n_target <- 85L
  n_control <- 170L
  n <- n_target + n_control
  source <- c(rep(1L, n_target), rep(0L, n_control))
  age <- rnorm(n, mean = 0.5 * source)
  marker <- rbinom(n, 1L, plogis(-0.2 + 0.3 * source))
  event_time <- rexp(n, rate = 0.07 * exp(0.4 * age + 0.5 * marker))
  censor_time <- rexp(n, rate = 0.025)
  observed_time <- pmin(event_time, censor_time, 24)
  status <- as.integer(event_time <= pmin(censor_time, 24))
  sample_data <- data.frame(
    months = observed_time, died = status,
    cohort = source, age = age, marker = marker
  )
  previous_random_state <- .Random.seed

  fit <- fit_mec_cox(
    data = sample_data, time = "months", event = "died",
    source = "cohort", covariates = c("age", "marker"),
    ps_learner = "glm", survival_learner = "cox",
    n_folds = 3L, n_landmarks = 2L, seed = 42L
  )
  expect_identical(.Random.seed, previous_random_state)

  expect_s3_class(fit, "mec_cox_fit")
  expect_equal(length(fit$weights), n)
  expect_true(all(fit$weights > 0))
  expect_equal(fit$weights[source == 1L], rep(1, n_target))
  expect_equal(sum(fit$weights[source == 0L]), n_target,
               tolerance = 1e-5)
  expect_equal(sum(fit$base_weights[source == 0L]), n_target,
               tolerance = 1e-9)
  expect_equal(length(fit$ps_oof), n)
  expect_equal(nrow(fit$basis), n)
  expect_equal(length(unique(fit$fold_id)), 3L)
  expect_true(is.finite(fit$theta) && is.finite(fit$se))
  expect_true(fit$se > 0)
  expect_true(fit$calibration$converged)

  target_total <- colSums(fit$basis[source == 1L, , drop = FALSE])
  weighted_total <- colSums(
    fit$basis[source == 0L, , drop = FALSE] *
      fit$weights[source == 0L]
  )
  expect_equal(weighted_total, target_total, tolerance = 1e-5)

  # The fit's coefficient must solve the same Breslow weighted Cox problem.
  direct <- survival::coxph(
    survival::Surv(months, died) ~ cohort,
    data = sample_data, weights = fit$weights, ties = "breslow"
  )
  expect_equal(unname(fit$theta), unname(stats::coef(direct)[1L]),
               tolerance = 1e-7)

  requested <- mecCox:::.mec_oof_survival_at_times(
    fit, fit$landmark_times
  )
  expect_equal(requested, unname(fit$survival_oof),
               ignore_attr = TRUE, tolerance = 1e-8)
})

test_that("MLP propensity and RSF survival run together", {
  skip_if_not_installed("brulee")
  skip_if_not_installed("torch")
  skip_if_not_installed("ranger")
  if (!torch::torch_is_installed()) {
    skip("The Torch runtime is not installed")
  }

  set.seed(301)
  n <- 120L
  source <- rep(c(0L, 1L), each = n / 2L)
  sample_data <- data.frame(
    months = rexp(n, rate = 0.07) + 0.01,
    died = rbinom(n, 1L, 0.8),
    cohort = source,
    group = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
    marker = rnorm(n)
  )
  fit <- fit_mec_cox(
    sample_data, "months", "died", "cohort", c("group", "marker"),
    ps_learner = "mlp", survival_learner = "rsf",
    n_folds = 3L, n_landmarks = 3L,
    mlp_epochs = 5L, rsf_num_trees = 50L,
    rsf_auto_tune = FALSE, seed = 301L
  )
  expect_s3_class(fit, "mec_cox_fit")
  expect_true(is.finite(fit$se))
  expect_true(all(fit$weights > 0))
  expect_lt(fit$calibration$max_mean_balance_residual, 1e-6)

  repeated <- fit_mec_cox(
    sample_data, "months", "died", "cohort", c("group", "marker"),
    ps_learner = "mlp", survival_learner = "rsf",
    n_folds = 3L, n_landmarks = 3L,
    mlp_epochs = 5L, rsf_num_trees = 50L,
    rsf_auto_tune = FALSE, seed = 301L
  )
  expect_equal(repeated$weights, fit$weights, tolerance = 1e-10)
  expect_equal(repeated$ps_oof, fit$ps_oof, tolerance = 1e-10)
})

test_that("MEC-Cox score and Cox fit agree for near-tied event times", {
  set.seed(211)
  n <- 100L
  months <- seq(0.5, 12, length.out = n)
  months[25:26] <- c(3, 3 + 1e-12)
  sample_data <- data.frame(
    months = months,
    died = rep(1L, n),
    cohort = sample(rep(c(0L, 1L), each = n / 2L)),
    marker = rnorm(n)
  )
  fit <- fit_mec_cox(
    sample_data, "months", "died", "cohort", "marker",
    n_folds = 3L, n_landmarks = 1L, seed = 211L
  )
  residual_score <- mecCox:::.mec_cox_score(
    fit$theta, fit$analysis_data, fit$weights
  )$score
  expect_lt(abs(residual_score), 1e-6)
})

test_that("RSF tuning rewards correctly ordered mortality risk", {
  months <- 1:12
  died <- rep(1L, length(months))
  correct_risk <- rev(months)
  wrong_risk <- months
  expect_gt(mecCox:::.mec_risk_concordance(
    months, died, correct_risk
  ), 0.99)
  expect_lt(mecCox:::.mec_risk_concordance(
    months, died, wrong_risk
  ), 0.01)
})

test_that("BART propensity works with Cox prognostic calibration", {
  skip_if_not_installed("dbarts")
  set.seed(404)
  n <- 120L
  sample_data <- data.frame(
    months = rexp(n, rate = 0.06) + 0.01,
    died = rbinom(n, 1L, 0.8),
    cohort = rep(c(0L, 1L), each = n / 2L),
    marker = rnorm(n),
    group = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  fit <- fit_mec_cox(
    sample_data, "months", "died", "cohort", c("marker", "group"),
    ps_learner = "bart", survival_learner = "cox",
    n_folds = 3L, n_landmarks = 2L,
    bart_num_trees = 20L, bart_posterior_draws = 40L,
    bart_burn_in = 20L, seed = 404L
  )
  expect_true(all(is.finite(fit$ps_oof)))
  expect_true(all(fit$ps_oof > 0 & fit$ps_oof < 1))
  expect_true(is.finite(fit$theta) && is.finite(fit$se))
  expect_lt(fit$calibration$max_mean_balance_residual, 1e-6)
})

test_that("RSF tuning selects finite fold-specific models", {
  skip_if_not_installed("ranger")
  set.seed(902)
  n <- 120L
  sample_data <- data.frame(
    months = rexp(n, rate = 0.08) + 0.01,
    died = rbinom(n, 1L, 0.8),
    cohort = rep(c(0L, 1L), each = n / 2L),
    marker = rnorm(n),
    group = rbinom(n, 1L, 0.5)
  )
  fit <- fit_mec_cox(
    sample_data, "months", "died", "cohort", c("marker", "group"),
    ps_learner = "glm", survival_learner = "rsf",
    n_folds = 3L, n_landmarks = 2L,
    rsf_num_trees = 50L, rsf_auto_tune = TRUE, seed = 902L
  )
  selected <- lapply(fit$outcome_fits, `[[`, "tuning")
  expect_true(all(vapply(selected, function(item) {
    is.finite(item$mtry) && is.finite(item$min_node_size)
  }, logical(1))))
  expect_lt(fit$calibration$max_mean_balance_residual, 1e-6)
})

test_that("MEC-Cox rejects invalid folds and bad landmark specifications", {
  tiny <- data.frame(
    time = c(1, 2, 3, 4, 5, 6),
    event = c(1, 1, 0, 1, 0, 1),
    source = c(1, 1, 1, 0, 0, 0),
    x = c(0, 1, 0, 1, 0, 1)
  )
  expect_error(
    fit_mec_cox(tiny, "time", "event", "source", "x", n_folds = 4L),
    "at least `n_folds`"
  )
  expect_error(
    fit_mec_cox(tiny, "time", "event", "source", "x",
                n_folds = 2L, landmarks = c(2, 2)),
    "duplicates"
  )
})
