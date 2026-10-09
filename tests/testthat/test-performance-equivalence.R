test_that("score-only evaluation preserves tied-event scores and derivatives", {
  prepared <- data.frame(
    time = c(1, 1, 2, 2, 2 + 1e-12, 3, 4, 4),
    delta = c(1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L),
    A = c(0L, 1L, 1L, 0L, 1L, 0L, 0L, 1L)
  )
  weights <- c(0.3, 1, 1, 1.8, 1, 0.7, 1.2, 1)
  for (score_function in list(mecCox:::.att_cox_score,
                             mecCox:::.mec_cox_score)) {
    for (theta in c(-0.5, 0, 0.4)) {
      full <- score_function(theta, prepared, weights)
      score_only <- score_function(theta, prepared, weights,
                                   individual = FALSE)
      expect_identical(score_only$score, full$score)
      expect_named(score_only, "score")
    }
    # Check the perturbations in both the Cox and weight-estimation directions.
    evaluate <- function(parameter, individual) {
      candidate_weights <- weights * exp(parameter[2L] * (1 - prepared$A))
      score_function(parameter[1L], prepared, candidate_weights,
                     individual = individual)$score
    }
    full_derivative <- mecCox:::.mec_central_difference(
      function(parameter) evaluate(parameter, TRUE), c(0.2, 0.1)
    )
    score_derivative <- mecCox:::.mec_central_difference(
      function(parameter) evaluate(parameter, FALSE), c(0.2, 0.1)
    )
    expect_identical(score_derivative, full_derivative)
  }
})

test_that("reusing BART predictions preserves RSF estimates and RNG state", {
  skip_if_not_installed("dbarts")
  skip_if_not_installed("ranger")
  set.seed(630)
  n <- 150L
  sample_data <- data.frame(
    months = rexp(n, rate = 0.08) + 0.01,
    died = rbinom(n, 1L, 0.8),
    cohort = rep(c(0L, 1L), each = n / 2L),
    marker = rnorm(n),
    group = rbinom(n, 1L, 0.5)
  )
  arguments <- list(
    data = sample_data, time = "months", event = "died", source = "cohort",
    covariates = c("marker", "group"), ps_learner = "bart",
    n_folds = 3L, n_landmarks = 1L, seed = 630L,
    ps_trim = c(0.45, 0.55), bart_num_trees = 5L,
    bart_posterior_draws = 12L, bart_burn_in = 6L,
    rsf_num_trees = 20L, rsf_auto_tune = FALSE
  )
  cox_fit <- do.call(fit_mec_cox, c(arguments, list(survival_learner = "cox")))
  set.seed(808)
  fresh <- do.call(fit_mec_cox, c(arguments, list(survival_learner = "rsf")))
  fresh_rng <- .Random.seed

  local_mocked_bindings(
    .mec_propensity_bart = function(...) stop("BART must not be refitted"),
    .package = "mecCox"
  )
  set.seed(808)
  reused <- do.call(fit_mec_cox, c(arguments, list(
    survival_learner = "rsf", ps_predictions = cox_fit$ps_oof
  )))
  expect_identical(.Random.seed, fresh_rng)
  for (component in c("ps_oof", "fold_id", "survival_oof", "baseline_weights",
                      "weights", "theta", "se")) {
    expect_identical(reused[[component]], fresh[[component]])
  }
  expect_false(fresh$ps_predictions_supplied)
  expect_true(reused$ps_predictions_supplied)
  expect_true(all(vapply(reused$outcome_fits, function(item) {
    is.null(item$fit$survival) && is.null(item$fit$chf)
  }, logical(1))))
  predictions <- mecCox:::.mec_oof_survival_at_times(
    reused, reused$landmark_times
  )
  expect_equal(predictions, unname(reused$survival_oof), ignore_attr = TRUE,
               tolerance = 1e-12)
})

test_that("supplied propensity vectors are checked before fitting", {
  sample_data <- data.frame(
    time = 1:6, event = c(1L, 1L, 0L, 1L, 0L, 1L),
    source = rep(c(0L, 1L), each = 3L), x = c(0, 1, 0, 1, 0, 1)
  )
  invalid <- list(
    rep(0.5, 5L), rep("0.5", 6L), rep(TRUE, 6L), matrix(0.5, 6L, 1L),
    c(NA_real_, rep(0.5, 5L)), c(Inf, rep(0.5, 5L)),
    c(0, rep(0.5, 5L)), c(1, rep(0.5, 5L)),
    c(-0.1, rep(0.5, 5L)), c(1.1, rep(0.5, 5L))
  )
  for (probabilities in invalid) {
    expect_error(fit_mec_cox(
      sample_data, "time", "event", "source", "x", n_folds = 2L,
      n_landmarks = 1L, ps_predictions = probabilities
    ), "`ps_predictions` must be")
  }
})

test_that("supplied propensities receive the normal trimming rule", {
  data(example_external_controls)
  n <- nrow(example_external_controls)
  probabilities <- rep(c(0.001, 0.999, 0.5), length.out = n)
  fit <- fit_mec_cox(
    example_external_controls, "time", "event", "source",
    c("age", "sex", "marker"), n_folds = 3L, n_landmarks = 1L,
    seed = 12L, ps_trim = c(0.1, 0.9), ps_predictions = probabilities
  )
  expect_identical(fit$ps_oof, pmin(pmax(probabilities, 0.1), 0.9))
})

test_that("the original-study preset tunes small learners and records fold choices", {
  skip_if_not_installed("dbarts")
  skip_if_not_installed("ranger")
  set.seed(917)
  n <- 240L
  sample_data <- data.frame(
    time = rexp(n, 0.08) + 0.01,
    event = rbinom(n, 1L, 0.8),
    source = rep(c(0L, 1L), each = n / 2L)
  )
  covariates <- paste0("x", seq_len(10L))
  sample_data[covariates] <- lapply(covariates, function(variable) rnorm(n))
  arguments <- list(
    data = sample_data, time = "time", event = "event", source = "source",
    covariates = covariates, ps_learner = "bart", survival_learner = "rsf",
    n_folds = 3L, n_landmarks = 1L, seed = 917L,
    ps_trim = c(0.2, 0.8), nuisance_settings = "original_study"
  )
  set.seed(181)
  fresh <- do.call(fit_mec_cox, arguments)
  fresh_rng <- .Random.seed
  expect_true(fresh$bart_auto_tune)
  expect_identical(fresh$nuisance_settings, "original_study")
  expect_true(is.finite(fresh$theta) && is.finite(fresh$se))
  expect_lt(fresh$calibration$max_mean_balance_residual, 1e-6)

  for (fold in seq_len(3L)) {
    propensity <- fresh$propensity_tuning[[fold]]
    expect_equal(propensity$tuning$candidates$num_trees, c(25, 50, 100))
    expect_equal(propensity$tuning$tuning_rows, 100)
    expect_equal(propensity$tuning$training_rows, 70)
    expect_equal(propensity$tuning$validation_rows, 30)
    expect_equal(propensity$parameters$posterior_draws, 100)
    expect_equal(propensity$parameters$burn_in, 50)
    expect_equal(propensity$parameters$shrinkage, 2)
    expect_true(propensity$parameters$num_trees %in% c(25, 50, 100))
    expect_equal(propensity$tuning$seed, 917 + fold)

    outcome <- fresh$outcome_fits[[fold]]
    expect_equal(sort(unique(outcome$tuning$candidates$mtry)), c(3, 4, 5))
    expect_equal(sort(unique(outcome$tuning$candidates$min_node_size)),
                 c(15, 30, 50))
    # Each training fold has only 80 controls, so tuning uses all 80.
    expect_equal(outcome$tuning$tuning_rows, 80)
    expect_true(outcome$tuning$training_rows %in% c(39, 40))
    expect_equal(outcome$tuning$candidate_num_trees, 100)
    expect_equal(outcome$parameters$num_trees, 100)
    expect_identical(outcome$parameters$splitrule, "extratrees")
    expect_false(outcome$parameters$replace)
    expect_equal(outcome$parameters$sample_fraction, 0.632)
    expect_equal(outcome$parameters$num_random_splits, 1)
    expect_equal(outcome$parameters$seed, 917 + fold)
  }

  local_mocked_bindings(
    .mec_propensity_bart = function(...) stop("BART must not be refitted"),
    .package = "mecCox"
  )
  set.seed(181)
  reused <- do.call(fit_mec_cox, c(arguments, list(ps_predictions = fresh$ps_oof)))
  expect_identical(.Random.seed, fresh_rng)
  for (component in c("ps_oof", "fold_id", "survival_oof", "baseline_weights",
                      "weights", "theta", "se")) {
    expect_identical(reused[[component]], fresh[[component]])
  }
})

test_that("original-study quick controls and tuning arguments are respected", {
  skip_if_not_installed("dbarts")
  skip_if_not_installed("ranger")
  data(example_external_controls)
  arguments <- list(
    data = example_external_controls, time = "time", event = "event",
    source = "source", covariates = c("age", "sex", "marker"),
    ps_learner = "bart", survival_learner = "rsf", n_folds = 3L,
    n_landmarks = 1L, seed = 31L, nuisance_settings = "original_study",
    bart_auto_tune = FALSE, rsf_auto_tune = FALSE, bart_num_trees = 5L,
    bart_posterior_draws = 12L, bart_burn_in = 6L, rsf_num_trees = 20L
  )
  fit <- do.call(fit_mec_cox, arguments)
  for (fold in seq_len(3L)) {
    expect_null(fit$propensity_tuning[[fold]]$tuning)
    expect_equal(fit$propensity_tuning[[fold]]$parameters$num_trees, 5)
    expect_equal(fit$propensity_tuning[[fold]]$parameters$posterior_draws, 12)
    expect_equal(fit$propensity_tuning[[fold]]$parameters$burn_in, 6)
    expect_null(fit$outcome_fits[[fold]]$tuning)
    expect_equal(fit$outcome_fits[[fold]]$parameters$num_trees, 20)
    expect_identical(fit$outcome_fits[[fold]]$parameters$splitrule, "extratrees")
  }
  arguments$bart_auto_tune <- NA
  expect_error(do.call(fit_mec_cox, arguments), "`bart_auto_tune` must be")
  arguments$bart_auto_tune <- FALSE
  for (invalid in list(0, -1, 1.5, NA, Inf, TRUE, c(50, 100))) {
    arguments$ml_tune_n <- invalid
    expect_error(do.call(fit_mec_cox, arguments), "`ml_tune_n` must be")
  }
})

test_that("original-study RSF tuning falls back when the subset lacks factor levels", {
  skip_if_not_installed("ranger")
  training <- data.frame(
    time = seq_len(220L), delta = 1L, A = 0L,
    x = factor(c(rep("a", 100L), rep("b", 120L)))
  )
  expect_warning(
    result <- mecCox:::.mec_predict_rsf_survival(
      training, training[c(1L, 220L), ], "x", landmarks = 80,
      seed = 91L, num_trees = 100L, min_node_size = 15L,
      auto_tune = TRUE, nuisance_settings = "original_study"
    ), "tuning subset has inadequate factor levels"
  )
  expect_equal(result$parameters$num_trees, 300)
  expect_equal(result$parameters$min_node_size, 15)
  expect_equal(result$parameters$mtry, 1)
  expect_true(all(is.finite(result$survival)))
  expect_true(all(result$survival >= 0 & result$survival <= 1))
})
