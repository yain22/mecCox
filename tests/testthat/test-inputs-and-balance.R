testthat::test_that("unweighted comparison agrees with the Breslow Cox reference", {
  data(example_external_controls, package = "mecCox")
  fit <- fit_unweighted_cox(
    example_external_controls, "time", "event", "source"
  )
  reference <- survival::coxph(
    survival::Surv(time, event) ~ source,
    data = example_external_controls,
    ties = "breslow", robust = FALSE
  )

  testthat::expect_equal(fit$theta, unname(stats::coef(reference)[1]),
                         tolerance = 1e-12)
  testthat::expect_equal(fit$se,
                         sqrt(unname(stats::vcov(reference)[1, 1])),
                         tolerance = 1e-12)
})

testthat::test_that("malformed survival and source inputs fail clearly", {
  data(example_external_controls, package = "mecCox")
  invalid <- example_external_controls
  invalid$event[1] <- 2L
  testthat::expect_error(
    fit_unweighted_cox(invalid, "time", "event", "source"),
    "event.*0 and 1"
  )

  invalid <- example_external_controls
  invalid$source <- 1L
  testthat::expect_error(
    fit_unweighted_cox(invalid, "time", "event", "source"),
    "Both source groups"
  )

  invalid <- example_external_controls
  invalid$event[invalid$source == 0L] <- 0L
  testthat::expect_error(
    fit_unweighted_cox(invalid, "time", "event", "source"),
    "event.*each source group"
  )

  invalid <- example_external_controls
  invalid$age[1] <- NA_real_
  testthat::expect_error(
    fit_att_ipw_cox(invalid, "time", "event", "source",
                    c("age", "sex", "marker")),
    "age.*finite"
  )
})

testthat::test_that("MEC calibration balances the actual prediction basis", {
  data(example_external_controls, package = "mecCox")
  fit <- fit_mec_cox(
    example_external_controls,
    "time", "event", "source",
    c("age", "sex", "marker"),
    n_folds = 4, n_landmarks = 4,
    seed = 1001
  )

  testthat::expect_lt(
    fit$calibration$max_mean_balance_residual, 1e-7
  )
  risk_balance <- predicted_risk_balance(
    fit, times = fit$landmark_times
  )
  calibrated_rows <- risk_balance$weighting == "mec"
  testthat::expect_lt(max(abs(risk_balance$gap[calibrated_rows])), 1e-7)

  integer_balance <- predicted_risk_balance(fit, times = c(3, 6, 12))
  testthat::expect_equal(sort(unique(integer_balance$time)), c(3, 6, 12))
  testthat::expect_true(all(is.finite(integer_balance$gap)))

  covariates <- covariate_balance(fit)
  testthat::expect_true(all(c("unweighted", "source_odds", "mec") %in%
                              covariates$weighting))
  binary_rows <- covariates$variable_type == "binary"
  continuous_rows <- covariates$variable_type == "continuous"
  testthat::expect_equal(
    covariates$reported_difference[binary_rows],
    covariates$difference[binary_rows]
  )
  testthat::expect_equal(
    covariates$reported_difference[continuous_rows],
    covariates$standardized_difference[continuous_rows]
  )
})
