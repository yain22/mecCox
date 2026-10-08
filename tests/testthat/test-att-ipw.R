make_att_example <- function() {
  set.seed(20261008)
  sample_size <- 180L
  age_score <- stats::rnorm(sample_size)
  marker <- stats::rbinom(sample_size, 1L, 0.45)
  source_probability <- stats::plogis(
    -0.25 + 0.40 * age_score + 0.50 * marker
  )
  source <- stats::rbinom(sample_size, 1L, source_probability)

  event_time <- stats::rexp(
    sample_size,
    rate = 0.45 * exp(0.20 * source + 0.30 * age_score)
  )
  censor_time <- stats::rexp(sample_size, rate = 0.13)
  observed_time <- round(pmin(event_time, censor_time), 2) + 0.01

  data.frame(
    follow_up = observed_time,
    death = as.integer(event_time <= censor_time),
    cohort = source,
    age_score = age_score,
    marker = marker
  )
}

testthat::test_that("ATT-IPW uses one Cox estimate with three variance methods", {
  example <- make_att_example()
  fit <- fit_att_ipw_cox(
    example,
    time = "follow_up",
    event = "death",
    source = "cohort",
    covariates = c("age_score", "marker")
  )

  testthat::expect_s3_class(fit, "att_ipw_cox_fit")
  testthat::expect_named(fit$se, c("naive", "robust", "corrected"))
  testthat::expect_equal(fit$hr, exp(fit$theta), tolerance = 1e-12)
  testthat::expect_true(all(is.finite(fit$se) & fit$se > 0))
  testthat::expect_equal(fit$weights[example$cohort == 1L],
                         rep(1, sum(example$cohort == 1L)))
  testthat::expect_equal(
    sum(fit$weights[example$cohort == 0L]),
    sum(example$cohort == 1L),
    tolerance = 1e-10
  )
  testthat::expect_equal(fit$diagnostics$clipped_probability_count, 0L)

  independent_data <- data.frame(
    time = example$follow_up,
    event = example$death,
    source = example$cohort,
    weight = fit$weights
  )
  naive_cox <- survival::coxph(
    survival::Surv(time, event) ~ source,
    data = independent_data,
    weights = weight,
    robust = FALSE,
    ties = "breslow",
    control = survival::coxph.control(timefix = FALSE)
  )
  robust_cox <- survival::coxph(
    survival::Surv(time, event) ~ source,
    data = independent_data,
    weights = weight,
    robust = TRUE,
    ties = "breslow",
    control = survival::coxph.control(timefix = FALSE)
  )

  testthat::expect_equal(fit$theta, unname(stats::coef(naive_cox)[1]),
                         tolerance = 1e-10)
  testthat::expect_equal(fit$se[["naive"]],
                         sqrt(unname(stats::vcov(naive_cox)[1, 1])),
                         tolerance = 1e-10)
  testthat::expect_equal(fit$se[["robust"]],
                         sqrt(unname(stats::vcov(robust_cox)[1, 1])),
                         tolerance = 1e-8)
  testthat::expect_equal(
    fit$conf_int$hazard_lower,
    exp(fit$theta - stats::qnorm(0.975) * fit$se),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
})

testthat::test_that("the corrected sandwich is reproducible without clipping", {
  example <- make_att_example()
  fit <- fit_att_ipw_cox(
    example, "follow_up", "death", "cohort",
    c("age_score", "marker")
  )

  # Regression values are compared with the original manuscript's
  # fit_att_ipw_all_variances() on this fixed tied-event example.
  testthat::expect_equal(fit$theta, 0.2240078140326363,
                         tolerance = 1e-8)
  testthat::expect_equal(fit$se[["corrected"]], 0.1685415907587880,
                         tolerance = 1e-8)
})

testthat::test_that("clipping changes weights but not the fitted GLM score", {
  example <- make_att_example()
  fit <- fit_att_ipw_cox(
    example, "follow_up", "death", "cohort",
    c("age_score", "marker"),
    ps_clip = c(0.45, 0.55)
  )

  testthat::expect_gt(fit$diagnostics$clipped_probability_count, 0L)
  testthat::expect_true(all(is.finite(fit$se)))
  testthat::expect_equal(
    fit$ps,
    unname(stats::predict(fit$ps_fit, type = "response")),
    tolerance = 1e-12
  )
  testthat::expect_true(all(fit$ps_clipped >= 0.45 &
                              fit$ps_clipped <= 0.55))
  testthat::expect_equal(
    sum(fit$weights[example$cohort == 0L]),
    sum(example$cohort == 1L),
    tolerance = 1e-10
  )
})
