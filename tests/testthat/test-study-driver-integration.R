# Exercise the real driver/worker plumbing with a cheap stand-in for a study
# fit. The separate study-engine tests verify the scientific function snapshot.
study_driver_directory <- testthat::test_path("..", "..", "inst", "reproduce")
if (!dir.exists(study_driver_directory)) {
  study_driver_directory <- system.file("reproduce", package = "mecCox")
}
study_driver_directory <- normalizePath(study_driver_directory, mustWork = TRUE)

study_driver_fixture <- function(scenario) {
  execution <- new.env(parent = baseenv())
  source(file.path(study_driver_directory, "simulation_helpers.R"),
         local = execution)
  source(file.path(study_driver_directory, "original_study_helpers.R"),
         local = execution)
  expressions <- as.list(parse(file.path(study_driver_directory,
                                         paste0(scenario, ".R"))))
  for (expression in expressions) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        is.call(expression[[3L]]) &&
        identical(expression[[3L]][[1L]], as.name("function"))) {
      eval(expression, execution)
    }
  }
  key <- if (scenario == "scenario1") "2" else "none"
  config <- execution$get_original_study_config(scenario, key)
  execution$study_configs <- stats::setNames(list(config), key)
  execution$study_streams <- stats::setNames(
    list(execution$original_study_streams(config, replicates = 3L)), key
  )
  execution$study_reference <- execution$original_study_metadata(scenario)
  execution$study_engine <- new.env(parent = baseenv())
  execution$study_engine$method_labels <- c(
    "ATT-IPW Cox: Naive model-based", "ATT-IPW Cox: Lin-Wei", "ATT-IPW Cox: Shu",
    "MEC-Cox: PS=bart, OR=cox, L=5, K=10, kl"
  )
  if (scenario == "scenario2") {
    execution$study_engine$method_labels <- c(
      execution$study_engine$method_labels,
      "MEC-Cox: PS=bart, OR=rsf, L=5, K=10, kl"
    )
  }
  execution$study_engine$fit_one_rep <- function(n1, n0, rep_id, theta_true,
                                                  trim, auto_tune_seed) {
    estimate <- stats::rnorm(1L)
    data.frame(
      rep = rep_id, n1 = n1, n0 = n0, Method = method_labels,
      theta_true = theta_true, Estimate = estimate, SE = 0.1,
      CI_L = estimate - 0.196, CI_U = estimate + 0.196,
      ESS = 100, Rel_ESS = 0.5, W_CV = 1, W_Min = 0.1, W_Max = 2,
      Cal_Grad = 1e-8, Cal_Converged = TRUE, Error = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  environment(execution$study_engine$fit_one_rep) <- execution$study_engine
  execution$design <- list(
    replications = 3L, treated_sizes = config$sample_size_n1[1:2],
    control_multipliers = 2L, control_multiplier = config$n0_over_n1
  )
  execution$settings <- data.frame(setting_id = 1L, setting_key = key,
                                    setting = "No nonlinearity",
                                    kappa_pi = 0, kappa_m = 0)
  execution$options <- list(quick_run = TRUE, cores = 1L, replications = 3L)
  execution
}

testthat::test_that("driver dispatch selects the study grid stream and keeps raw intervals", {
  testthat::skip_if_not_installed("rngtools")
  for (scenario in c("scenario1", "scenario2")) {
    execution <- study_driver_fixture(scenario)
    plan <- execution$study_streams[[1L]]
    config <- execution$study_configs[[1L]]
    # In the five-size grid, the second replicate at the second size is job 7.
    testthat::expect_equal(plan$grid$job_id[
      plan$grid$n1 == 250L & plan$grid$rep_id == 2L], 7L)
    arguments <- list(treated_count = 250L, replicate = 2L,
                       target = 999, design = execution$design)
    if (scenario == "scenario1") arguments$multiplier <- 2L
    else arguments$setting <- execution$settings
    result <- do.call(execution$run_replication, arguments)
    reference <- execution$fit_original_study_rep(
      execution$study_engine, config, 250L, 2L, plan$streams[[7L]]
    )
    testthat::expect_identical(result$estimate, reference$Estimate)
    testthat::expect_identical(result$standard_error, reference$SE)
    testthat::expect_identical(result$ci_lower, reference$CI_L)
    testthat::expect_identical(result$ci_upper, reference$CI_U)
    testthat::expect_identical(result$target, reference$theta_true)
    testthat::expect_true(all(result$original_rng_job == 7L))
    testthat::expect_identical(result$original_method, reference$Method)
    testthat::expect_identical(result$Cal_Converged, reference$Cal_Converged)
    testthat::expect_identical(result$method[1:3],
                               c("Naive", "Robust sandwich", "Corrected sandwich"))
    expected_mec <- if (scenario == "scenario1") "MEC-Cox" else {
      c("MEC-Cox (BART/Cox)", "MEC-Cox (BART/RSF)")
    }
    testthat::expect_identical(result$method[-(1:3)], expected_mec)
  }
})

testthat::test_that("summaries use supplied confidence limits and finite fit criteria", {
  testthat::skip_if_not_installed("rngtools")
  for (scenario in c("scenario1", "scenario2")) {
    execution <- study_driver_fixture(scenario)
    if (scenario == "scenario1") {
      rows <- execution$make_result(2L, 200L, 1:3, "MEC-Cox",
                                    c(1.95998, 0, 0), c(0, 0, NA_real_), c(1, 0, 1))
    } else {
      rows <- execution$make_result(execution$settings, 200L, 1:3,
                                    "MEC-Cox (BART/Cox)", c(1.95998, 0, 0),
                                    execution$design, c(0, 0, NA_real_), c(1, 0, 1))
    }
    # The first target is covered by 1.96 but not by qnorm(.975); a finite
    # zero-SE row also follows the study's existing inclusion criterion.
    result <- if (scenario == "scenario1") execution$summarize_results(rows) else {
      execution$summarize_scenario2_results(rows, execution$settings)
    }
    testthat::expect_equal(result$successful, 2L)
    testthat::expect_equal(result$failed, 1L)
    testthat::expect_equal(result$coverage, 1)
  }
})

testthat::test_that("bundled study streams give identical serial and PSOCK driver output", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rngtools")
  for (scenario in c("scenario1", "scenario2")) {
    execution <- study_driver_fixture(scenario)
    run <- function() {
      if (scenario == "scenario1") {
        execution$run_scenario1(execution$design, execution$options)
      } else {
        execution$run_scenario2(execution$design, execution$settings,
                                execution$options)
      }
    }
    invisible(capture.output(serial <- suppressMessages(run())))
    execution$options$cores <- 2L
    invisible(capture.output(concurrent <- suppressMessages(run())))
    testthat::expect_identical(serial$status, "completed")
    testthat::expect_identical(concurrent$status, "completed")
    testthat::expect_identical(concurrent$replications, serial$replications)
    testthat::expect_identical(concurrent$summary, serial$summary)
    testthat::expect_identical(concurrent$targets, serial$targets)
    testthat::expect_identical(concurrent$metadata$study_engine,
                               execution$study_reference)
    testthat::expect_equal(concurrent$metadata$execution$completed_replications, 6L)
    per_job <- unique(concurrent$replications[c("n1", "replicate", "original_rng_job")])
    testthat::expect_equal(per_job$original_rng_job, c(1, 6, 11, 2, 7, 12))
  }
})

testthat::test_that("Scenario 1 preserves completed datasets after interruption", {
  testthat::skip_if_not_installed("rngtools")
  execution <- study_driver_fixture("scenario1")
  execution$run_replication <- function(multiplier, treated_count, replicate,
                                         target, design) {
    if (treated_count == 250L && replicate == 2L) {
      stop(structure(list(message = "test interruption", call = NULL),
                     class = c("interrupt", "condition")))
    }
    make_result(multiplier, treated_count, replicate, "MEC-Cox", target,
                  estimate = 0, standard_error = 1)
  }
  environment(execution$run_replication) <- execution
  invisible(capture.output(answer <- suppressMessages(execution$run_scenario1(
    execution$design, execution$options
  ))))
  testthat::expect_identical(answer$status, "interrupted")
  testthat::expect_equal(answer$replications$n1, c(200, 200, 200, 250))
  testthat::expect_identical(answer$replications$replicate, c(1:3, 1L))
  testthat::expect_equal(answer$summary$replications, c(3, 1))
  testthat::expect_equal(answer$metadata$execution$completed_replications, 4L)
  testthat::expect_equal(answer$metadata$execution$planned_replications, 6L)
})
