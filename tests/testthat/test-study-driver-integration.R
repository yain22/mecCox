# Exercise the real driver/worker plumbing with a cheap stand-in for a study
# fit. The separate study-engine tests verify the scientific function snapshot.
study_driver_directory <- testthat::test_path("..", "..", "inst", "reproduce")
if (!dir.exists(study_driver_directory)) {
  study_driver_directory <- system.file("reproduce", package = "mecCox")
}
study_driver_directory <- normalizePath(study_driver_directory, mustWork = TRUE)

source(testthat::test_path("helper-final-driver-fixtures.R"), local = TRUE)
study_driver_fixture <- final_driver_test_fixture

testthat::test_that("driver dispatch selects the study grid stream and keeps raw intervals", {
  for (scenario in c("scenario1", "scenario2")) {
    execution <- study_driver_fixture(scenario)
    reference_bundle <- execution$final_reference$bundle
    job_id <- if (scenario == "scenario1") 7L else 15007L
    # A quick subset retains the full frozen grid's global ID and state.
    testthat::expect_equal(nrow(reference_bundle$grid), 30000L)
    arguments <- list(treated_count = 250L, replicate = 2L,
                       target = 999, design = execution$design)
    if (scenario == "scenario1") arguments$multiplier <- 2L
    else arguments$setting <- execution$settings
    result <- do.call(execution$run_replication, arguments)
    reference_run <- execution$final_reference$api$fit_final_study_job(
      job_id, reference_bundle)
    reference <- execution$final_reference$api$bind_final_study_target(
      reference_run, reference_bundle)
    testthat::expect_identical(result$estimate, reference$Estimate)
    testthat::expect_identical(result$standard_error, reference$SE)
    testthat::expect_identical(result$ci_lower, reference$CI_L)
    testthat::expect_identical(result$ci_upper, reference$CI_U)
    testthat::expect_identical(result$target, reference$theta_true)
    testthat::expect_true(all(result$original_rng_job == job_id))
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
    offset <- if (scenario == "scenario1") 0L else 15000L
    testthat::expect_equal(per_job$original_rng_job, offset + c(1, 6, 11, 2, 7, 12))
    testthat::expect_identical(concurrent$diagnostics, serial$diagnostics)
  }
})

testthat::test_that("Scenario 1 preserves completed datasets after interruption", {
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
