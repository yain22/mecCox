source(testthat::test_path("helper-final-driver-fixtures.R"), local = TRUE)

testthat::test_that("final drivers use frozen IDs for every setting and subset", {
  execution <- final_driver_test_fixture("scenario2")
  reference <- execution$final_reference
  for (setting in c("none", "mild", "severe")) {
    for (n1 in c(200L, 300L, 400L)) {
      jobs <- lapply(c(1L, 2L, 1000L), function(rep) {
        execution$final_driver_job(reference, "figure6", setting, n1, rep)
      })
      expected <- reference$bundle$grid[
        reference$bundle$grid$figure == "figure6" &
          reference$bundle$grid$setting == setting &
          reference$bundle$grid$n1 == n1 &
          reference$bundle$grid$rep_id %in% c(1L, 2L, 1000L), ]
      testthat::expect_identical(vapply(jobs, function(x) x$id, integer(1)),
                                 expected$id)
    }
  }
  for (count in c(1L, 2L, 1000L)) {
    testthat::expect_invisible(execution$check_final_driver_replications(count))
  }
  for (count in list(0, -1, 1.5, NA_real_, Inf, 1001, 1e20, "2", c(1, 2))) {
    testthat::expect_error(execution$check_final_driver_replications(count),
                           "1 to 1000")
  }
  testthat::expect_error(execution$final_driver_job(
    reference, "figure6", "severe", 300L, 1001L), "frozen final grid")
})

testthat::test_that("changing quick subset length preserves common runs and diagnostics", {
  for (scenario in c("scenario1", "scenario2")) {
    execution <- final_driver_test_fixture(scenario)
    run <- function(count) {
      execution$design$replications <- count
      execution$options$replications <- count
      invisible(capture.output(answer <- suppressMessages(
        if (scenario == "scenario1") execution$run_scenario1(
          execution$design, execution$options) else execution$run_scenario2(
            execution$design, execution$settings, execution$options))))
      answer
    }
    one <- run(1L)
    three <- run(3L)
    subset <- three$replications[three$replications$replicate == 1L, ]
    rownames(subset) <- NULL
    testthat::expect_identical(one$replications, subset)
    ids <- vapply(one$diagnostics, function(x) x$job$id, integer(1))
    matched <- Filter(function(x) x$job$id %in% ids, three$diagnostics)
    testthat::expect_identical(one$diagnostics, matched)
    testthat::expect_true(all(vapply(one$diagnostics, function(x)
      all(is.na(x$results$theta_true)), logical(1))))
    testthat::expect_identical(one$diagnostics[[1L]]$warnings,
                               "retained fixture warning")
    testthat::expect_true(all(is.finite(one$replications$target)))
  }
})

testthat::test_that("fatal datasets retain all expected rows and complete diagnostics", {
  for (scenario in c("scenario1", "scenario2")) {
    execution <- final_driver_test_fixture(scenario)
    api <- execution$final_reference$api
    api$successful_stub <- api$fit_final_study_job
    api$fit_final_study_job <- function(job_id, bundle, capture_full = FALSE) {
      run <- successful_stub(job_id, bundle, capture_full)
      run$results[c("Estimate", "SE", "CI_L", "CI_U")] <- NA_real_
      run$results$Error <- "Fatal dataset error: fixture failure"
      run$fatal_error <- "fixture failure"
      run$estimator_trace <- run$math_trace <- list()
      run
    }
    environment(api$fit_final_study_job) <- api
    execution$design$replications <- 1L
    execution$design$treated_sizes <- 200L
    invisible(capture.output(answer <- suppressMessages(
      if (scenario == "scenario1") execution$run_scenario1(
        execution$design, execution$options) else execution$run_scenario2(
          execution$design, execution$settings, execution$options))))
    expected <- if (scenario == "scenario1") 4L else 5L
    testthat::expect_identical(answer$status, "completed")
    testthat::expect_equal(nrow(answer$replications), expected)
    testthat::expect_true(all(answer$summary$failed == 1L))
    testthat::expect_true(all(answer$summary$coverage_denominator == 0L))
    testthat::expect_true(all(is.na(answer$summary$coverage)))
    testthat::expect_identical(answer$diagnostics[[1L]]$fatal_error,
                               "fixture failure")
    testthat::expect_true(all(answer$replications$fatal_error == "fixture failure"))
  }
})

testthat::test_that("finite-fit summaries keep zero SE and expose invalid CI denominators", {
  for (scenario in c("scenario1", "scenario2")) {
    execution <- final_driver_test_fixture(scenario)
    rows <- if (scenario == "scenario1") execution$make_result(
      2L, 200L, 1:6, "MEC-Cox", 0, 0, c(0, .1, .1, .1, .1, NA)) else
      execution$make_result(execution$settings, 200L, 1:6, "MEC-Cox (BART/Cox)",
                            0, execution$design, 0, c(0, .1, .1, .1, .1, NA))
    rows$ci_lower <- c(0, .1, 1, NA, -Inf, -1)
    rows$ci_upper <- c(0, .2, -1, .2, Inf, 1)
    summary <- if (scenario == "scenario1") execution$summarize_results(rows) else
      execution$summarize_scenario2_results(rows, execution$settings)
    testthat::expect_equal(summary$successful, 5L)
    testthat::expect_equal(summary$failed, 1L)
    testthat::expect_equal(summary$coverage_numerator, 1L)
    testthat::expect_equal(summary$coverage_denominator, 2L)
    testthat::expect_equal(summary$coverage, .5)
    testthat::expect_equal(summary$valid_and_covers_per_attempt, 1 / 6)
    testthat::expect_equal(summary$bias, 0)
    testthat::expect_equal(summary$rmse, 0)
  }
})

testthat::test_that("one socket worker serializes the final API closure and frozen bundle", {
  testthat::skip_on_cran()
  cluster <- parallel::makePSOCKcluster(1L)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  for (scenario in c("scenario1", "scenario2")) {
    execution <- final_driver_test_fixture(scenario)
    arguments <- list(treated_count = 250L, replicate = 2L, target = 999,
                       design = execution$design)
    if (scenario == "scenario1") arguments$multiplier <- 2L
    else arguments$setting <- execution$settings
    serial <- do.call(execution$run_replication, arguments)
    socket <- parallel::clusterCall(cluster, function(f, arguments) {
      do.call(f, arguments)
    }, execution$run_replication, arguments)[[1L]]
    testthat::expect_identical(socket, serial)
    testthat::expect_identical(attr(socket, "final_diagnostics"),
                               attr(serial, "final_diagnostics"))
  }
})
