breast_rng_directory <- testthat::test_path("..", "..", "inst", "reproduce")
if (!file.exists(file.path(breast_rng_directory, "breast_cancer_helpers.R"))) {
  breast_rng_directory <- system.file("reproduce", package = "mecCox")
}
breast_rng_helpers <- new.env(parent = baseenv())
source(file.path(breast_rng_directory, "breast_cancer_helpers.R"),
       local = breast_rng_helpers)
for (expression in as.list(parse(file.path(breast_rng_directory,
                                          "breast_cancer.R")))) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.call(expression[[3L]]) &&
      identical(expression[[3L]][[1L]], as.name("function"))) {
    eval(expression, envir = breast_rng_helpers)
  }
}

testthat::test_that("breast fold assignment advances the study random stream", {
  withr::local_rng_version("4.0.0")
  withr::local_seed(701L)
  data <- breast_rng_helpers$prepare_breast_cancer_data()
  fold <- breast_rng_helpers$.breast_folds(data$A, 10L, 20260427L)
  # Reference draw obtained by executing make_stratified_folds() from the
  # original study script on these cohorts. Restoring the pre-fold stream,
  # or resetting the seed before fitting the first MLP, changes this draw.
  testthat::expect_identical(sample.int(100000L, 1L), 57365L)
  testthat::expect_equal(as.integer(table(fold, data$A)[, 1L]),
                         c(rep(265L, 3L), rep(264L, 7L)))
  testthat::expect_equal(as.integer(table(fold, data$A)[, 2L]),
                         c(rep(25L, 6L), rep(24L, 4L)))
})

testthat::test_that("DL/RSF learners follow the original cross-fold RNG sequence", {
  testthat::skip_if_not_installed("brulee")
  testthat::skip_if_not_installed("torch")
  testthat::skip_if_not_installed("ranger")
  withr::local_rng_version("4.0.0")
  withr::local_seed(702L)
  data <- breast_rng_helpers$prepare_breast_cancer_data()
  covariates <- c("age", "meno", "size_cat", "grade", "log_nodes",
                  "log_pgr", "log_er")
  initialization_draws <- integer()
  testthat::local_mocked_bindings(brulee_mlp = function(...) {
    # brulee obtains its initialization seed from R's active random stream.
    initialization_draws <<- c(initialization_draws, sample.int(100000L, 1L))
    # Further fitting can consume random numbers. The original RSF tuning
    # explicitly seeds its split before the next fold's neural fit.
    runif(17L)
    list(rng_test_learner = "mlp")
  }, .package = "brulee")
  testthat::local_mocked_bindings(
    torch_get_num_threads = function() 1L,
    torch_get_rng_state = function() raw(),
    torch_set_num_threads = function(...) invisible(NULL),
    torch_set_rng_state = function(...) invisible(NULL),
    .package = "torch")
  testthat::local_mocked_bindings(ranger = function(...) {
    # The study supplies an explicit ranger seed; forest fitting does not
    # consume the R stream that the following brulee call uses.
    list(rng_test_learner = "rsf")
  }, .package = "ranger")
  testthat::local_mocked_bindings(predict = function(object, new_data = NULL,
                                                    data = NULL, ...) {
    if (identical(object$rng_test_learner, "mlp")) {
      return(data.frame(.pred_1 = rep(0.1, nrow(new_data))))
    }
    list(unique.death.times = c(1, 3, 5),
         survival = outer(stats::plogis(-data$age / 100), c(0.9, 0.8, 0.7)))
  }, .package = "stats")
  # Stop after all learners have run: this test concerns their RNG ordering,
  # independently of calibration, Cox estimation, or neural optimization.
  fit_helpers <- new.env(parent = breast_rng_helpers)
  fit_helpers$fit_breast_mec <- breast_rng_helpers$fit_breast_mec
  environment(fit_helpers$fit_breast_mec) <- fit_helpers
  fit_helpers$.breast_calibrate <- function(...) stop("RNG checkpoint reached")
  before <- .Random.seed
  testthat::expect_error(suppressMessages(fit_helpers$fit_breast_mec(
    data, covariates, ps_learner = "mlp", survival_learner = "rsf",
    seed = 20260427L)), "RNG checkpoint reached", fixed = TRUE)
  # Independently recorded from the original study's make_stratified_folds()
  # and tune_rsf() functions on the public cohorts, with forest fitting
  # mocked as above. Neither the current helper nor a per-fold reseeding
  # scheme was used to generate this reference sequence.
  testthat::expect_identical(initialization_draws,
    c(57365L, 34146L, 99178L, 29247L, 34168L,
      21880L, 86180L, 74012L, 51635L, 57335L))
  testthat::expect_identical(.Random.seed, before)
})
