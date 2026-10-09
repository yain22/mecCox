study_helper_path <- testthat::test_path("..", "..", "inst", "reproduce",
                                        "original_study_helpers.R")
if (!file.exists(study_helper_path)) {
  study_helper_path <- system.file("reproduce", "original_study_helpers.R",
                                  package = "mecCox")
}
study_helpers <- new.env(parent = baseenv())
source(study_helper_path, local = study_helpers)

test_that("study engines contain function definitions and isolated lookups", {
  for (scenario in c("scenario1", "scenario2")) {
    engine <- study_helpers$get_original_study_engine(scenario)
    metadata <- study_helpers$original_study_metadata(scenario)
    expressions <- parse(file.path(dirname(study_helper_path), "study_reference",
                                    metadata$engine_basename))
    expect_length(expressions, 40L)
    expect_true(all(vapply(expressions, function(expression) {
      is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        is.call(expression[[3L]]) &&
        identical(expression[[3L]][[1L]], as.name("function"))
    }, logical(1))))
    expect_identical(parent.env(engine), asNamespace("stats"))
    expect_true(all(metadata$function_names %in% ls(engine, all.names = TRUE)))
    expect_false(exists("run_simulation_parallel", envir = engine,
                        inherits = FALSE))
    expect_false(identical(engine$fit_mec_cox, mecCox::fit_mec_cox))
  }
})

test_that("study settings retain full precision and the complete designs", {
  for (ratio in 2:4) {
    config <- study_helpers$get_original_study_config("scenario1", ratio)
    expect_identical(config$n0_over_n1, as.double(ratio))
    expect_identical(config$sample_size_n1, c(200, 250, 300, 350, 400))
    expect_identical(config$sample_size_grid$n0, config$sample_size_n1 * ratio)
    expect_identical(config$lambda0, 0.00008)
    expect_identical(config$censor_rate, 0.0008)
    expect_identical(config$beta_cond, log(0.7))
    expect_equal(config$Kfold_mec, 10)
    expect_equal(config$n_landmarks, 5)
    expect_equal(config$R, 1000)
  }
  expected <- list(none=c(0,0), mild=c(1,2), severe=c(2,5))
  for (setting in names(expected)) {
    config <- study_helpers$get_original_study_config("scenario2", setting)
    expect_identical(c(config$ps_nonlinearity, config$or_nonlinearity),
                     expected[[setting]])
    expect_identical(config$mec_or_method, c("cox", "rsf"))
    expect_identical(config$mec_ps_method, "bart")
    expect_equal(config$n0_over_n1, 4)
    expect_equal(config$auto_tune_n, 100)
    expect_true(config$auto_tune_ml)
    expect_match(attr(config, "study_reference")$source_sha256, "^[0-9a-f]{64}$")
  }
  expect_identical(study_helpers$get_original_study_config("scenario1", 2)$theta_true,
                   -0.22377865980151568)
  expect_identical(study_helpers$get_original_study_config("scenario2", "none")$theta_true,
                   -0.21345570006135195)
})

test_that("study streams match saved doRNG state and preserve caller RNG", {
  skip_if_not_installed("rngtools")
  config <- study_helpers$get_original_study_config("scenario2", "none")
  set.seed(681)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  jobs <- study_helpers$original_study_streams(config, 2L)
  expect_identical(.Random.seed, seed_before)
  expect_identical(RNGkind(), kind_before)
  expect_equal(jobs$grid$ss, rep(1:5, 2))
  expect_equal(jobs$grid$rep_id, rep(1:2, each = 5))
  expect_equal(jobs$grid$n1, rep(config$sample_size_n1, 2))
  # First job's state is a small provenance fixture from the saved study.
  saved_first <- c(10407L, -1213238679L, 1924592406L, 435078815L,
                   -1427496876L, -512481467L, -1756957886L)
  expect_identical(jobs$streams[[1L]], saved_first)
  expect_identical(jobs$streams[[2L]], parallel::nextRNGStream(saved_first))
  longer <- study_helpers$original_study_streams(config, 3L)
  expect_identical(jobs$streams, longer$streams[seq_len(10L)])
  expect_error(study_helpers$original_study_streams(config, 0), "positive integer")
})

test_that("a study job matches saved estimates and survives worker serialization", {
  skip_if_not_installed("rngtools")
  skip_if_not_installed("dbarts")
  skip_if_not_installed("ranger")
  # One dataset only: this verifies integration, not the complete simulation.
  config <- study_helpers$get_original_study_config("scenario2", "none")
  engine <- study_helpers$get_original_study_engine("scenario2")
  engine <- unserialize(serialize(engine, NULL))
  job <- study_helpers$original_study_streams(config, 1L)
  set.seed(868)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  result <- study_helpers$fit_original_study_rep(
    engine, config, 200L, 1L, job$streams[[1L]]
  )
  expect_identical(.Random.seed, seed_before)
  expect_identical(RNGkind(), kind_before)
  expect_true(all(is.na(result$Error)))
  expect_length(result$Estimate, 5L)
  expect_true(all(is.finite(result$Estimate)))
  expect_true(all(is.finite(result$SE) & result$SE > 0))
  # The numerical fixture applies to its recorded software environment.
  # Other supported versions are checked for successful, finite fits above.
  fixture_environment <- .Platform$OS.type == "windows" &&
    getRversion() == "4.5.1" &&
    utils::packageVersion("survival") == "3.8.3" &&
    utils::packageVersion("dbarts") == "0.9.32" &&
    utils::packageVersion("ranger") == "0.17.0"
  if (fixture_environment) {
    expect_equal(result$Estimate, c(rep(-0.34043567356217552, 3L),
      -0.37360440031192166, -0.39550247440702435), tolerance = 1e-8)
    expect_equal(result$SE, c(0.104636188098269953, 0.101167020439408645,
      0.074733454715128364, 0.066590198257541006, 0.075857650323297041),
      tolerance = 1e-8)
  }
  expect_error(study_helpers$fit_original_study_rep(
    engine, config, 999L, 1L, job$streams[[1L]]), "original sample sizes")
  expect_error(study_helpers$fit_original_study_rep(
    engine, config, 200L, 1L, c(1,2,3)), "L'Ecuyer-CMRG")
})
