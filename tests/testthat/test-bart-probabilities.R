test_that("BART propensity averages probabilities with a nonzero offset", {
  skip_if_not_installed("dbarts")
  set.seed(721)
  training <- data.frame(A = c(rep(0L, 80L), rep(1L, 20L)),
                         X1 = rnorm(100), X2 = rnorm(100))
  validation <- training[c(1:8, 81:88), ]
  xvars <- c("X1", "X2")
  offset <- qnorm(mean(training$A))
  seed <- 912L
  posterior <- withr::with_seed(seed, dbarts::bart(
    x.train = as.matrix(training[, xvars]), y.train = training$A,
    x.test = as.matrix(validation[, xvars]), ntree = 20L,
    ndpost = 40L, nskip = 20L, k = 2, binaryOffset = offset,
    verbose = FALSE, keeptrees = FALSE, keeptrainfits = FALSE,
    nthread = 1L, seed = seed
  ))
  expected <- as.numeric(colMeans(pnorm(posterior$yhat.test)))
  actual <- mecCox:::.mec_propensity_bart(
    training, validation, xvars, seed = seed, num_trees = 20L,
    posterior_draws = 40L, burn_in = 20L, shrinkage = 2
  )
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_true(all(is.finite(actual) & actual > 0 & actual < 1))
  # These data distinguish both an extra offset and a CDF of latent means.
  expect_gt(max(abs(expected - colMeans(pnorm(
    posterior$yhat.test + offset
  )))), 0.05)
  expect_gt(max(abs(expected - pnorm(colMeans(posterior$yhat.test)))), 1e-5)
})

test_that("pre-final compatibility snapshots retain their probability means", {
  skip_if_not_installed("dbarts")
  helper_path <- testthat::test_path("..", "..", "inst", "reproduce",
                                    "original_study_helpers.R")
  if (!file.exists(helper_path)) {
    helper_path <- system.file("reproduce", "original_study_helpers.R",
                               package = "mecCox")
  }
  helpers <- new.env(parent = baseenv())
  source(helper_path, local = helpers)
  set.seed(721)
  training <- data.frame(A = c(rep(0L, 80L), rep(1L, 20L)),
                         X1 = rnorm(100), X2 = rnorm(100))
  validation <- training[c(1:8, 81:88), ]
  xvars <- c("X1", "X2")
  posterior <- withr::with_seed(912, dbarts::bart(
    x.train = as.matrix(training[, xvars]), y.train = training$A,
    x.test = as.matrix(validation[, xvars]), ntree = 20L,
    ndpost = 40L, nskip = 20L, binaryOffset = qnorm(mean(training$A)),
    verbose = FALSE, keeptrees = FALSE, nthread = 1
  ))
  expected <- pmin(pmax(as.numeric(colMeans(pnorm(posterior$yhat.test))),
                        0.1), 0.9)
  for (scenario in c("scenario1", "scenario2")) {
    engine <- helpers$get_original_study_engine(scenario)
    actual <- withr::with_seed(912, engine$estimate_ps_ml(
      training, xvars, ps_method = "bart", newdata = validation,
      bart_ntree = 20L, bart_ndpost = 40L, bart_nskip = 20L,
      trim = c(0.1, 0.9), auto_tune = FALSE
    ))
    expect_equal(actual$ehat, expected, tolerance = 1e-12)
  }
})

test_that("API and pre-final snapshots handle mean-only BART output", {
  skip_if_not_installed("dbarts")
  latent_means <- c(-0.3, 0.4)
  local_mocked_bindings(
    bart = function(...) list(yhat.test.mean = latent_means),
    .package = "dbarts"
  )
  training <- data.frame(A = c(0L, 0L, 0L, 1L), X1 = 1:4)
  validation <- training[1:2, ]
  expected <- pnorm(latent_means)
  actual <- mecCox:::.mec_propensity_bart(
    training, validation, "X1", seed = 912L, num_trees = 20L,
    posterior_draws = 40L, burn_in = 20L, shrinkage = 2
  )
  expect_equal(actual, expected)
  helper_path <- testthat::test_path("..", "..", "inst", "reproduce",
                                    "original_study_helpers.R")
  if (!file.exists(helper_path)) {
    helper_path <- system.file("reproduce", "original_study_helpers.R",
                               package = "mecCox")
  }
  helpers <- new.env(parent = baseenv())
  source(helper_path, local = helpers)
  for (scenario in c("scenario1", "scenario2")) {
    engine <- helpers$get_original_study_engine(scenario)
    actual <- engine$estimate_ps_ml(
      training, "X1", ps_method = "bart", newdata = validation,
      auto_tune = FALSE
    )
    expect_equal(actual$ehat, expected)
  }
})
