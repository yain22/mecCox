final_helper_path <- testthat::test_path("..", "..", "inst", "reproduce",
                                        "final_study_helpers.R")
if (!file.exists(final_helper_path)) {
  final_helper_path <- system.file("reproduce", "final_study_helpers.R",
                                  package = "mecCox")
}
final_helpers <- new.env(parent = baseenv())
source(final_helper_path, local = final_helpers)

test_that("final bundle preserves the independent full grid and caller RNG", {
  skip_if_not_installed("digest")
  withr::local_seed(734)
  before <- .Random.seed
  bundle <- final_helpers$load_final_study_bundle()
  expect_identical(.Random.seed, before)
  expect_equal(nrow(bundle$grid), 30000L)
  expect_equal(length(bundle$streams), 30000L)
  expect_identical(as.integer(names(bundle$streams)), bundle$grid$id)
  expect_false(anyDuplicated(vapply(bundle$streams, paste, character(1),
                                    collapse = ",")) > 0L)
  expect_identical(as.integer(table(bundle$grid$figure)), c(15000L, 15000L))
  expect_equal(unname(table(bundle$grid$config_id, bundle$grid$n1)),
               matrix(1000L, 6L, 5L))
  for (config in bundle$configs) {
    expect_equal(config$n_landmarks, 5L)
    expect_equal(config$Kfold_mec, 10L)
    expect_true(is.na(config$theta_true))
  }
  severe <- bundle$configs[["figure6_severe"]]
  expect_identical(c(severe$ps_nonlinearity, severe$or_nonlinearity), c(1.5, 4))
  expect_equal(severe$n0_over_n1, 4)
  none_targets <- subset(bundle$targets,
                         figure == "Figure5" | setting == "none")
  expect_length(unique(none_targets$numerical_target), 1L)
  expect_true(all(is.finite(bundle$targets$mean_MCSE) &
                    bundle$targets$mean_MCSE > 0))
})

test_that("final study BART preserves mean-first conversion without a second offset", {
  skip_if_not_installed("digest")
  skip_if_not_installed("dbarts")
  skip_if_not_installed("ranger")
  training <- data.frame(A = c(0L, 0L, 0L, 1L), X1 = 1:4)
  validation <- training[1:2, ]
  draws <- matrix(c(-2, 0.4, -0.8, 2.4), nrow = 2L)
  calls <- list()
  local_mocked_bindings(bart = function(...) {
    calls[[length(calls) + 1L]] <<- list(...)
    list(yhat.test = draws)
  }, .package = "dbarts")
  expected <- pnorm(colMeans(draws))
  expect_gt(max(abs(expected - colMeans(pnorm(draws)))), 0.01)
  expect_gt(max(abs(expected - pnorm(colMeans(draws) + qnorm(0.25)))), 0.1)
  withr::local_seed(918)
  for (figure in c("figure5", "figure6")) {
    engine <- final_helpers$get_final_study_engine(figure)
    actual <- engine$estimate_ps_ml(
      training, "X1", ps_method = "bart", newdata = validation,
      auto_tune = TRUE, tune_seed = 918L, trim = c(0.01, 0.99)
    )
    expect_equal(actual$ehat, expected, tolerance = 0)
    args <- calls[[length(calls)]]
    expect_equal(args$binaryOffset, qnorm(0.25))
    expect_identical(args$ntree, 100L)
    expect_identical(args$ndpost, 300L)
    expect_identical(args$nskip, 150L)
    expect_equal(args$k, 2)
    expect_equal(args$nthread, 1)
  }
})

test_that("final study output-range and probability-field rules are retained", {
  skip_if_not_installed("digest")
  skip_if_not_installed("dbarts")
  skip_if_not_installed("ranger")
  training <- data.frame(A = c(0L, 0L, 0L, 1L), X1 = 1:4)
  validation <- training[1:2, ]
  returned <- list(yhat.test.mean = c(0.2, 0.8))
  local_mocked_bindings(bart = function(...) returned, .package = "dbarts")
  engine <- final_helpers$get_final_study_engine("figure6")
  predict <- function() engine$estimate_ps_ml(
    training, "X1", ps_method = "bart", newdata = validation,
    auto_tune = FALSE, trim = c(0.01, 0.99)
  )$ehat
  expect_equal(predict(), c(0.2, 0.8))
  returned <- list(yhat.test.mean = c(-0.8, 0.8))
  expect_equal(predict(), pnorm(c(-0.8, 0.8)))
  returned <- list(yhat.test.mean = c(-0.6, 0.4),
                   yhat.test = matrix(c(-3, 1, 2, 4), nrow = 2L))
  expect_equal(predict(), pnorm(c(-0.6, 0.4)))
  returned <- list(prob.test.mean = c(0.3, 0.7), yhat.test.mean = c(-0.8, 0.8))
  expect_equal(predict(), c(0.3, 0.7))
})

test_that("final KL cross-derivative includes weighted events and risk sets", {
  # Deterministic tied/censored data; no nuisance model or simulation fit.
  data <- data.frame(A = c(1, 0, 0, 1, 0, 1, 0, 0),
                     time = c(2, 2, 4, 5, 5, 7, 8, 9),
                     delta = c(1, 1, 0, 1, 1, 0, 1, 0))
  controls <- which(data$A == 0)
  H <- cbind(1, seq(-1, 1, length.out = 5L), c(0.3, 0.7, -0.2, 0.5, -0.4))
  lambda <- c(0.12, -0.07, 0.09)
  baseline <- c(0.7, 1.2, 0.9, 1.6, 1.1)
  weights <- function(value) {
    answer <- rep(1, nrow(data))
    answer[controls] <- baseline * exp(as.vector(H %*% value))
    answer
  }
  theta <- -0.3
  # Independently evaluate the Breslow score by event time. In particular,
  # each perturbed weight enters both event totals and the risk-set mean.
  score <- function(value) {
    W <- weights(value)
    sum(vapply(sort(unique(data$time[data$delta == 1])), function(t) {
      event <- data$time == t & data$delta == 1
      risk <- data$time >= t
      risk_weight <- W[risk] * exp(theta * data$A[risk])
      average <- sum(risk_weight * data$A[risk]) / sum(risk_weight)
      sum(W[event] * data$A[event]) - sum(W[event]) * average
    }, numeric(1)))
  }
  expect_true(all(is.finite(weights(lambda)) & weights(lambda) > 0))
  directions <- cbind(diag(3L), c(0.4, -0.7, 0.2))
  h <- 1e-5
  numerical <- apply(directions, 2L, function(v) {
    (score(lambda + h * v) - score(lambda - h * v)) / (2 * h)
  })
  for (figure in c("figure5", "figure6")) {
    engine <- new.env(parent = asNamespace("stats"))
    sys.source(file.path(dirname(final_helper_path), "final_reference",
                         paste0(figure, "_functions.R")), engine)
    ce <- engine$cox_score_eta_fixed_weight(theta, data, weights(lambda))
    expect_equal(ce$U, score(lambda), tolerance = 1e-12)
    assignments <- Filter(function(x) {
      is.call(x) && identical(x[[1L]], as.name("<-")) &&
        identical(x[[2L]], as.name("D_theta_lambda"))
    }, as.list(body(engine$fit_mec_cox)))
    expect_length(assignments, 1L)
    # Evaluate the actual frozen fitting expression, not a duplicated formula.
    b <- eval(assignments[[1L]][[3L]], list(
      divergence = "kl", eta_theta = ce$eta, idx0 = controls, H0 = H
    ))
    expect_equal(as.vector(b %*% directions), numerical, tolerance = 2e-8)
  }
})
