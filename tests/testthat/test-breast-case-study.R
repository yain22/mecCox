breast_directory <- testthat::test_path("..", "..", "inst", "reproduce")
if (!file.exists(file.path(breast_directory, "breast_cancer.R"))) {
  breast_directory <- system.file("reproduce", package = "mecCox")
}
breast_helpers <- new.env(parent = baseenv())
source(file.path(breast_directory, "breast_cancer_helpers.R"), local = breast_helpers)
source(file.path(breast_directory, "simulation_helpers.R"), local = breast_helpers)
# Read the driver's pure functions without launching its full DL/RSF analysis.
for (expression in as.list(parse(file.path(breast_directory, "breast_cancer.R")))) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.call(expression[[3L]]) &&
      identical(expression[[3L]][[1L]], as.name("function"))) {
    eval(expression, envir = breast_helpers)
  }
}

testthat::test_that("public cohorts, endpoint, and unsupported category are preserved", {
  data <- breast_helpers$prepare_breast_cancer_data()
  testthat::expect_equal(nrow(data), 2889L)
  testthat::expect_equal(table(data$A), table(c(rep(1L, 246L), rep(0L, 2643L))),
                         ignore_attr = TRUE)
  testthat::expect_false(anyNA(data))
  testthat::expect_true(all(data$time > 0 & data$time <= 5))
  testthat::expect_equal(sum(data$A == 1L & data$grade == "1"), 33L)
  testthat::expect_equal(sum(data$A == 0L & data$grade == "1"), 0L)

  original <- new.env(parent = baseenv())
  utils::data("cancer", package = "survival", envir = original)
  gb <- original$gbsg[original$gbsg$hormon == 1L, ]
  testthat::expect_equal(data$time[data$A == 1L], pmin(gb$rfstime / 365.25, 5))
  testthat::expect_equal(data$delta[data$A == 1L],
                         as.integer(gb$status == 1L & gb$rfstime <= 5 * 365.25))
})

testthat::test_that("case fits support both table definitions and exact calibration", {
  testthat::skip_if_not_installed("MASS")
  data <- breast_helpers$prepare_breast_cancer_data()
  covariates <- c("age", "meno", "size_cat", "grade", "log_nodes", "log_pgr", "log_er")
  unweighted <- fit_unweighted_cox(data, "time", "delta", "A")
  ipw <- fit_att_ipw_cox(data, "time", "delta", "A", covariates)
  mec <- breast_helpers$fit_breast_mec(data, covariates,
                                      ps_learner = "glm", survival_learner = "cox")
  testthat::expect_equal(round(unweighted$theta, 3), -0.037)
  testthat::expect_equal(round(ipw$theta, 3), -0.564)
  testthat::expect_true(all(mec$weights > 0))
  testthat::expect_equal(mec$weights[data$A == 1L], rep(1, 246))
  testthat::expect_lt(abs(sum(mec$weights[data$A == 0L]) - 246), 1e-5)
  residual <- colSums(mec$basis[data$A == 0L, , drop = FALSE] *
                       mec$weights[data$A == 0L]) / 246 -
    colMeans(mec$basis[data$A == 1L, , drop = FALSE])
  testthat::expect_lt(max(abs(residual)), 1e-8)
  testthat::expect_true(all(vapply(mec$learner_diagnostics, function(fold) {
    "grade3" %in% fold$aliased_coefficients
  }, logical(1))))

  # The formatting test reuses this fit for the ML slot; a full DL/RSF run is
  # checked separately before release, without adding neural training to CI.
  fits <- list(unweighted = unweighted, ipw = ipw, glm_cox = mec, dl_rsf = mec)
  estimates <- breast_helpers$breast_cancer_estimate_table(fits)
  balance <- breast_helpers$breast_cancer_balance_table(data, fits)
  testthat::expect_equal(nrow(estimates), 6L)
  testthat::expect_equal(estimates$log_HR[2:4], rep(ipw$theta, 3L))
  testthat::expect_equal(estimates$HR, exp(estimates$log_HR))
  testthat::expect_true(all(estimates$CI_lower < estimates$HR &
                            estimates$CI_upper > estimates$HR))
  testthat::expect_equal(dim(balance), c(16L, 6L))
  grade1 <- balance[balance$Covariate == "grade 1", 3:6]
  testthat::expect_equal(as.numeric(grade1), rep(33 / 246, 4L))
  testthat::expect_equal(balance["ATT-IPW"][[1L]][15L],
                         sum(ipw$weights[data$A == 0L])^2 /
                           sum(ipw$weights[data$A == 0L]^2))
  testthat::expect_equal(round(balance[["Unweighted"]][1L], 3), 0.268)
  testthat::expect_equal(round(balance[["Unweighted"]][9L], 3), 1.158)

  testthat::skip_if_not_installed("kableExtra")
  testthat::skip_if_not_installed("htmltools")
  output <- breast_helpers$build_breast_cancer_report(estimates, balance,
                    list(seed = 20260427L, n_folds = 10L, n_landmarks = 20L))
  testthat::expect_named(output$tables, c("estimates", "balance"))
  for (table in output$tables) testthat::expect_s3_class(table, "kableExtra")
  scratch <- tempfile("breast viewer ")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  previous <- setwd(scratch)
  on.exit(setwd(previous), add = TRUE)
  viewed <- NULL
  shown <- breast_helpers$show_simulation_report(output$report, viewer = function(path) {
    testthat::expect_true(startsWith(normalizePath(path), normalizePath(tempdir())))
    viewed <<- paste(readLines(path, warn = FALSE), collapse = "\n")
  })
  testthat::expect_true(shown)
  testthat::expect_match(viewed, "Table S2", fixed = TRUE)
  testthat::expect_match(viewed, "Table S3", fixed = TRUE)
  testthat::expect_match(viewed, "Corrected sandwich", fixed = TRUE)
  testthat::expect_match(viewed, "Grade 1 occurs", fixed = TRUE)
  testthat::expect_length(list.files(scratch, all.files = TRUE, no.. = TRUE), 0L)
})
