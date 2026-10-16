# Supplemental reproduction contracts and caller-state regression tests.
# These tests do not fit a model, calculate a target, or write a figure.

supp234_test_path <- function(name) {
  p <- testthat::test_path("..", "..", "inst", "reproduce", name)
  if (!file.exists(p)) p <- system.file("reproduce", name, package = "mecCox")
  stopifnot(nzchar(p), file.exists(p))
  p
}

supp234_test_api <- function() {
  e <- new.env(parent = globalenv())
  sys.source(supp234_test_path("supplement_S2_S4_helpers.R"), e)
  e
}

supp234_test_bundle <- function(api) {
  api$load_supp234_bundle(supp234_test_path("final_supp_S2_S4_reference"))
}

test_that("supplemental files only define functions when sourced", {
  withr::local_seed(421)
  paths <- vapply(c("supplement_S2_S4_helpers.R", "run_supplement_S2_S4.R",
    "redraw_supplement_S2_S4.R"), function(name)
      normalizePath(supp234_test_path(name), mustWork = TRUE), character(1))
  folder <- withr::local_tempdir()
  withr::local_dir(folder)
  snapshot <- function() list(seed = .Random.seed, kind = RNGkind(), wd = getwd(),
    locale = Sys.getlocale(), options = options(), env = Sys.getenv(), search = search(),
    namespaces = sort(loadedNamespaces()), files = list.files(folder, recursive = TRUE, all.files = TRUE))
  for (path in paths) {
    ast <- as.list(parse(path, keep.source = FALSE))
    expect_true(all(vapply(ast, function(x) is.call(x) &&
      identical(x[[1L]], as.name("<-")) && is.call(x[[3L]]) &&
      identical(x[[3L]][[1L]], as.name("function")), logical(1))))
    e <- new.env(parent = globalenv())
    before <- snapshot()
    sys.source(path, e)
    expect_identical(snapshot(), before)
    expect_true(length(ls(e)) > 0L)
    expect_true(all(vapply(mget(ls(e), e), is.function, logical(1))))
  }
})

test_that("supplemental bundle retains six designs and independent saved streams", {
  skip_if_not_installed("digest")
  api <- supp234_test_api()
  withr::local_seed(424)
  before <- .Random.seed
  b <- supp234_test_bundle(api)
  expect_identical(.Random.seed, before)
  expect_identical(names(b$configs), c("S2_p50", "S2_p200", "S3_ratio2", "S3_ratio4", "S4_ratio2", "S4_ratio4"))
  expect_identical(b$grid$id, seq_len(30000L))
  expect_identical(names(b$streams), as.character(b$grid$id))
  cell_counts <- table(b$grid$config_id, b$grid$n1)
  expect_identical(dim(cell_counts), c(6L, 5L))
  expect_true(all(cell_counts == 1000L))
  expect_true(all(vapply(b$streams, function(x) is.integer(x) && length(x) == 7L &&
    !anyNA(x) && x[1L] %% 10000L == 407L, logical(1))))
  expect_equal(anyDuplicated(vapply(b$streams, paste, character(1), collapse = ",")), 0L)
  expect_true(all(vapply(b$configs, function(x) x$Kfold_mec == 5 && is.na(x$theta_true), logical(1))))
  expect_equal(vapply(b$configs, `[[`, numeric(1), "p"),
    c(S2_p50 = 50, S2_p200 = 200, S3_ratio2 = 50, S3_ratio4 = 50, S4_ratio2 = 50, S4_ratio4 = 50))
  for (k in c("S2_p50", "S2_p200")) {
    expect_identical(b$configs[[k]]$mec_ps_method, "glm_lasso")
    expect_identical(b$configs[[k]]$mec_or_method, "lasso_cox_lp")
    expect_equal(b$configs[[k]]$n0_over_n1, 4)
  }
  for (k in c("S3_ratio2", "S3_ratio4", "S4_ratio2", "S4_ratio4")) {
    expect_identical(b$configs[[k]]$mec_ps_method, "glm")
    expect_identical(b$configs[[k]]$mec_or_method, "cox")
    expect_equal(b$configs[[k]]$n_landmarks, 5)
    expect_equal(b$configs[[k]]$n0_over_n1, if (endsWith(k, "2")) 2 else 4)
    expect_equal(b$configs[[k]]$censor_rate, if (startsWith(k, "S4")) .0016 else .0008)
  }
  expect_identical(b$targets$config_id, names(b$configs))
  expect_true(all(b$targets$attempted == 100 & b$targets$n1_super == 150000 & b$targets$n0_super == 300000))
  expect_length(unique(b$targets$numerical_target[1:4]), 1L)
  expect_length(unique(b$targets$numerical_target[5:6]), 1L)
  expect_false(identical(b$targets$numerical_target[1L], b$targets$numerical_target[5L]))
  expected_hash <- c(runtime_configs.rds = "b858d09f772bf8cbb5a3f7ce5340ea0b62b823d299878d6341dc2a2dcebabd7e",
    job_grid.rds = "d5bdfdb4ca4b22b768482bb54e414996e7db49246088a5f1ade88dd3918e602a",
    rng_states.rds = "f0728ba43a72c48cfbe7234cacfbc85da7375bc5b774d8e1bae20d872b533764",
    targets.rds = "b51c9504f7a822d43055b6f22f4ad1dfc9dcd9f699bba0dc111ffd6b17e2c4ec")
  for (n in names(expected_hash)) expect_identical(api$supp234_sha(file.path(b$directory, n)), unname(expected_hash[n]))
})

test_that("corrupt reference bytes and traversal entries are rejected", {
  skip_if_not_installed("digest")
  api <- supp234_test_api(); b <- supp234_test_bundle(api)
  copy <- withr::local_tempdir()
  files <- list.files(b$directory, full.names = TRUE)
  expect_true(all(file.copy(files, copy)))
  path <- file.path(copy, "runtime_configs.rds")
  con <- file(path, "ab"); writeBin(as.raw(0L), con); close(con)
  expect_error(api$load_supp234_bundle(copy))
  malicious <- withr::local_tempdir()
  utils::write.csv(data.frame(file = "../outside.rds", sha256 = paste(rep("0", 64), collapse = "")),
    file.path(malicious, "SHA256.csv"), row.names = FALSE)
  expect_error(api$supp234_check_hashes(malicious))
})

test_that("supplemental engine namespaces preserve family-specific interfaces", {
  for (p in c("digest", "MASS", "dplyr", "survival")) skip_if_not_installed(p)
  api <- supp234_test_api(); b <- supp234_test_bundle(api)
  before <- search()
  s2 <- api$supp234_plain_engine("S2", b)
  s34 <- api$supp234_plain_engine("S34", b)
  expect_identical(search(), before)
  expect_true(exists("fit_mec_cox_from_nuisance", s2, inherits = FALSE))
  expect_false(exists("fit_mec_cox_from_nuisance", s34, inherits = FALSE))
  expect_true("ps_lasso_lambda_rule" %in% names(formals(s2$fit_one_rep)))
  expect_false("ps_lasso_lambda_rule" %in% names(formals(s34$fit_one_rep)))
  expect_identical(get("filter", s2), dplyr::filter)
  expect_identical(get("coxph", s2), survival::coxph)
})

test_that("original eligibility retains zero SE and invalid CI diagnostics", {
  for (p in c("digest", "MASS", "dplyr", "survival")) skip_if_not_installed(p)
  api <- supp234_test_api(); b <- supp234_test_bundle(api)
  rows <- data.frame(config_id = "synthetic", figure = "S2", n1 = 200, n0 = 800,
    Method = "MEC-Cox synthetic", theta_true = 0, Estimate = c(0, 0, 0, 0, NA_real_, 0),
    SE = c(0, .1, .1, .1, .1, Inf), CI_L = c(-1, 1, NA, -Inf, -1, -1),
    CI_U = c(1, -1, .1, Inf, 1, 1), ESS = 100, Rel_ESS = .5, W_CV = .2,
    W_Min = .1, W_Max = 2, Cal_Grad = 0, Cal_Converged = TRUE)
  counts <- api$supp234_original_group_counts(rows)
  expect_equal(counts$Attempted, 6L)
  expect_equal(counts$Trials, 4L)
  for (family in c("S2", "S34")) {
    e <- api$supp234_plain_engine(family, b)
    s <- e$summarize_simulation(rows)
    expect_equal(s$Trials, 4L)
    expect_true(is.na(s$Coverage_95))
    expect_equal(s$Bias, 0)
    expect_equal(s$RMSE, 0)
    zero <- rows; zero$Estimate[] <- NA_real_; zero$SE[] <- Inf
    expect_equal(nrow(e$summarize_simulation(zero)), 0L)
    expect_equal(api$supp234_original_group_counts(zero)$Attempted, 6L)
    expect_equal(api$supp234_original_group_counts(zero)$Trials, 0L)
  }
})

test_that("public fitting wrapper restores caller state on mocked success and error", {
  skip_if_not_installed("digest"); skip_if_not_installed("glmnet")
  api <- supp234_test_api(); b <- supp234_test_bundle(api)
  # Exercise the actual wrapper's restoration without running scientific code.
  api$supp234_check_fit_environment <- function(bundle) invisible(TRUE)
  api$supp234_context <- function() new.env(parent = globalenv())
  api$sys.source <- function(file, envir, keep.source = FALSE) {
    envir$study_one <- function(id, m, directory, derivative, logging, capture_full) {
      stopifnot(identical(derivative, "analytic"), isTRUE(logging), isTRUE(capture_full),
        getOption("ranger.num.threads") == 1L, Sys.getenv("OMP_NUM_THREADS") == "1")
      set.seed(12)
      if (fail) stop("deliberate no-model fixture error")
      list(id = id, fixture_only = TRUE)
    }
    environment(envir$study_one) <- api
  }
  withr::local_seed(902)
  withr::local_envvar(c(OMP_NUM_THREADS = "7", R_RANGER_NUM_THREADS = NA_character_))
  withr::local_options(list(ranger.num.threads = 3L))
  snapshot <- function() list(seed = get0(".Random.seed", .GlobalEnv), kind = RNGkind(),
    env = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
      "RCPP_PARALLEL_NUM_THREADS", "R_RANGER_NUM_THREADS"), unset = NA_character_),
    option = getOption("ranger.num.threads"), locale = Sys.getlocale(), cwd = getwd())
  for (seed_present in c(TRUE, FALSE)) {
    if (!seed_present) rm(".Random.seed", envir = .GlobalEnv)
    before <- snapshot()
    api$fail <- FALSE
    expect_true(api$fit_supp234_job(1L, b, capture_full = TRUE)$fixture_only)
    expect_identical(snapshot(), before)
    api$fail <- TRUE
    expect_error(api$fit_supp234_job(1L, b, capture_full = TRUE), "deliberate no-model")
    expect_identical(snapshot(), before)
  }
})

test_that("plotting locale is restored on success and a deliberate error", {
  api <- supp234_test_api()
  withr::local_seed(903)
  previous <- Sys.getlocale("LC_CTYPE")
  withr::defer(Sys.setlocale("LC_CTYPE", previous))
  suppressWarnings(Sys.setlocale("LC_CTYPE", "C"))
  before <- list(locale = Sys.getlocale(), seed = .Random.seed)
  probe <- tryCatch(api$supp234_plot_locale(), error = identity)
  if (inherits(probe, "error")) {
    expect_identical(list(locale = Sys.getlocale(), seed = .Random.seed), before)
    skip("The platform has no supported UTF-8 LC_CTYPE")
  }
  expect_true(isTRUE(l10n_info()[["UTF-8"]]))
  api$supp234_restore_plot_locale(probe)
  expect_identical(list(locale = Sys.getlocale(), seed = .Random.seed), before)
  # This fails after entering the scoped UTF-8 locale but before package/style
  # loading or a graphics device. The real on.exit restoration must still run.
  expect_error(api$supp234_plot_style("deliberately_invalid", list()))
  expect_identical(list(locale = Sys.getlocale(), seed = .Random.seed), before)
})
