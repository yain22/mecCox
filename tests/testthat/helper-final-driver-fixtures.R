# Cheap final-API stand-ins: these fixtures never call a scientific fit.
final_driver_test_directory <- function() {
  directory <- testthat::test_path("..", "..", "inst", "reproduce")
  if (!dir.exists(directory)) directory <- system.file("reproduce", package = "mecCox")
  normalizePath(directory, mustWork = TRUE)
}

final_driver_test_fixture <- function(scenario) {
  directory <- final_driver_test_directory()
  execution <- new.env(parent = baseenv())
  source(file.path(directory, "simulation_helpers.R"), local = execution)
  source(file.path(directory, "final_driver_helpers.R"), local = execution)
  expressions <- as.list(parse(file.path(directory, paste0(scenario, ".R"))))
  for (expression in expressions) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        is.call(expression[[3L]]) &&
        identical(expression[[3L]][[1L]], as.name("function"))) eval(expression, execution)
  }
  reference <- execution$load_final_driver_reference(directory)
  reference$api$fit_final_study_job <- function(job_id, bundle, capture_full = FALSE) {
    job <- bundle$grid[bundle$grid$id == job_id, , drop = FALSE]
    had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
    old_kind <- RNGkind()
    on.exit({
      do.call(RNGkind, as.list(old_kind))
      if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv)
      else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    assign(".Random.seed", bundle$streams[[as.character(job_id)]], .GlobalEnv)
    estimate <- stats::rnorm(1L)
    rows <- data.frame(
      rep = job$rep_id, n1 = job$n1, n0 = job$n0,
      Method = bundle$method_labels[[job$figure]], theta_true = NA_real_,
      Estimate = estimate, SE = 0.1, CI_L = estimate - 0.196,
      CI_U = estimate + 0.196, ESS = 100, Rel_ESS = 0.5,
      W_CV = 1, W_Min = 0.1, W_Max = 2, Cal_Grad = 1e-8,
      Cal_Converged = TRUE, Error = NA_character_, stringsAsFactors = FALSE
    )
    ors <- if (job$figure == "figure5") "cox" else c("cox", "rsf")
    list(job = job, results = rows, fatal_error = NA_character_,
         warnings = "retained fixture warning", data_hash = paste0("fixture-", job$id),
         PS_trace = list(list(ids = 1:3, p = c(.1, .2, .3))),
         nuisance_trace = list(list(H_hash = "fixture-H")),
         calibration_trace = list(list(converged = TRUE)),
         estimator_trace = lapply(ors, function(or) list(
           or_method = or, raw_sandwich_SE = .1, fallback_used = FALSE)),
         math_trace = lapply(ors, function(or) list(
           or_method = or, positive_finite_weights = TRUE,
           exponential_weights_identical = TRUE, solve_path = "direct")),
         solve_trace = list(list(path = "direct")),
         OR_tuning_trace = list(list(rule = "unchanged")),
         R_rng_after = get(".Random.seed", .GlobalEnv, inherits = FALSE))
  }
  environment(reference$api$fit_final_study_job) <- reference$api
  execution$final_reference <- reference
  key <- if (scenario == "scenario1") "2" else "none"
  config_id <- if (scenario == "scenario1") "figure5_ratio2" else "figure6_none"
  config <- reference$bundle$configs[[config_id]]
  execution$study_configs <- stats::setNames(list(config), key)
  execution$study_reference <- reference$bundle$provenance
  execution$design <- list(
    replications = 3L, treated_sizes = config$sample_size_grid$n1[1:2],
    control_multipliers = 2L, control_multiplier = config$n0_over_n1
  )
  execution$settings <- data.frame(setting_id = 1L, setting_key = "none",
                                    setting = "none", kappa_pi = 0, kappa_m = 0)
  execution$options <- list(quick_run = TRUE, cores = 1L, replications = 3L)
  execution
}
