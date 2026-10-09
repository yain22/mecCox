# Run the second simulation using the bundled main-paper study engine.
# Run after installing mecCox, dbarts, and ranger:
#   Rscript path/to/mecCox/inst/reproduce/scenario2.R
# The default is up to 20 workers; --cores=1 runs sequentially.
# --quick exercises all three settings with a deliberately reduced workload.
# In the R console or RStudio, edit the three settings below, then source this
# file (or run it from the editor). Rscript arguments override these settings.
# Results stay in R, with formatted tables in the RStudio Viewer and plots in
# the active graphics window. The script does not save result files.

quick_run <- FALSE
cores <- 20L
replications <- 1000L

# source() records the current file in an `ofile` frame. Selected lines in an
# editor have no such frame, so also look in the working directory and package.
source_files <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
}, character(1))
source_files <- rev(source_files[nzchar(source_files)])
script_option <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                      value = TRUE)
script_files <- c(source_files, sub("^--file=", "", script_option))
helper_candidates <- c(
  file.path(dirname(script_files), "simulation_helpers.R"),
  "simulation_helpers.R",
  file.path("inst", "reproduce", "simulation_helpers.R"),
  file.path("mecCox", "inst", "reproduce", "simulation_helpers.R"),
  system.file("reproduce", "simulation_helpers.R", package = "mecCox")
)
helper_candidates <- unique(helper_candidates[
  nzchar(helper_candidates) & file.exists(helper_candidates)
])
helper_candidates <- helper_candidates[vapply(dirname(helper_candidates), function(path) {
  file.exists(file.path(path, "original_study_helpers.R")) &&
    dir.exists(file.path(path, "study_reference"))
}, logical(1))]
if (!length(helper_candidates)) {
  stop("Cannot find the complete simulation helpers. Install the current mecCox ",
       "package, or keep simulation_helpers.R, original_study_helpers.R, and ",
       "study_reference beside this script.", call. = FALSE)
}
source(helper_candidates[1L], local = TRUE)
source(file.path(dirname(helper_candidates[1L]), "original_study_helpers.R"), local = TRUE)
if (!exists("build_simulation_report", mode = "function", inherits = FALSE)) {
  stop("Update mecCox or download the complete current reproduce folder.", call. = FALSE)
}

# A sourced file must not interpret arguments belonging to an outer R script.
arguments <- if (interactive() || length(source_files)) {
  character()
} else {
  commandArgs(trailingOnly = TRUE)
}
options <- parse_simulation_arguments(
  arguments, quick_run = quick_run, cores = cores, replications = replications
)
check_simulation_display_packages()

required_packages <- c("mecCox", "survival", "MASS", "rngtools", "dbarts", "ranger")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Install these packages before starting: ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}

# Bundled configurations supply the main-paper study targets and fitting settings.
study_engine <- get_original_study_engine("scenario2")
study_reference <- original_study_metadata("scenario2")
study_keys <- c("none", "mild", "severe")
study_configs <- stats::setNames(lapply(study_keys, function(key) {
  get_original_study_config("scenario2", key)
}), study_keys)
# Create the complete study grid before selecting quick-check jobs.
# Its order is sample size first within each replicate, as in expand.grid(ss, rep_id).
study_streams <- lapply(study_configs, original_study_streams,
                        replicates = options$replications)
study_config <- study_configs[[1L]]
design <- list(
  seed = study_config$seed,
  replications = options$replications,
  treated_sizes = study_config$sample_size_n1,
  control_multiplier = study_config$n0_over_n1,
  covariate_count = study_config$p,
  folds = study_config$Kfold_mec,
  landmarks = study_config$n_landmarks,
  weibull_scale = study_config$lambda0,
  weibull_shape = study_config$eta_shape,
  censoring_rate = study_config$censor_rate,
  conditional_log_hr = study_config$beta_cond,
  engine = "Bundled main-paper study engine",
  engine_source = study_reference$source_basename,
  engine_source_sha256 = study_reference$source_sha256,
  engine_snapshot_sha256 = study_reference$engine_sha256
)
settings <- data.frame(
  setting_id = seq_along(study_keys), setting_key = study_keys,
  setting = c("No nonlinearity", "Mild nonlinearity", "Severe nonlinearity"),
  kappa_pi = vapply(study_configs, function(config) config$ps_nonlinearity, numeric(1)),
  kappa_m = vapply(study_configs, function(config) config$or_nonlinearity, numeric(1)),
  stringsAsFactors = FALSE, row.names = NULL
)

if (options$quick_run) {
  design$treated_sizes <- study_config$sample_size_n1[1L]

  message("Quick check: fewer datasets and sample-size cells; study learner ",
          "settings, folds, landmarks, and reference targets are retained.")
}

make_result <- function(setting, treated_count, replicate, method, target,
                        design, estimate = NA_real_,
                        standard_error = NA_real_, error = NA_character_) {
  data.frame(
    setting = setting$setting,
    kappa_pi = setting$kappa_pi,
    kappa_m = setting$kappa_m,
    n1 = treated_count,
    n0 = design$control_multiplier * treated_count,
    replicate = replicate,
    method = method,
    target = target,
    estimate = estimate,
    standard_error = standard_error,
    ci_lower = estimate - 1.96 * standard_error,
    ci_upper = estimate + 1.96 * standard_error,
    error = error,
    stringsAsFactors = FALSE
  )
}

run_replication <- function(setting, treated_count, replicate, target, design) {
  key <- as.character(setting$setting_key)
  plan <- study_streams[[key]]
  job <- which(plan$grid$n1 == treated_count & plan$grid$rep_id == replicate)
  if (length(job) != 1L) stop("The requested run is not in the study grid.")
  raw <- fit_original_study_rep(study_engine, study_configs[[key]],
                                n1 = treated_count, rep_id = replicate,
                                rng_stream = plan$streams[[job]])
  methods <- rep(NA_character_, nrow(raw))
  methods[raw$Method == "ATT-IPW Cox: Naive model-based"] <- "Naive"
  methods[raw$Method == "ATT-IPW Cox: Lin-Wei"] <- "Robust sandwich"
  methods[raw$Method == "ATT-IPW Cox: Shu"] <- "Corrected sandwich"
  methods[grepl("^MEC-Cox:", raw$Method) & grepl("OR=cox", raw$Method, fixed = TRUE)] <- "MEC-Cox (BART/Cox)"
  methods[grepl("^MEC-Cox:", raw$Method) & grepl("OR=rsf", raw$Method, fixed = TRUE)] <- "MEC-Cox (BART/RSF)"
  keep <- !is.na(methods)
  raw <- raw[keep, , drop = FALSE]
  methods <- methods[keep]
  if (!nrow(raw)) stop("The study engine returned no requested method rows.")
  answer <- make_result(setting, treated_count, raw$rep, methods,
                        raw$theta_true, design, raw$Estimate, raw$SE, raw$Error)
  # Keep the study engine's confidence intervals and fit diagnostics.
  answer$ci_lower <- raw$CI_L
  answer$ci_upper <- raw$CI_U
  answer$original_method <- raw$Method
  answer$original_rng_job <- plan$grid$job_id[job]
  for (name in c("ESS", "Rel_ESS", "W_CV", "W_Min", "W_Max", "Cal_Grad", "Cal_Converged")) {
    answer[[name]] <- raw[[name]]
  }
  answer
}

summarize_scenario2_results <- function(results, settings) {
  groups <- split(results,
                  interaction(results$setting, results$n1, results$method,
                              drop = TRUE))
  summaries <- lapply(groups, function(group) {
    valid <- is.finite(group$estimate) &
      is.finite(group$standard_error)
    successful <- group[valid, , drop = FALSE]
    error <- successful$estimate - successful$target
    result <- group[1L, c("setting", "kappa_pi", "kappa_m", "n1", "n0",
                         "method")]
    result$replications <- nrow(group)
    result$successful <- nrow(successful)
    result$failed <- nrow(group) - nrow(successful)
    result$coverage <- if (length(error)) mean(
      successful$target >= successful$ci_lower & successful$target <= successful$ci_upper
    ) else NA_real_
    result$bias <- if (length(error)) mean(error) else NA_real_
    result$rmse <- if (length(error)) sqrt(mean(error^2)) else NA_real_
    result
  })
  answer <- do.call(rbind, summaries)
  rownames(answer) <- NULL
  methods <- c("Naive", "Robust sandwich", "Corrected sandwich",
               "MEC-Cox (BART/Cox)", "MEC-Cox (BART/RSF)")
  answer[order(match(answer$setting, settings$setting), answer$n1,
               match(answer$method, methods)), ]
}

run_scenario2 <- function(design, settings, options) {
  started <- proc.time()[["elapsed"]]
  worker_count <- choose_worker_count(options$cores, design$replications)
  planned_cells <- nrow(settings) * length(design$treated_sizes)
  planned_replications <- planned_cells * design$replications
  state <- new.env(parent = emptyenv())
  state$rows <- vector("list", planned_replications)
  state$completed <- 0L
  state$status <- "completed"
  state$condition <- NULL
  cluster <- NULL
  on.exit(stop_simulation_cluster(cluster), add = TRUE)
  export_names <- c("make_result", "run_replication", "fit_original_study_rep",
                    "study_engine", "study_configs", "study_streams")
  execution <- list(
    requested_cores = options$cores,
    workers = worker_count,
    backend = if (worker_count == 1L) "serial" else "PSOCK",
    rng_kind = RNGkind()
  )
  message(sprintf("Scenario 2 uses %d worker(s); %d requested.",
                  worker_count, options$cores))
  message(sprintf(
    "Workload: %d settings x %d sample sizes = %d cells; %d runs per cell = %d datasets in total.",
    nrow(settings), length(design$treated_sizes), planned_cells,
    design$replications, planned_replications
  ))
  targets <- settings
  targets$target <- vapply(study_configs[as.character(settings$setting_key)],
                            function(config) config$theta_true, numeric(1))
  tryCatch({
    cluster <- start_simulation_cluster(worker_count, export_names,
                                        envir = environment(run_replication))
    for (index in seq_len(nrow(settings))) {
      message(sprintf("%s: reference log-hazard ratio %.6f",
                      settings$setting[index], targets$target[index]))
    }

    cell_index <- 0L
    for (setting_index in seq_len(nrow(settings))) {
      setting <- settings[setting_index, , drop = FALSE]
      target <- targets$target[setting_index]
      for (treated_count in design$treated_sizes) {
        cell_index <- cell_index + 1L
        cell_completed <- 0L
        last_progress <- proc.time()[["elapsed"]]
        message(sprintf("Running cell %d/%d: %s, n1=%d, n0=%d (%d runs)",
                        cell_index, planned_cells, setting$setting,
                        treated_count, design$control_multiplier * treated_count,
                        design$replications))
        record_result <- function(replicate, value) {
          # Store by task ID, not completion order, to keep serial and parallel
          # output identical. Only this short state update defers interrupts.
          suspendInterrupts({
            task <- (cell_index - 1L) * design$replications + replicate
            state$rows[task] <- list(value)
            state$completed <- state$completed + 1L
            cell_completed <<- cell_completed + 1L
          })
          now <- proc.time()[["elapsed"]]
          if (cell_completed == 1L || cell_completed == design$replications ||
              now - last_progress >= 5) {
            message(sprintf(
              "Completed %d/%d in cell %d/%d; %d/%d datasets overall; elapsed %.1f min (replicate %d).",
              cell_completed, design$replications, cell_index, planned_cells,
              state$completed, planned_replications, (now - started) / 60,
              replicate
            ))
            last_progress <<- now
          }
        }
        run_simulation_replications(
          design$replications, run_replication,
          arguments = list(setting = setting, treated_count = treated_count,
                           target = target, design = design),
          cluster = cluster, on_result = record_result
        )
      }
    }
  }, interrupt = function(condition) {
    state$status <- "interrupted"
    state$condition <- conditionMessage(condition)
  }, error = function(condition) {
    state$status <- "failed"
    state$condition <- conditionMessage(condition)
  })

  completed_rows <- Filter(Negate(is.null), state$rows)
  if (length(completed_rows)) {
    results <- do.call(rbind, completed_rows)
    rownames(results) <- NULL
    summary <- summarize_scenario2_results(results, settings)
    print(summary, row.names = FALSE, digits = 4)
    if (any(summary$failed > 0L)) {
      message("Some fits failed; inspect the error column in the individual simulation results.")
    }
  } else {
    results <- make_result(settings[1L, , drop = FALSE],
                           design$treated_sizes[1L], 1L, "", NA_real_, design)[0, ]
    summary <- results[, c("setting", "kappa_pi", "kappa_m", "n1", "n0",
                           "method"), drop = FALSE]
    for (name in c("replications", "successful", "failed")) {
      summary[[name]] <- integer()
    }
    for (name in c("coverage", "bias", "rmse")) summary[[name]] <- numeric()
  }
  execution$status <- state$status
  execution$completed_replications <- state$completed
  execution$planned_replications <- planned_replications
  execution$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
  if (!is.null(state$condition)) execution$condition <- state$condition
  metadata <- list(design = design, settings = settings,
                   quick_run = options$quick_run, targets = targets,
                   study_engine = study_reference,
                   saved_configurations = study_configs,
                   rng_order = "Sample size varies fastest within each replicate",
                   execution = execution, session = utils::sessionInfo())
  message(sprintf("Scenario 2 %s: %d/%d completed datasets retained in memory.",
                  state$status, state$completed, planned_replications))
  if (state$status == "failed") message("Execution error: ", state$condition)
  invisible(list(status = state$status, replications = results, summary = summary,
                 targets = targets, metadata = metadata))
}

scenario2_results <- run_scenario2(design, settings, options)
scenario2_summary <- scenario2_results$summary
scenario2_replications <- scenario2_results$replications
scenario2_targets <- scenario2_results$targets

# Keep completed runs in R, including partial results after an interruption.
scenario2_display <- if (nrow(scenario2_replications)) {
  build_simulation_report(scenario2_results, "Scenario 2")
} else {
  list(tables = list(), report = NULL)
}
scenario2_tables <- scenario2_display$tables
scenario2_report <- scenario2_display$report
if (!is.null(scenario2_report)) {
  show_simulation_report(scenario2_report)
  plot_scenario2_results(scenario2_summary, settings)
}
