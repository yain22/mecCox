# Run the first simulation using the bundled main-paper study engine.
# Run from any directory after installing mecCox:
#   Rscript path/to/mecCox/inst/reproduce/scenario1.R
# A short code-path check is available with --quick; it is not a paper result.
# Simulation runs use up to 20 workers by default; use --cores=1 for a serial run.
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
script_option <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
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

required_packages <- c("mecCox", "survival", "MASS", "rngtools")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Install these packages before starting: ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}

# Bundled configurations supply the main-paper study targets and fitting settings.
study_engine <- get_original_study_engine("scenario1")
study_reference <- original_study_metadata("scenario1")
study_keys <- as.character(c(2L, 3L, 4L))
study_configs <- stats::setNames(lapply(study_keys, function(key) {
  get_original_study_config("scenario1", key)
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
  control_multipliers = as.integer(study_keys),
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

if (options$quick_run) {
  design$treated_sizes <- study_config$sample_size_n1[1L]
  design$control_multipliers <- 2L
  message("Quick check: fewer datasets and sample-size cells; study learner ",
          "settings, folds, landmarks, and reference targets are retained.")
}

make_result <- function(multiplier, treated_count, replicate, method,
                        target, estimate = NA_real_, standard_error = NA_real_,
                        error = NA_character_) {
  data.frame(
    ratio = paste0("1:", multiplier),
    n1 = treated_count,
    n0 = multiplier * treated_count,
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

run_replication <- function(multiplier, treated_count, replicate, target, design) {
  key <- as.character(multiplier)
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
  methods[grepl("^MEC-Cox:", raw$Method)] <- "MEC-Cox"
  keep <- !is.na(methods)
  raw <- raw[keep, , drop = FALSE]
  methods <- methods[keep]
  if (!nrow(raw)) stop("The study engine returned no requested method rows.")
  answer <- make_result(multiplier, treated_count, raw$rep, methods,
                        raw$theta_true, raw$Estimate, raw$SE, raw$Error)
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

summarize_results <- function(results) {
  groups <- split(results,
                  interaction(results$ratio, results$n1, results$method,
                              drop = TRUE))
  summaries <- lapply(groups, function(group) {
    valid <- is.finite(group$estimate) &
      is.finite(group$standard_error)
    successful <- group[valid, , drop = FALSE]
    error <- successful$estimate - successful$target
    result <- group[1L, c("ratio", "n1", "n0", "method")]
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
  answer[order(answer$ratio, answer$n1,
               match(answer$method,
                     c("Naive", "Robust sandwich", "Corrected sandwich",
                       "MEC-Cox"))), ]
}

run_scenario1 <- function(design, options) {
  started <- proc.time()[["elapsed"]]
  worker_count <- choose_worker_count(options$cores, design$replications)
  planned_cells <- length(design$control_multipliers) * length(design$treated_sizes)
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
  message(sprintf("Scenario 1 uses %d worker(s); %d requested.",
                  worker_count, options$cores))
  message(sprintf(
    "Workload: %d ratios x %d sample sizes = %d cells; %d runs per cell = %d datasets in total.",
    length(design$control_multipliers), length(design$treated_sizes), planned_cells,
    design$replications, planned_replications
  ))
  targets <- data.frame(
    ratio = paste0("1:", design$control_multipliers),
    target = vapply(study_configs[as.character(design$control_multipliers)],
                     function(config) config$theta_true, numeric(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
  tryCatch({
    cluster <- start_simulation_cluster(worker_count, export_names,
                                        envir = environment(run_replication))
    cell_index <- 0L
    for (ratio_index in seq_along(design$control_multipliers)) {
      multiplier <- design$control_multipliers[ratio_index]
      target <- targets$target[ratio_index]
      message(sprintf("Ratio 1:%d: reference log-hazard ratio %.6f", multiplier, target))
      for (treated_count in design$treated_sizes) {
        cell_index <- cell_index + 1L
        cell_completed <- 0L
        last_progress <- proc.time()[["elapsed"]]
        message(sprintf("Running cell %d/%d: ratio 1:%d, n1=%d, n0=%d (%d runs)",
                        cell_index, planned_cells, multiplier,
                        treated_count, multiplier * treated_count,
                        design$replications))
        record_result <- function(replicate, value) {
          # Retain task order even when workers finish out of order.
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
          arguments = list(multiplier = multiplier, treated_count = treated_count,
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
    summary <- summarize_results(results)
    print(summary, row.names = FALSE, digits = 4)
    if (any(summary$failed > 0L)) {
      message("Some fits failed; inspect the error column in the individual simulation results.")
    }
  } else {
    results <- make_result(design$control_multipliers[1L],
                           design$treated_sizes[1L], 1L, "", NA_real_)[0, ]
    summary <- results[, c("ratio", "n1", "n0", "method"), drop = FALSE]
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
  metadata <- list(design = design, quick_run = options$quick_run, targets = targets,
                   study_engine = study_reference,
                   saved_configurations = study_configs,
                   rng_order = "Sample size varies fastest within each replicate",
                   execution = execution, session = utils::sessionInfo())
  message(sprintf("Scenario 1 %s: %d/%d completed datasets retained in memory.",
                  state$status, state$completed, planned_replications))
  if (state$status == "failed") message("Execution error: ", state$condition)
  invisible(list(status = state$status, replications = results, summary = summary,
                 targets = targets, metadata = metadata))
}

scenario1_results <- run_scenario1(design, options)
scenario1_summary <- scenario1_results$summary
scenario1_replications <- scenario1_results$replications
scenario1_targets <- scenario1_results$targets

# Keep completed runs in R, including partial results after an interruption.
scenario1_display <- if (nrow(scenario1_replications)) {
  build_simulation_report(scenario1_results, "Scenario 1")
} else {
  list(tables = list(), report = NULL)
}
scenario1_tables <- scenario1_display$tables
scenario1_report <- scenario1_display$report
if (!is.null(scenario1_report)) {
  show_simulation_report(scenario1_report)
  plot_scenario1_results(scenario1_summary)
}
