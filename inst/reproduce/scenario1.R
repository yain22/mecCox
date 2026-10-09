# Run the paper's first simulation experiment with the public mecCox API.
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
if (!length(helper_candidates)) {
  stop("Cannot find simulation_helpers.R. Keep it beside scenario1.R and ",
       "use source('path/to/scenario1.R'), or install the current mecCox package.",
       call. = FALSE)
}
source(helper_candidates[1L], local = TRUE)
if (!exists("build_simulation_report", mode = "function", inherits = FALSE)) {
  stop("Update mecCox or keep the current simulation_helpers.R beside ",
       "scenario1.R before starting the simulation.", call. = FALSE)
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

suppressPackageStartupMessages({
  library(mecCox)
  library(survival)
})

quick_run <- options$quick_run

design <- list(
  seed = 20260427L,
  replications = options$replications,
  treated_sizes = c(200L, 250L, 300L, 350L, 400L),
  control_multipliers = c(2L, 3L, 4L),
  covariate_count = 50L,
  folds = 10L,
  landmarks = 5L,
  super_treated = 30000L,
  super_controls = 60000L,
  weibull_scale = 0.00008,
  weibull_shape = 2,
  censoring_rate = 0.0008,
  conditional_log_hr = log(0.70)
)

if (quick_run) {
  design$treated_sizes <- 200L
  design$control_multipliers <- 2L
  design$super_treated <- 2000L
  design$super_controls <- 4000L
  message("Quick check: reduced workload and target sample. Do not cite as a paper result.")
}

# The first five covariates drive source membership. The true probability is
# clipped for overlap, exactly as in the manuscript's simulation setup.
source_probability <- function(covariates) {
  linear_predictor <- -0.2 +
    0.75 * covariates[, 1L] + 0.75 * covariates[, 2L] +
    0.65 * covariates[, 3L] + 0.65 * covariates[, 4L] +
    0.55 * covariates[, 5L]
  probability <- stats::plogis(linear_predictor)
  pmin(pmax(probability, 0.02), 0.98)
}

# Sample from X | A until the requested two cohort sizes have been reached.
# Keeping the -0.2 source-model intercept fixed is important: the ratios are
# imposed by sampling, not by retuning the source mechanism.
draw_source_covariates <- function(treated_count, control_count, dimension) {
  treated <- matrix(numeric(), nrow = 0L, ncol = dimension)
  controls <- matrix(numeric(), nrow = 0L, ncol = dimension)

  while (nrow(treated) < treated_count || nrow(controls) < control_count) {
    batch_count <- max(4000L, 4L * (treated_count + control_count))
    candidates <- matrix(stats::rnorm(batch_count * dimension),
                         nrow = batch_count, ncol = dimension)
    source <- stats::rbinom(batch_count, 1L,
                            source_probability(candidates))
    treated <- rbind(treated, candidates[source == 1L, , drop = FALSE])
    controls <- rbind(controls, candidates[source == 0L, , drop = FALSE])
  }

  treated <- treated[seq_len(treated_count), , drop = FALSE]
  controls <- controls[seq_len(control_count), , drop = FALSE]
  colnames(treated) <- colnames(controls) <- paste0("X", seq_len(dimension))
  list(treated = treated, controls = controls)
}

# The first ten covariates have linear prognostic effects; X11--X50 are noise.
control_log_hazard <- function(covariates) {
  coefficients <- numeric(ncol(covariates))
  coefficients[1:10] <- log(c(1.75, 1.75, 1.60, 1.60, 1.50,
                            rep(1.25, 5L)))
  as.vector(covariates %*% coefficients)
}

draw_observed_data <- function(treated_count, control_count, design) {
  source_covariates <- draw_source_covariates(
    treated_count, control_count, design$covariate_count
  )
  treated <- source_covariates$treated
  controls <- source_covariates$controls

  draw_event_time <- function(log_hazard) {
    uniform <- stats::runif(length(log_hazard))
    (-log(uniform) /
       (design$weibull_scale * exp(log_hazard)))^(1 / design$weibull_shape)
  }

  treated_event_time <- draw_event_time(
    control_log_hazard(treated) + design$conditional_log_hr
  )
  control_event_time <- draw_event_time(control_log_hazard(controls))
  treated_censor_time <- stats::rexp(treated_count, design$censoring_rate)
  control_censor_time <- stats::rexp(control_count, design$censoring_rate)

  treated_data <- data.frame(
    source = 1L,
    time = pmin(treated_event_time, treated_censor_time),
    event = as.integer(treated_event_time <= treated_censor_time),
    treated,
    check.names = FALSE
  )
  control_data <- data.frame(
    source = 0L,
    time = pmin(control_event_time, control_censor_time),
    event = as.integer(control_event_time <= control_censor_time),
    controls,
    check.names = FALSE
  )
  rbind(treated_data, control_data)
}

# The conditional log-HR is not the benchmark after covariate marginalization.
# Approximate the ATT-weighted marginal Cox projection with a superpopulation
# and the *true* source odds, using the same censoring law as the simulation.
compute_reference_target <- function(design) {
  set.seed(design$seed)
  data <- draw_observed_data(
    design$super_treated, design$super_controls, design
  )
  covariates <- as.matrix(data[paste0("X", seq_len(design$covariate_count))])
  odds <- source_probability(covariates) /
    (1 - source_probability(covariates))
  control_rows <- data$source == 0L
  weights <- rep(1, nrow(data))
  weights[control_rows] <- design$super_treated * odds[control_rows] /
    sum(odds[control_rows])

  fit <- survival::coxph(
    survival::Surv(time, event) ~ source,
    data = data, weights = weights, ties = "breslow", robust = FALSE,
    control = survival::coxph.control(timefix = FALSE)
  )
  unname(stats::coef(fit)["source"])
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
    error = error,
    stringsAsFactors = FALSE
  )
}

run_replication <- function(multiplier, treated_count, replicate, target,
                            design) {
  # Separate, deterministic seeds keep a rerun stable even if grid order changes.
  replicate_seed <- design$seed + 100000L * multiplier +
    100L * treated_count + replicate
  set.seed(replicate_seed)
  data <- draw_observed_data(treated_count,
                             multiplier * treated_count, design)
  covariates <- paste0("X", seq_len(design$covariate_count))

  ipw <- tryCatch(
    fit_att_ipw_cox(data, "time", "event", "source", covariates,
                    ps_clip = c(0.01, 0.99)),
    error = function(condition) condition
  )
  ipw_labels <- c(naive = "Naive", robust = "Robust sandwich",
                  corrected = "Corrected sandwich")
  rows <- vector("list", length(ipw_labels) + 1L)
  if (inherits(ipw, "error")) {
    for (index in seq_along(ipw_labels)) {
      rows[[index]] <- make_result(
        multiplier, treated_count, replicate, ipw_labels[index],
        target, error = conditionMessage(ipw)
      )
    }
  } else {
    for (index in seq_along(ipw_labels)) {
      variance_name <- names(ipw_labels)[index]
      rows[[index]] <- make_result(
        multiplier, treated_count, replicate, ipw_labels[index], target,
        estimate = ipw$theta, standard_error = ipw$se[variance_name]
      )
    }
  }

  mec <- tryCatch(
    fit_mec_cox(
      data, "time", "event", "source", covariates,
      ps_learner = "glm", survival_learner = "cox",
      n_folds = design$folds, n_landmarks = design$landmarks,
      seed = replicate_seed, ps_trim = c(0.01, 0.99)
    ),
    error = function(condition) condition
  )
  mec_index <- length(rows)
  if (inherits(mec, "error")) {
    rows[[mec_index]] <- make_result(
      multiplier, treated_count, replicate, "MEC-Cox", target,
      error = conditionMessage(mec)
    )
  } else {
    rows[[mec_index]] <- make_result(
      multiplier, treated_count, replicate, "MEC-Cox", target,
      estimate = mec$theta, standard_error = mec$se
    )
  }
  do.call(rbind, rows)
}

summarize_results <- function(results) {
  groups <- split(results,
                  interaction(results$ratio, results$n1, results$method,
                              drop = TRUE))
  summaries <- lapply(groups, function(group) {
    valid <- is.finite(group$estimate) &
      is.finite(group$standard_error) & group$standard_error > 0
    successful <- group[valid, , drop = FALSE]
    error <- successful$estimate - successful$target
    result <- group[1L, c("ratio", "n1", "n0", "method")]
    result$replications <- nrow(group)
    result$successful <- nrow(successful)
    result$failed <- nrow(group) - nrow(successful)
    result$coverage <- if (length(error)) mean(
      abs(error) <= stats::qnorm(0.975) * successful$standard_error
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
  worker_count <- choose_worker_count(options$cores, design$replications)
  message(sprintf("Using %d worker(s); requested %d.",
                  worker_count, options$cores))
  message("Computing the ATT Cox-projection reference target ...")
  target <- compute_reference_target(design)
  message(sprintf("Reference log-hazard ratio: %.6f", target))
  execution <- list(requested_cores = options$cores, workers = worker_count,
                    backend = if (worker_count > 1L) "PSOCK" else "serial",
                    rng_kind = RNGkind())
  metadata <- list(design = design, quick_run = options$quick_run, target = target,
                   execution = execution, session = utils::sessionInfo())

  worker_functions <- c("source_probability", "draw_source_covariates",
                        "control_log_hazard", "draw_observed_data",
                        "make_result", "run_replication")
  cluster <- start_simulation_cluster(
    worker_count, worker_functions, envir = environment(run_replication)
  )
  on.exit({
    if (!is.null(cluster)) parallel::stopCluster(cluster)
  }, add = TRUE)

  result_cells <- list()
  cell_index <- 0L
  for (multiplier in design$control_multipliers) {
    for (treated_count in design$treated_sizes) {
      message(sprintf("Running n1=%d, n0=%d (%d runs)",
                      treated_count, multiplier * treated_count,
                      design$replications))
      cell_rows <- run_simulation_replications(
        design$replications, run_replication,
        arguments = list(multiplier = multiplier, treated_count = treated_count,
                         target = target, design = design),
        cluster = cluster
      )
      cell_results <- do.call(rbind, cell_rows)
      message("Completed cell.")
      cell_index <- cell_index + 1L
      result_cells[[cell_index]] <- cell_results
    }
  }

  results <- do.call(rbind, result_cells)
  summary <- summarize_results(results)

  print(summary, row.names = FALSE, digits = 4)
  if (any(summary$failed > 0L)) {
    warning("Some fits failed; inspect the error column in the individual simulation results.",
            call. = FALSE)
  }
  invisible(list(replications = results, summary = summary, metadata = metadata))
}

scenario1_results <- run_scenario1(design, options)
scenario1_summary <- scenario1_results$summary
scenario1_replications <- scenario1_results$replications

# Keep every run in R; display the complete aggregated results.
scenario1_display <- build_simulation_report(scenario1_results, "Scenario 1")
scenario1_tables <- scenario1_display$tables
scenario1_report <- scenario1_display$report
show_simulation_report(scenario1_report)
plot_scenario1_results(scenario1_summary)
