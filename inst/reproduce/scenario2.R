# Reproduce the second simulation experiment with the public mecCox API.
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
if (!length(helper_candidates)) {
  stop("Cannot find simulation_helpers.R. Keep it beside scenario2.R and ",
       "use source('path/to/scenario2.R'), or install the current mecCox package.",
       call. = FALSE)
}
source(helper_candidates[1L], local = TRUE)
if (!exists("build_simulation_report", mode = "function", inherits = FALSE)) {
  stop("Update mecCox or keep the current simulation_helpers.R beside ",
       "scenario2.R before starting the simulation.", call. = FALSE)
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

required_packages <- c("dbarts", "ranger")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Scenario 2 requires the optional packages: ",
       paste(missing_packages, collapse = ", "),
       ". Install them before starting the simulation.", call. = FALSE)
}

design <- list(
  seed = 20260427L,
  replications = options$replications,
  treated_sizes = c(200L, 250L, 300L, 350L, 400L),
  control_multiplier = 4L,
  covariate_count = 10L,
  folds = 10L,
  landmarks = 5L,
  super_treated = 30000L,
  super_controls = 60000L,
  weibull_scale = 0.00008,
  weibull_shape = 2,
  censoring_rate = 0.0008,
  conditional_log_hr = log(0.70),
  bart_num_trees = 100L,
  bart_posterior_draws = 1000L,
  bart_burn_in = 500L,
  bart_shrinkage = 2,
  rsf_num_trees = 500L,
  rsf_min_node_size = 15L,
  rsf_auto_tune = TRUE
)

settings <- data.frame(
  setting_id = 1:3,
  setting = c("No nonlinearity", "Mild nonlinearity", "Severe nonlinearity"),
  kappa_pi = c(0, 1, 2),
  kappa_m = c(0, 2, 5),
  stringsAsFactors = FALSE
)

if (options$quick_run) {
  design$treated_sizes <- 200L
  design$super_treated <- 2000L
  design$super_controls <- 4000L
  design$bart_num_trees <- 25L
  design$bart_posterior_draws <- 50L
  design$bart_burn_in <- 25L
  design$rsf_num_trees <- 100L
  design$rsf_auto_tune <- FALSE
  message("Quick check: all three settings, reduced workload, target ",
          "sample, and tree workloads. Do not cite as a paper result.")
}

# The source mechanism includes the manuscript's nonlinear component. Cohort
# sizes are imposed by sampling from X | A, leaving the intercept at -0.2.
source_probability <- function(covariates, kappa_pi) {
  x1 <- covariates[, 1L]
  x2 <- covariates[, 2L]
  x3 <- covariates[, 3L]
  x4 <- covariates[, 4L]
  x5 <- covariates[, 5L]

  linear <- 0.75 * x1 + 0.75 * x2 + 0.65 * x3 +
    0.65 * x4 + 0.55 * x5
  nonlinear <- 0.70 * sin(1.25 * x1) + 0.45 * (x2^2 - 1) -
    0.55 * (as.numeric(x3 > 0) - 0.5) + 0.35 * x4 * x5 +
    0.25 * (cos(x1 + x2) - exp(-1))
  probability <- stats::plogis(-0.2 + linear + kappa_pi * nonlinear)
  pmin(pmax(probability, 0.02), 0.98)
}

draw_source_covariates <- function(treated_count, control_count, dimension,
                                   kappa_pi) {
  treated <- matrix(numeric(), nrow = 0L, ncol = dimension)
  controls <- matrix(numeric(), nrow = 0L, ncol = dimension)

  while (nrow(treated) < treated_count || nrow(controls) < control_count) {
    batch_count <- max(4000L, 4L * (treated_count + control_count))
    candidates <- matrix(stats::rnorm(batch_count * dimension),
                         nrow = batch_count, ncol = dimension)
    source <- stats::rbinom(batch_count, 1L,
                            source_probability(candidates, kappa_pi))
    if (nrow(treated) < treated_count) {
      treated <- rbind(treated, candidates[source == 1L, , drop = FALSE])
    }
    if (nrow(controls) < control_count) {
      controls <- rbind(controls, candidates[source == 0L, , drop = FALSE])
    }
  }

  treated <- treated[seq_len(treated_count), , drop = FALSE]
  controls <- controls[seq_len(control_count), , drop = FALSE]
  colnames(treated) <- colnames(controls) <- paste0("X", seq_len(dimension))
  list(treated = treated, controls = controls)
}

control_log_hazard <- function(covariates, kappa_m) {
  coefficients <- log(c(1.75, 1.75, 1.60, 1.60, 1.50, rep(1.25, 5L)))
  linear <- as.vector(covariates %*% coefficients)
  x1 <- covariates[, 1L]
  x2 <- covariates[, 2L]
  x3 <- covariates[, 3L]
  x4 <- covariates[, 4L]
  x5 <- covariates[, 5L]
  nonlinear <- 0.45 * sin(x2) + 0.35 * (x3^2 - 1) +
    0.30 * (as.numeric(x4 > 0) - 0.5) + 0.25 * x1 * x5 +
    0.20 * (cos(x2 + x5) - exp(-1))
  linear + kappa_m * nonlinear
}

draw_observed_data <- function(treated_count, control_count, setting, design) {
  source_covariates <- draw_source_covariates(
    treated_count, control_count, design$covariate_count, setting$kappa_pi
  )
  treated <- source_covariates$treated
  controls <- source_covariates$controls

  draw_event_time <- function(log_hazard) {
    uniform <- stats::runif(length(log_hazard))
    (-log(uniform) /
       (design$weibull_scale * exp(log_hazard)))^(1 / design$weibull_shape)
  }

  treated_event_time <- draw_event_time(
    control_log_hazard(treated, setting$kappa_m) + design$conditional_log_hr
  )
  control_event_time <- draw_event_time(
    control_log_hazard(controls, setting$kappa_m)
  )
  treated_censor_time <- stats::rexp(treated_count, design$censoring_rate)
  control_censor_time <- stats::rexp(control_count, design$censoring_rate)

  treated_data <- data.frame(
    source = 1L,
    time = pmin(treated_event_time, treated_censor_time),
    event = as.integer(treated_event_time <= treated_censor_time),
    treated, check.names = FALSE
  )
  control_data <- data.frame(
    source = 0L,
    time = pmin(control_event_time, control_censor_time),
    event = as.integer(control_event_time <= control_censor_time),
    controls, check.names = FALSE
  )
  rbind(treated_data, control_data)
}

# Each nonlinearity setting has its own marginal Cox projection. True source
# odds, rather than any fitted learner, define the superpopulation benchmark.
compute_reference_target <- function(setting, design) {
  set.seed(design$seed)
  data <- draw_observed_data(
    design$super_treated, design$super_controls, setting, design
  )
  covariates <- as.matrix(data[paste0("X", seq_len(design$covariate_count))])
  probability <- source_probability(covariates, setting$kappa_pi)
  odds <- probability / (1 - probability)
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
    error = error,
    stringsAsFactors = FALSE
  )
}

run_replication <- function(setting, treated_count, replicate, target, design) {
  # A replication owns its seed, independently of worker count and scheduling.
  replicate_seed <- design$seed + 100000L * setting$setting_id +
    100L * treated_count + replicate
  set.seed(replicate_seed)
  data <- draw_observed_data(
    treated_count, design$control_multiplier * treated_count, setting, design
  )
  covariates <- paste0("X", seq_len(design$covariate_count))

  # The conventional comparators retain logistic main-effect propensity scores.
  ipw <- tryCatch(
    fit_att_ipw_cox(data, "time", "event", "source", covariates,
                    ps_clip = c(0.01, 0.99)),
    error = function(condition) condition
  )
  ipw_labels <- c(naive = "Naive", robust = "Robust sandwich",
                  corrected = "Corrected sandwich")
  rows <- vector("list", length(ipw_labels) + 2L)
  for (index in seq_along(ipw_labels)) {
    if (inherits(ipw, "error")) {
      rows[[index]] <- make_result(
        setting, treated_count, replicate, ipw_labels[index], target, design,
        error = conditionMessage(ipw)
      )
    } else {
      variance_name <- names(ipw_labels)[index]
      rows[[index]] <- make_result(
        setting, treated_count, replicate, ipw_labels[index], target, design,
        estimate = ipw$theta, standard_error = ipw$se[variance_name]
      )
    }
  }

  # Both variants use identical folds, BART controls, and seeds. Their source
  # predictions therefore agree; only the control-survival learner changes.
  survival_learners <- c(cox = "MEC-Cox (BART/Cox)", rsf = "MEC-Cox (BART/RSF)")
  for (index in seq_along(survival_learners)) {
    learner <- names(survival_learners)[index]
    mec <- tryCatch(
      fit_mec_cox(
        data, "time", "event", "source", covariates,
        ps_learner = "bart", survival_learner = learner,
        n_folds = design$folds, n_landmarks = design$landmarks,
        seed = replicate_seed, ps_trim = c(0.01, 0.99),
        bart_num_trees = design$bart_num_trees,
        bart_posterior_draws = design$bart_posterior_draws,
        bart_burn_in = design$bart_burn_in,
        bart_shrinkage = design$bart_shrinkage,
        rsf_num_trees = design$rsf_num_trees,
        rsf_min_node_size = design$rsf_min_node_size,
        rsf_auto_tune = design$rsf_auto_tune
      ),
      error = function(condition) condition
    )
    row_index <- length(ipw_labels) + index
    if (inherits(mec, "error")) {
      rows[[row_index]] <- make_result(
        setting, treated_count, replicate, survival_learners[index],
        target, design, error = conditionMessage(mec)
      )
    } else {
      rows[[row_index]] <- make_result(
        setting, treated_count, replicate, survival_learners[index],
        target, design, estimate = mec$theta, standard_error = mec$se
      )
    }
  }
  do.call(rbind, rows)
}

summarize_scenario2_results <- function(results, settings) {
  groups <- split(results,
                  interaction(results$setting, results$n1, results$method,
                              drop = TRUE))
  summaries <- lapply(groups, function(group) {
    valid <- is.finite(group$estimate) &
      is.finite(group$standard_error) & group$standard_error > 0
    successful <- group[valid, , drop = FALSE]
    error <- successful$estimate - successful$target
    result <- group[1L, c("setting", "kappa_pi", "kappa_m", "n1", "n0",
                         "method")]
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
  methods <- c("Naive", "Robust sandwich", "Corrected sandwich",
               "MEC-Cox (BART/Cox)", "MEC-Cox (BART/RSF)")
  answer[order(match(answer$setting, settings$setting), answer$n1,
               match(answer$method, methods)), ]
}

run_scenario2 <- function(design, settings, options) {
  worker_count <- choose_worker_count(options$cores, design$replications)
  export_names <- c("source_probability", "draw_source_covariates",
                    "control_log_hazard", "draw_observed_data", "make_result",
                    "run_replication")
  cluster <- start_simulation_cluster(worker_count, export_names,
                                      envir = environment(run_replication))
  if (!is.null(cluster)) {
    on.exit(parallel::stopCluster(cluster), add = TRUE)
  }
  execution <- list(
    requested_cores = options$cores,
    workers = worker_count,
    backend = if (is.null(cluster)) "serial" else "PSOCK",
    rng_kind = RNGkind()
  )
  message(sprintf("Scenario 2 uses %d worker(s); %d requested.",
                  worker_count, options$cores))

  message("Computing the three ATT Cox-projection reference targets ...")
  targets <- settings
  targets$target <- vapply(seq_len(nrow(settings)), function(index) {
    target <- compute_reference_target(settings[index, , drop = FALSE], design)
    message(sprintf("%s: reference log-hazard ratio %.6f",
                    settings$setting[index], target))
    target
  }, numeric(1))
  metadata <- list(design = design, settings = settings,
                   quick_run = options$quick_run, targets = targets,
                   execution = execution, session = utils::sessionInfo())

  result_cells <- list()
  cell_index <- 0L
  for (setting_index in seq_len(nrow(settings))) {
    setting <- settings[setting_index, , drop = FALSE]
    target <- targets$target[setting_index]
    for (treated_count in design$treated_sizes) {
      control_count <- design$control_multiplier * treated_count
      message(sprintf("Running %s: n1=%d, n0=%d (%d replications)",
                      setting$setting, treated_count, control_count,
                      design$replications))
      cell_rows <- run_simulation_replications(
        design$replications, run_replication,
        arguments = list(setting = setting, treated_count = treated_count,
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
  summary <- summarize_scenario2_results(results, settings)
  print(summary, row.names = FALSE, digits = 4)
  if (any(summary$failed > 0L)) {
    warning("Some fits failed; inspect the error column in the replication results.",
            call. = FALSE)
  }
  invisible(list(replications = results, summary = summary,
                 targets = targets, metadata = metadata))
}

scenario2_results <- run_scenario2(design, settings, options)
scenario2_summary <- scenario2_results$summary
scenario2_replications <- scenario2_results$replications
scenario2_targets <- scenario2_results$targets

# Keep every replication in R; display the complete aggregated results.
scenario2_display <- build_simulation_report(scenario2_results, "Scenario 2")
scenario2_tables <- scenario2_display$tables
scenario2_report <- scenario2_display$report
show_simulation_report(scenario2_report)
plot_scenario2_results(scenario2_summary, settings)
