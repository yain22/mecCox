# Reproduce the paper's first simulation experiment with the public mecCox API.
# Run from any directory after installing mecCox:
#   Rscript path/to/mecCox/inst/reproduce/scenario1.R --output=simulation1-output
# A short code-path check is available with --quick; it is not a paper result.

suppressPackageStartupMessages({
  library(mecCox)
  library(survival)
})

arguments <- commandArgs(trailingOnly = TRUE)
quick_run <- "--quick" %in% arguments
output_option <- grep("^--output=", arguments, value = TRUE)
unknown <- arguments[!grepl("^(--quick|--output=.+)$", arguments)]
if (length(unknown) || length(output_option) > 1L) {
  stop("Use only --quick and one --output=directory argument.", call. = FALSE)
}

output_directory <- if (length(output_option)) {
  sub("^--output=", "", output_option)
} else {
  "scenario1-output"
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

design <- list(
  seed = 20260427L,
  replications = 1000L,
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
  design$replications <- 2L
  design$treated_sizes <- 200L
  design$control_multipliers <- 2L
  design$super_treated <- 2000L
  design$super_controls <- 4000L
  message("Quick check: reduced replications and target sample. Do not cite as a paper result.")
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

plot_results <- function(summary, output_file) {
  methods <- c("Naive", "Robust sandwich", "Corrected sandwich", "MEC-Cox")
  colors <- c("gray25", "#0072B2", "#009E73", "#D55E00")
  line_types <- c(3L, 2L, 4L, 1L)
  symbols <- c(4L, 1L, 2L, 16L)
  metrics <- c(coverage = "Coverage", bias = "Bias", rmse = "RMSE")
  ratios <- unique(summary$ratio)

  grDevices::pdf(output_file, width = 14, height = 3.1 * length(ratios) + 2.5)
  panel_count <- 3L * length(ratios)
  panel_layout <- matrix(seq_len(panel_count), ncol = 3L, byrow = TRUE)
  panel_layout <- rbind(panel_layout, rep(panel_count + 1L, 3L))
  graphics::layout(panel_layout,
                   heights = c(rep(1, length(ratios)), 0.20))
  old <- graphics::par(mar = c(4, 4.2, 3, 1),
                       oma = c(0, 0, 1.5, 0), las = 1)
  on.exit({
    graphics::par(old)
    grDevices::dev.off()
  }, add = TRUE)
  panel_index <- 0L

  for (ratio in ratios) {
    for (metric in names(metrics)) {
      panel_index <- panel_index + 1L
      panel_data <- summary[summary$ratio == ratio, , drop = FALSE]
      observed <- panel_data[[metric]][is.finite(panel_data[[metric]])]
      limits <- if (length(observed)) range(observed) else c(0, 1)
      if (diff(limits) < 1e-8) limits <- limits + c(-0.01, 0.01)
      if (metric == "coverage") limits <- range(c(limits, 0.95))
      if (metric == "bias") limits <- range(c(limits, 0))
      graphics::plot(range(panel_data$n1), limits, type = "n",
                     xlab = expression(n[1]), ylab = metrics[metric],
                     main = sprintf("(%s) %s, n1:n0 = %s",
                                    letters[panel_index], metrics[metric], ratio))
      if (metric == "coverage") graphics::abline(h = 0.95, col = "gray70")
      if (metric == "bias") graphics::abline(h = 0, col = "gray70")
      for (index in seq_along(methods)) {
        method_data <- panel_data[panel_data$method == methods[index], ]
        method_data <- method_data[order(method_data$n1), ]
        graphics::lines(method_data$n1, method_data[[metric]],
                        col = colors[index], lty = line_types[index], lwd = 1.5)
        graphics::points(method_data$n1, method_data[[metric]],
                         col = colors[index], pch = symbols[index])
      }
    }
  }
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend("center", legend = methods,
                   col = colors, lty = line_types, pch = symbols,
                   horiz = TRUE, bty = "n", cex = 0.85)
}

message("Computing the ATT Cox-projection reference target ...")
target <- compute_reference_target(design)
message(sprintf("Reference log-hazard ratio: %.6f", target))
saveRDS(list(design = design, quick_run = quick_run, target = target,
             session = utils::sessionInfo()),
        file.path(output_directory, "run_metadata.rds"))

result_cells <- list()
cell_index <- 0L
for (multiplier in design$control_multipliers) {
  for (treated_count in design$treated_sizes) {
    message(sprintf("Running n1=%d, n0=%d (%d replications)",
                    treated_count, multiplier * treated_count,
                    design$replications))
    cell_rows <- vector("list", design$replications)
    for (replicate in seq_len(design$replications)) {
      cell_rows[[replicate]] <- run_replication(
        multiplier, treated_count, replicate, target, design
      )
    }
    cell_results <- do.call(rbind, cell_rows)
    checkpoint <- sprintf("checkpoint_n1-%d_n0-%d.csv",
                          treated_count, multiplier * treated_count)
    utils::write.csv(cell_results,
                     file.path(output_directory, checkpoint), row.names = FALSE)
    message("Completed cell; checkpoint: ", checkpoint)
    cell_index <- cell_index + 1L
    result_cells[[cell_index]] <- cell_results
  }
}

results <- do.call(rbind, result_cells)
summary <- summarize_results(results)
utils::write.csv(results, file.path(output_directory, "replications.csv"),
                 row.names = FALSE)
utils::write.csv(summary, file.path(output_directory, "summary.csv"),
                 row.names = FALSE)
plot_results(summary, file.path(output_directory, "scenario1.pdf"))

print(summary, row.names = FALSE, digits = 4)
if (any(summary$failed > 0L)) {
  warning("Some fits failed; inspect the error column in replications.csv.",
          call. = FALSE)
}
message("Results written to: ", normalizePath(output_directory))
