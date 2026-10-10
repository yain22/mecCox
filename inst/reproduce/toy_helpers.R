# Shared, standalone engine for the two manuscript toy examples.
# The methods below intentionally retain the manuscript's specialized binary-X1
# source model, normalized KL calibration, and independent weighted Cox fits.

parse_toy_arguments <- function(arguments = character(), quick_run = FALSE,
                                replications = 10000L, seed = 20261007L,
                                output_dir = NULL) {
  allowed <- arguments == "--quick" |
    grepl("^--replications=", arguments) | grepl("^--output-dir=", arguments)
  if (any(!allowed)) {
    stop("Unknown argument(s): ", paste(arguments[!allowed], collapse = ", "),
         call. = FALSE)
  }
  for (prefix in c("--replications=", "--output-dir=")) {
    if (sum(startsWith(arguments, prefix)) > 1L) {
      stop("Supply ", prefix, " only once.", call. = FALSE)
    }
  }
  if (!is.logical(quick_run) || length(quick_run) != 1L || is.na(quick_run)) {
    stop("quick_run must be TRUE or FALSE.", call. = FALSE)
  }
  quick_run <- quick_run || "--quick" %in% arguments
  count_argument <- arguments[startsWith(arguments, "--replications=")]
  if (length(count_argument)) {
    replications <- suppressWarnings(as.numeric(sub("^--replications=", "", count_argument)))
  } else if (quick_run && identical(as.numeric(replications), 10000)) {
    replications <- 20L
  }
  if (!is.numeric(replications) || length(replications) != 1L ||
      !is.finite(replications) || replications < 2 ||
      replications != floor(replications) || replications > .Machine$integer.max) {
    stop("replications must be an integer of at least 2.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) {
    stop("seed must be a nonnegative integer supported by set.seed().", call. = FALSE)
  }
  directory_argument <- arguments[startsWith(arguments, "--output-dir=")]
  if (length(directory_argument)) output_dir <- sub("^--output-dir=", "", directory_argument)
  if (!is.null(output_dir) && (!is.character(output_dir) ||
      length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir))) {
    stop("output_dir must be NULL or a nonempty directory path.", call. = FALSE)
  }
  list(quick_run = quick_run, replications = as.integer(replications),
       seed = as.integer(seed), output_dir = output_dir)
}

check_toy_packages <- function() {
  packages <- c("survival", "ggplot2", "patchwork")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install these packages before running the toy examples: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

toy_one_replication <- function(design) {
  n_treated <- design$n_treated
  n_controls <- design$n_controls
  t_star <- design$t_star
  baseline_hazard <- design$baseline_hazard
  censor_hazard <- design$censor_hazard
  log_hr_x1 <- design$log_hr_x1
  log_hr_x2 <- design$log_hr_x2

  # Keep the paper's exact draw order, including X2 when its coefficient is zero.
  x1_t <- rbinom(n_treated, 1, 0.7)
  x1_c <- rbinom(n_controls, 1, 0.3)
  x2_t <- rnorm(n_treated)
  x2_c <- rnorm(n_controls)
  n_t_by_x1 <- tabulate(x1_t + 1L, nbins = 2L)
  n_c_by_x1 <- tabulate(x1_c + 1L, nbins = 2L)
  if (any(n_c_by_x1 == 0L) || any(n_t_by_x1 == 0L)) stop("Empty source stratum")

  # These stratum ratios are the exact fitted odds from logistic source ~ X1.
  d_c <- n_t_by_x1[x1_c + 1L] / n_c_by_x1[x1_c + 1L]
  stopifnot(abs(sum(d_c) - n_treated) < 1e-10,
            abs(sum(d_c * x1_c) / n_treated - mean(x1_t)) < 1e-10)
  psi <- function(x1, x2) {
    exp(-baseline_hazard * t_star * exp(log_hr_x1 * x1 + log_hr_x2 * x2))
  }
  psi_t <- psi(x1_t, x2_t)
  psi_c <- psi(x1_c, x2_c)
  basis_c <- cbind(x1_c, psi_c)
  target <- c(mean(x1_t), mean(psi_t))
  baseline_residual <- drop(crossprod(basis_c, d_c)) / n_treated - target
  if (log_hr_x2 == 0) stopifnot(max(abs(baseline_residual)) < 1e-12)

  # The intercept constraint is enforced by normalization; preserve X1 balance
  # and add exact balance of the oracle control-survival score at t_star.
  calibrated_weights <- function(lambda) {
    eta <- drop(basis_c %*% lambda)
    u <- d_c * exp(eta - max(eta))
    n_treated * u / sum(u)
  }
  balance <- function(lambda) {
    drop(crossprod(basis_c, calibrated_weights(lambda))) / n_treated - target
  }
  lambda <- c(0, 0)
  calibration_iterations <- 0L
  for (iteration in seq_len(40L)) {
    w <- calibrated_weights(lambda)
    residual <- drop(crossprod(basis_c, w)) / n_treated - target
    # Test feasibility before Newton: the X1-only oracle basis is redundant.
    if (max(abs(residual)) < 1e-12) break
    centered <- sweep(basis_c, 2, target + residual)
    jacobian <- crossprod(centered, centered * w) / n_treated
    step <- drop(solve(jacobian, residual))
    step_scale <- 1
    repeat {
      candidate <- lambda - step_scale * step
      if (max(abs(balance(candidate))) < max(abs(residual))) break
      step_scale <- step_scale / 2
      if (step_scale < 1e-8) stop("Calibration line search failed")
    }
    lambda <- candidate
    calibration_iterations <- calibration_iterations + 1L
  }
  if (max(abs(balance(lambda))) > 1e-9) stop("Calibration did not converge")
  w_c <- calibrated_weights(lambda)
  stopifnot(all(is.finite(d_c)), all(is.finite(w_c)), all(d_c > 0), all(w_c > 0))

  x1 <- c(x1_t, x1_c)
  x2 <- c(x2_t, x2_c)
  a <- c(rep(1L, n_treated), rep(0L, n_controls))
  # Both cohorts share the same event hazard: the target log hazard ratio is zero.
  event_time <- rexp(length(a), rate = baseline_hazard * exp(log_hr_x1 * x1 + log_hr_x2 * x2))
  censor_time <- rexp(length(a), rate = censor_hazard)
  followup <- pmin(event_time, censor_time, t_star)
  event <- as.integer(event_time <= censor_time & event_time <= t_star)
  fit_cox <- function(weights) {
    warning_count <- 0L
    fit <- withCallingHandlers(
      survival::coxph(survival::Surv(followup, event) ~ a, weights = weights,
                      ties = "breslow", robust = FALSE),
      warning = function(warning) {
        warning_count <<- warning_count + 1L
        invokeRestart("muffleWarning")
      }
    )
    c(theta = unname(stats::coef(fit)[1]), warnings = warning_count)
  }
  # Both fits are computed independently, even when the calibrated weights agree.
  fit_ipw <- fit_cox(c(rep(1, n_treated), d_c))
  fit_mec <- fit_cox(c(rep(1, n_treated), w_c))
  result <- c(
    theta_ipw = unname(fit_ipw["theta"]), theta_mec = unname(fit_mec["theta"]),
    imbalance_ipw = target[2] - sum(d_c * psi_c) / n_treated,
    imbalance_mec = target[2] - sum(w_c * psi_c) / n_treated,
    x1_imbalance_ipw = target[1] - sum(d_c * x1_c) / n_treated,
    x1_imbalance_mec = target[1] - sum(w_c * x1_c) / n_treated,
    x2_imbalance_ipw = mean(x2_t) - sum(d_c * x2_c) / n_treated,
    x2_imbalance_mec = mean(x2_t) - sum(w_c * x2_c) / n_treated,
    ess_ipw = sum(d_c)^2 / sum(d_c^2), ess_mec = sum(w_c)^2 / sum(w_c^2),
    cv_ipw = sd(d_c) / mean(d_c), cv_mec = sd(w_c) / mean(w_c),
    min_weight_ipw = min(d_c), min_weight_mec = min(w_c),
    max_weight_ipw = max(d_c), max_weight_mec = max(w_c),
    max_abs_weight_difference = max(abs(w_c - d_c)),
    theta_difference = unname(fit_mec["theta"] - fit_ipw["theta"]),
    max_abs_baseline_balance = max(abs(baseline_residual)),
    max_abs_calibrated_balance = max(abs(balance(lambda))),
    normalization_gap_ipw = sum(d_c) - n_treated,
    normalization_gap_mec = sum(w_c) - n_treated,
    calibration_iterations = calibration_iterations, calibration_checks = iteration,
    max_abs_lambda = max(abs(lambda)),
    cox_warnings_ipw = unname(fit_ipw["warnings"]),
    cox_warnings_mec = unname(fit_mec["warnings"])
  )
  stopifnot(all(is.finite(result)))
  result
}

build_toy_figure <- function(replications, scenario) {
  n_rep <- nrow(replications)
  no_gain <- identical(scenario, "toy_no_precision_gain")
  ipw_blue <- "#0072B2"
  mec_orange <- "#D55E00"
  method_levels <- c("ATT-IPW Cox", "MEC-Cox")
  method_colors <- c("ATT-IPW Cox" = ipw_blue, "MEC-Cox" = mec_orange)
  method_lines <- c("ATT-IPW Cox" = "solid", "MEC-Cox" = "dashed")
  theme_pub <- ggplot2::theme_bw(base_size = 14, base_family = "sans") +
    ggplot2::theme(
      aspect.ratio = 1, legend.position = "bottom", legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      panel.grid.major = ggplot2::element_line(color = "#e8e8e8", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0)
    )
  if (no_gain) theme_pub <- theme_pub + ggplot2::theme(legend.key.width = grid::unit(1.5, "cm"))
  gap <- replications[, "imbalance_ipw"]
  correction <- replications[, "theta_mec"] - replications[, "theta_ipw"]
  if (no_gain) {
    # Display floating-point residuals below 1e-10 as zero only in (a) and (b).
    # Keep the raw values in the result and the estimates used for panel (c).
    stopifnot(max(abs(gap)) < 1e-10,
              max(abs(replications[, "imbalance_mec"])) < 1e-10,
              max(abs(correction)) < 1e-10,
              max(replications[, "max_abs_weight_difference"]) < 1e-10)
    gap[] <- 0
    correction[] <- 0
    mass_data <- data.frame(method = factor(method_levels, levels = method_levels),
                            gap = c(0, 0), probability = c(1, 1))
    panel_a <- ggplot2::ggplot(mass_data) +
      ggplot2::geom_segment(ggplot2::aes(x = gap, xend = gap, y = 0, yend = probability,
                                        color = method, linetype = method),
                            linewidth = 1.1, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = method_colors) +
      ggplot2::scale_linetype_manual(values = method_lines) +
      ggplot2::scale_x_continuous(limits = c(-0.10, 0.10),
                                 breaks = c(-0.10, -0.05, 0, 0.05, 0.10)) +
      ggplot2::scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.25, 0.5, 0.75, 1),
                                 expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::labs(y = "Probability mass")
  } else {
    panel_a <- ggplot2::ggplot(data.frame(imbalance = gap), ggplot2::aes(x = imbalance)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 32,
                              fill = ipw_blue, alpha = 0.53, color = "white", linewidth = 0.2) +
      ggplot2::geom_vline(xintercept = 0, color = mec_orange, linewidth = 1.1) +
      ggplot2::labs(y = "Density")
  }
  panel_a <- panel_a + ggplot2::labs(title = "(a) Prognostic imbalance",
    x = expression("Empirical prognostic-score gap," ~ hat(Delta)[Psi](v))) +
    theme_pub + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 10))

  correction_data <- data.frame(score_imbalance = gap, estimator_change = correction)
  panel_b <- ggplot2::ggplot(correction_data,
                             ggplot2::aes(x = score_imbalance, y = estimator_change)) +
    ggplot2::geom_hline(yintercept = 0, color = "gray55", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, color = "gray55", linewidth = 0.4)
  if (no_gain) {
    panel_b <- panel_b + ggplot2::geom_point(color = ipw_blue, size = 1.5) +
      ggplot2::annotate("label", x = 0, y = 0.24,
        label = paste(format(n_rep, big.mark = ",", scientific = FALSE), "coincident points"),
        size = 3.6, family = "sans", fill = "white", linewidth = 0,
        label.padding = grid::unit(0.06, "lines"), label.r = grid::unit(0, "lines")) +
      ggplot2::scale_x_continuous(limits = c(-0.10, 0.10),
                                 breaks = c(-0.10, -0.05, 0, 0.05, 0.10)) +
      ggplot2::scale_y_continuous(limits = c(-0.35, 0.35),
                                 breaks = c(-0.3, -0.15, 0, 0.15, 0.3))
  } else {
    panel_b <- panel_b + ggplot2::geom_point(color = ipw_blue, alpha = 0.28, size = 0.75) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                           color = mec_orange, linewidth = 0.9)
  }
  panel_b <- panel_b + ggplot2::labs(title = "(b) Cox estimate correction",
    x = expression("Initial prognostic-score gap," ~ hat(Delta)[Psi](hat(d))),
    y = expression(hat(theta)[MEC] - hat(theta)[IPW])) +
    theme_pub + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 10))

  estimate_data <- data.frame(
    estimate = c(replications[, "theta_ipw"], replications[, "theta_mec"]),
    method = factor(rep(method_levels, each = n_rep), levels = method_levels)
  )
  panel_c <- ggplot2::ggplot(estimate_data, ggplot2::aes(x = estimate, color = method))
  if (no_gain) {
    density_max <- max(stats::density(replications[, "theta_ipw"])$y,
                       stats::density(replications[, "theta_mec"])$y)
    panel_c <- panel_c + ggplot2::aes(linetype = method) +
      ggplot2::geom_vline(xintercept = 0, color = "gray35", linetype = "dashed", linewidth = 0.55) +
      ggplot2::geom_density(linewidth = 1.05, key_glyph = "path") +
      ggplot2::scale_linetype_manual(values = method_lines) +
      ggplot2::coord_cartesian(xlim = c(-0.6, 0.6), ylim = c(0, density_max * 1.22))
  } else {
    panel_c <- panel_c + ggplot2::geom_density(linewidth = 1.05, key_glyph = "path") +
      ggplot2::geom_vline(xintercept = 0, color = "gray35", linetype = "dashed", linewidth = 0.55)
  }
  panel_c <- panel_c + ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::labs(title = "(c) ATT log hazard ratio", x = expression(hat(theta)[ATT]), y = "Density") +
    theme_pub
  panel_a + panel_b + panel_c + patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

run_toy_example <- function(scenario = c("toy_precision_gain", "toy_no_precision_gain"),
                            replications = 10000L, seed = 20261007L,
                            output_dir = NULL, display = interactive() || !is.null(grDevices::dev.list())) {
  scenario <- match.arg(scenario)
  config <- parse_toy_arguments(replications = replications, seed = seed, output_dir = output_dir)
  check_toy_packages()
  design <- list(n_treated = 200L, n_controls = 400L, t_star = 5,
                 baseline_hazard = 0.08, censor_hazard = 0.03, log_hr_x1 = 0.5,
                 log_hr_x2 = if (scenario == "toy_precision_gain") 1 else 0)
  # Explicit defaults reproduce the serial manuscript simulation stream.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(config$seed)
  started <- proc.time()["elapsed"]
  rows <- vector("list", config$replications)
  for (r in seq_len(config$replications)) {
    rows[[r]] <- toy_one_replication(design)
    if (r %% 1000L == 0L) {
      message(sprintf("%s: completed %d/%d replications (%.1f seconds).",
                      scenario, r, config$replications, proc.time()["elapsed"] - started))
    }
  }
  values <- do.call(rbind, rows)
  stopifnot(all(is.finite(values)), max(abs(values[, "imbalance_mec"])) < 1e-8,
            max(abs(values[, "x1_imbalance_mec"])) < 1e-8)
  summary <- data.frame(
    estimator = c("ATT-IPW Cox", "MEC-Cox (oracle score)"),
    mean_log_hr = colMeans(values[, c("theta_ipw", "theta_mec")]),
    empirical_sd = apply(values[, c("theta_ipw", "theta_mec")], 2, sd),
    rmse = sqrt(colMeans(values[, c("theta_ipw", "theta_mec")]^2)),
    mean_control_ess = colMeans(values[, c("ess_ipw", "ess_mec")]),
    mean_control_cv = colMeans(values[, c("cv_ipw", "cv_mec")]),
    mean_min_control_weight = colMeans(values[, c("min_weight_ipw", "min_weight_mec")]),
    mean_max_control_weight = colMeans(values[, c("max_weight_ipw", "max_weight_mec")]),
    row.names = NULL
  )
  diagnostics <- data.frame(
    scenario = scenario, seed = config$seed, replications = config$replications,
    variance_ratio_mec_ipw = var(values[, "theta_mec"]) / var(values[, "theta_ipw"]),
    max_abs_weight_difference = max(values[, "max_abs_weight_difference"]),
    max_abs_theta_difference = max(abs(values[, "theta_difference"])),
    max_abs_prognostic_gap_ipw = max(abs(values[, "imbalance_ipw"])),
    max_abs_prognostic_gap_mec = max(abs(values[, "imbalance_mec"])),
    max_abs_x1_gap_ipw = max(abs(values[, "x1_imbalance_ipw"])),
    max_abs_x1_gap_mec = max(abs(values[, "x1_imbalance_mec"])),
    max_abs_baseline_balance = max(values[, "max_abs_baseline_balance"]),
    max_abs_calibrated_balance = max(values[, "max_abs_calibrated_balance"]),
    max_abs_normalization_gap_ipw = max(abs(values[, "normalization_gap_ipw"])),
    max_abs_normalization_gap_mec = max(abs(values[, "normalization_gap_mec"])),
    mean_calibration_iterations = mean(values[, "calibration_iterations"]),
    max_calibration_iterations = max(values[, "calibration_iterations"]),
    max_calibration_checks = max(values[, "calibration_checks"]),
    max_abs_lambda = max(values[, "max_abs_lambda"]),
    total_cox_warnings_ipw = sum(values[, "cox_warnings_ipw"]),
    total_cox_warnings_mec = sum(values[, "cox_warnings_mec"])
  )
  figure <- build_toy_figure(values, scenario)
  metadata <- list(scenario = scenario, seed = config$seed,
                   replications = config$replications, design = design,
                   target_log_hr = 0, rng_kind = RNGkind(), session_info = utils::sessionInfo(),
                   figure_width = 10.6, figure_height = 4.75, png_dpi = 300,
                   plot_zero_tolerance = if (design$log_hr_x2 == 0) 1e-10 else NULL)
  result <- list(replications = values, summary = summary, diagnostics = diagnostics,
                  figure = figure, metadata = metadata)
  print(summary, row.names = FALSE, digits = 6)
  cat(sprintf("Variance ratio MEC/IPW: %.6f\n", diagnostics$variance_ratio_mec_ipw))
  cat(sprintf("Maximum absolute calibrated prognostic gap: %.3e\n",
              diagnostics$max_abs_prognostic_gap_mec))
  cat(sprintf("Cox-fit warnings (ATT-IPW / MEC): %d / %d\n",
              diagnostics$total_cox_warnings_ipw, diagnostics$total_cox_warnings_mec))
  # Creating the object does not open a device or leave Rplots.pdf in Rscript runs.
  if (isTRUE(display)) print(figure)
  if (!is.null(config$output_dir)) export_toy_results(result, config$output_dir)
  invisible(result)
}

export_toy_results <- function(result, output_dir) {
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("Cannot create output directory: ", output_dir, call. = FALSE)
  }
  prefix <- file.path(output_dir, result$metadata$scenario)
  for (name in c("replications", "summary", "diagnostics")) {
    suffix <- if (name == "replications") "replicates" else name
    utils::write.csv(result[[name]], paste0(prefix, "_", suffix, ".csv"), row.names = FALSE)
  }
  writeLines(capture.output(print(result$metadata$session_info)), paste0(prefix, "_sessionInfo.txt"))
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(paste0(prefix, ".pdf"), plot = result$figure, device = pdf_device,
                  width = 10.6, height = 4.75, units = "in", bg = "white", family = "sans")
  ggplot2::ggsave(paste0(prefix, ".png"), plot = result$figure,
                  width = 10.6, height = 4.75, units = "in", dpi = 300, bg = "white")
  message("Saved toy outputs in ", normalizePath(output_dir, winslash = "/"))
  invisible(result)
}
