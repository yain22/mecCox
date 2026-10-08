# Execution helpers shared by the two simulation reproduction scripts.
# The scientific generators and fitting calls remain in each scenario script.

parse_simulation_arguments <- function(default_output,
                                       arguments = commandArgs(trailingOnly = TRUE)) {
  output_option <- grep("^--output=", arguments, value = TRUE)
  cores_option <- grep("^--cores=", arguments, value = TRUE)
  unknown <- arguments[!grepl("^(--quick|--output=.+|--cores=[0-9]+)$",
                              arguments)]
  if (length(unknown) || length(output_option) > 1L ||
      length(cores_option) > 1L || sum(arguments == "--quick") > 1L) {
    stop("Use --quick, one --output=directory, and one --cores=positive_integer.",
         call. = FALSE)
  }

  cores <- 20L
  if (length(cores_option)) {
    value <- as.double(sub("^--cores=", "", cores_option))
    if (!is.finite(value) || value < 1 || value > .Machine$integer.max) {
      stop("--cores must be a positive integer.", call. = FALSE)
    }
    cores <- as.integer(value)
  }
  output_directory <- if (length(output_option)) {
    sub("^--output=", "", output_option)
  } else {
    default_output
  }

  list(quick_run = "--quick" %in% arguments,
       output_directory = output_directory, cores = cores)
}

choose_worker_count <- function(requested_cores, replications) {
  available_cores <- parallel::detectCores(logical = TRUE)
  workers <- min(requested_cores, replications)
  if (length(available_cores) == 1L && is.finite(available_cores) &&
      available_cores > 0L) {
    workers <- min(workers, available_cores)
  }
  as.integer(workers)
}

start_simulation_cluster <- function(worker_count, export_names,
                                     envir = parent.frame()) {
  if (worker_count == 1L) return(NULL)

  # Prevent each worker from also requesting a full set of BLAS/OpenMP threads.
  # Set these before starting R workers, and restore the parent environment.
  thread_settings <- c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                       MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
                       BLIS_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")
  old_settings <- Sys.getenv(names(thread_settings), unset = NA_character_)
  on.exit({
    unset <- is.na(old_settings)
    Sys.unsetenv(names(old_settings)[unset])
    if (any(!unset)) do.call(Sys.setenv, as.list(old_settings[!unset]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(thread_settings))

  # Socket workers also work on Windows, where forked workers are unavailable.
  cluster <- parallel::makePSOCKcluster(worker_count)
  initialized <- FALSE
  on.exit({
    if (!initialized) parallel::stopCluster(cluster)
  }, add = TRUE)
  parallel::clusterCall(cluster, function(library_paths, rng_kind) {
    .libPaths(library_paths)
    do.call(RNGkind, as.list(rng_kind))
    suppressPackageStartupMessages({
      library(mecCox)
      library(survival)
    })
    NULL
  }, .libPaths(), RNGkind())
  parallel::clusterExport(cluster, export_names, envir = envir)
  initialized <- TRUE
  cluster
}

run_simulation_replications <- function(replications, replication_function,
                                        arguments, cluster = NULL) {
  # Each scenario sets a seed inside replication_function from the cell and
  # replication indices. Random numbers therefore do not depend on the worker
  # that receives a task, even with load-balanced scheduling.
  run_one <- function(replicate, replication_function, arguments) {
    do.call(replication_function, c(arguments, list(replicate = replicate)))
  }
  if (is.null(cluster)) {
    lapply(seq_len(replications), run_one,
           replication_function = replication_function, arguments = arguments)
  } else {
    parallel::parLapplyLB(cluster, seq_len(replications), run_one,
                          replication_function = replication_function,
                          arguments = arguments, chunk.size = 1L)
  }
}

# Without an output_file, draw on the current graphics device.
plot_scenario1_results <- function(summary, output_file = NULL) {
  methods <- c("Naive", "Robust sandwich", "Corrected sandwich", "MEC-Cox")
  colors <- c("gray25", "#0072B2", "#009E73", "#D55E00")
  line_types <- c(3L, 2L, 4L, 1L)
  symbols <- c(4L, 1L, 2L, 16L)
  metrics <- c(coverage = "Coverage", bias = "Bias", rmse = "RMSE")
  ratios <- unique(summary$ratio)

  if (!is.null(output_file)) {
    grDevices::pdf(output_file, width = 14, height = 3.1 * length(ratios) + 2.5)
  }
  old <- graphics::par(no.readonly = TRUE)
  panel_count <- 3L * length(ratios)
  panel_layout <- matrix(seq_len(panel_count), ncol = 3L, byrow = TRUE)
  panel_layout <- rbind(panel_layout, rep(panel_count + 1L, 3L))
  graphics::layout(panel_layout,
                   heights = c(rep(1, length(ratios)), 0.20))
  graphics::par(mar = c(4, 4.2, 3, 1),
                       oma = c(0, 0, 1.5, 0), las = 1)
  on.exit({
    graphics::par(old)
    if (!is.null(output_file)) grDevices::dev.off()
  }, add = TRUE)
  device_size <- grDevices::dev.size("in")
  if (is.null(output_file) && (device_size[1L] < 5 || device_size[2L] < 4)) {
    # Keep all panels visible in a small RStudio plot pane.
    plot_scale <- min(1, device_size[1L] / 5, device_size[2L] / 4)
    graphics::par(mar = c(3, 3.2, 2.2, 0.6), oma = c(0, 0, 0.8, 0),
                  mgp = c(1.8, 0.55, 0),
                  cex = graphics::par("cex") * plot_scale)
  }
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

plot_scenario2_results <- function(summary, settings, output_file = NULL) {
  methods <- c("Naive", "Robust sandwich", "Corrected sandwich",
               "MEC-Cox (BART/Cox)", "MEC-Cox (BART/RSF)")
  colors <- c("gray25", "#0072B2", "#009E73", "#D55E00", "#CC79A7")
  line_types <- c(3L, 2L, 4L, 1L, 5L)
  symbols <- c(4L, 1L, 2L, 16L, 17L)
  metrics <- c(coverage = "Coverage", bias = "Bias", rmse = "RMSE")

  if (!is.null(output_file)) {
    grDevices::pdf(output_file, width = 14, height = 11.8)
  }
  old <- graphics::par(no.readonly = TRUE)
  panel_layout <- matrix(seq_len(9L), ncol = 3L, byrow = TRUE)
  panel_layout <- rbind(panel_layout, rep(10L, 3L))
  graphics::layout(panel_layout, heights = c(1, 1, 1, 0.25))
  graphics::par(mar = c(4, 4.2, 3, 1),
                       oma = c(0, 0, 1.5, 0), las = 1)
  on.exit({
    graphics::par(old)
    if (!is.null(output_file)) grDevices::dev.off()
  }, add = TRUE)
  device_size <- grDevices::dev.size("in")
  if (is.null(output_file) && (device_size[1L] < 5 || device_size[2L] < 4)) {
    # Keep all panels visible in a small RStudio plot pane.
    plot_scale <- min(1, device_size[1L] / 5, device_size[2L] / 4)
    graphics::par(mar = c(3, 3.2, 2.2, 0.6), oma = c(0, 0, 0.8, 0),
                  mgp = c(1.8, 0.55, 0),
                  cex = graphics::par("cex") * plot_scale)
  }
  panel_index <- 0L

  for (setting_name in settings$setting) {
    for (metric in names(metrics)) {
      panel_index <- panel_index + 1L
      panel_data <- summary[summary$setting == setting_name, , drop = FALSE]
      observed <- panel_data[[metric]][is.finite(panel_data[[metric]])]
      limits <- if (length(observed)) range(observed) else c(0, 1)
      if (diff(limits) < 1e-8) limits <- limits + c(-0.01, 0.01)
      if (metric == "coverage") limits <- range(c(limits, 0.95))
      if (metric == "bias") limits <- range(c(limits, 0))
      graphics::plot(range(panel_data$n1), limits, type = "n",
                     xlab = expression(n[1]), ylab = metrics[metric],
                     main = sprintf("(%s) %s: %s",
                                    letters[panel_index], setting_name,
                                    metrics[metric]))
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
                   horiz = TRUE, bty = "n", cex = 0.8)
}
