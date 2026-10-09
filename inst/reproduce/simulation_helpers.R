# Execution helpers shared by the two simulation reproduction scripts.
# The scientific generators and fitting calls remain in each scenario script.

parse_simulation_arguments <- function(arguments = character(),
                                       quick_run = FALSE, cores = 20L,
                                       replications = 1000L) {
  cores_option <- grep("^--cores=", arguments, value = TRUE)
  replications_option <- grep("^--replications=", arguments, value = TRUE)
  unknown <- arguments[!grepl("^(--quick|--cores=[0-9]+|--replications=[0-9]+)$",
                              arguments)]
  if (length(unknown) || length(replications_option) > 1L ||
      length(cores_option) > 1L || sum(arguments == "--quick") > 1L) {
    stop("Use --quick, one --cores=positive_integer, and ",
         "one --replications=positive_integer. Results are not saved to files.",
         call. = FALSE)
  }
  if (!is.logical(quick_run) || length(quick_run) != 1L || is.na(quick_run)) {
    stop("quick_run must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(cores_option)) {
    cores <- as.double(sub("^--cores=", "", cores_option))
  }
  cores <- simulation_positive_integer(cores, "cores")
  explicit_replications <- length(replications_option) > 0L
  if (explicit_replications) {
    replications <- as.double(sub("^--replications=", "", replications_option))
  }
  replications <- simulation_positive_integer(replications, "replications")
  quick_run <- quick_run || "--quick" %in% arguments
  # Keep --quick short unless the user has chosen a different count in the
  # configuration block or has supplied an explicit command-line count.
  if (quick_run && !explicit_replications && replications == 1000L) {
    replications <- 2L
  }

  list(quick_run = quick_run, cores = cores, replications = replications)
}

simulation_positive_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value != floor(value) || value > .Machine$integer.max) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  as.integer(value)
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

check_simulation_display_packages <- function() {
  packages <- c("kableExtra", "htmltools", "rstudioapi")
  installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (any(!installed)) {
    stop("The simulation tables require: ",
         paste(packages[!installed], collapse = ", "),
         ". Install these packages before starting the simulation.",
         call. = FALSE)
  }
  invisible(TRUE)
}

format_simulation_table <- function(data, title) {
  check_simulation_display_packages()
  table <- kableExtra::kbl(
    data, format = "html", digits = 4, row.names = FALSE,
    caption = title, col.names = gsub("_", " ", names(data), fixed = TRUE),
    escape = TRUE
  )
  table <- kableExtra::kable_styling(
    table, bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE, position = "left", font_size = 14
  )
  table <- kableExtra::row_spec(table, 0L, bold = TRUE,
                               color = "white", background = "#285579")
  if (nrow(data) > 20L) {
    table <- kableExtra::scroll_box(table, width = "100%", height = "520px")
  } else {
    table <- kableExtra::scroll_box(table, width = "100%")
  }
  table
}

build_simulation_report <- function(results, title) {
  check_simulation_display_packages()
  if (is.data.frame(results$targets)) {
    targets <- results$targets
    if ("target" %in% names(targets)) targets$target_hr <- exp(targets$target)
  } else {
    targets <- data.frame(target_log_hr = results$metadata$target,
                           target_hr = exp(results$metadata$target))
  }
  configuration_values <- c(list(quick_run = results$metadata$quick_run),
                             results$metadata$design,
                             results$metadata$execution)
  configuration <- data.frame(
    setting = names(configuration_values),
    value = vapply(configuration_values, function(value) {
      paste(as.character(value), collapse = ", ")
    }, character(1)),
    stringsAsFactors = FALSE
  )
  tables <- list(
    summary = format_simulation_table(results$summary,
                                      paste(title, "simulation summary")),
    targets = format_simulation_table(targets,
                                      paste(title, "reference target(s)")),
    configuration = format_simulation_table(configuration,
                                             paste(title, "run configuration"))
  )
  if (is.data.frame(results$metadata$settings)) {
    tables$settings <- format_simulation_table(results$metadata$settings,
                                                paste(title, "settings"))
  }
  content <- lapply(tables, function(table) {
    htmltools::tags$section(htmltools::HTML(as.character(table)))
  })
  report <- htmltools::browsable(htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$title(title),
      htmltools::tags$style(htmltools::HTML(paste(
        "body { font-family: Arial, sans-serif; color: #222; padding: 20px; }",
        "h1 { font-size: 24px; margin-bottom: 8px; }",
        "section { margin: 26px 0; }",
        "table { border-collapse: collapse; margin-bottom: 12px; }",
        "caption { font-size: 17px; font-weight: bold; text-align: left;",
        "padding: 10px 0; color: #222; }",
        "th, td { padding: 8px 12px; border-bottom: 1px solid #ddd;",
        "white-space: nowrap; }",
        "tbody tr:nth-child(even) { background: #f5f7fa; }",
        "tbody tr:hover { background: #eaf1f7; }"
      )))
    ),
    htmltools::tags$h1(title),
    htmltools::tags$p(
      "The complete aggregated results are shown below. Individual ",
      "replication records remain available in the R session."
    ),
    content
  ))
  list(tables = tables, report = report)
}

show_simulation_report <- function(report, viewer = NULL) {
  if (is.null(viewer)) {
    if (!interactive()) return(invisible(FALSE))
    viewer <- getOption("viewer")
    if (!is.function(viewer) && rstudioapi::isAvailable()) {
      viewer <- rstudioapi::viewer
    }
    if (!is.function(viewer)) {
      message("The formatted tables are retained in the *_tables and ",
              "*_report objects. Use RStudio to display them in its Viewer.")
      return(invisible(FALSE))
    }
  }
  if (!is.function(viewer)) {
    stop("viewer must be a function or NULL.", call. = FALSE)
  }
  # html_print uses a temporary HTML document for the Viewer, not a results
  # folder. The simulation does not write CSV, RDS, or PDF output.
  htmltools::html_print(report, viewer = viewer)
  invisible(TRUE)
}

simulation_graphics_available <- function() {
  if (grDevices::dev.cur() != 1L || interactive()) return(TRUE)
  message("Plots require an active graphics device. Source this script in ",
          "RStudio or a graphical R session to display them; no PDF is saved.")
  FALSE
}

# Draw on the active R graphics device; never open a file graphics device.
plot_scenario1_results <- function(summary) {
  if (!simulation_graphics_available()) return(invisible(FALSE))
  methods <- c("Naive", "Robust sandwich", "Corrected sandwich", "MEC-Cox")
  colors <- c("gray25", "#0072B2", "#009E73", "#D55E00")
  line_types <- c(3L, 2L, 4L, 1L)
  symbols <- c(4L, 1L, 2L, 16L)
  metrics <- c(coverage = "Coverage", bias = "Bias", rmse = "RMSE")
  ratios <- unique(summary$ratio)

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
  }, add = TRUE)
  device_size <- grDevices::dev.size("in")
  if (device_size[1L] < 5 || device_size[2L] < 4) {
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

plot_scenario2_results <- function(summary, settings) {
  if (!simulation_graphics_available()) return(invisible(FALSE))
  methods <- c("Naive", "Robust sandwich", "Corrected sandwich",
               "MEC-Cox (BART/Cox)", "MEC-Cox (BART/RSF)")
  colors <- c("gray25", "#0072B2", "#009E73", "#D55E00", "#CC79A7")
  line_types <- c(3L, 2L, 4L, 1L, 5L)
  symbols <- c(4L, 1L, 2L, 16L, 17L)
  metrics <- c(coverage = "Coverage", bias = "Bias", rmse = "RMSE")

  old <- graphics::par(no.readonly = TRUE)
  panel_layout <- matrix(seq_len(9L), ncol = 3L, byrow = TRUE)
  panel_layout <- rbind(panel_layout, rep(10L, 3L))
  graphics::layout(panel_layout, heights = c(1, 1, 1, 0.25))
  graphics::par(mar = c(4, 4.2, 3, 1),
                       oma = c(0, 0, 1.5, 0), las = 1)
  on.exit({
    graphics::par(old)
  }, add = TRUE)
  device_size <- grDevices::dev.size("in")
  if (device_size[1L] < 5 || device_size[2L] < 4) {
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
