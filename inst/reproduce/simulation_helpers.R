# Execution helpers shared by the two simulation scripts.
# Original study generators and fitting functions are loaded by the scenario
# scripts; this file handles execution and presentation.

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

stop_simulation_cluster <- function(cluster) {
  if (is.null(cluster)) return(invisible(NULL))
  # One broken socket must not prevent the other owned sockets from closing,
  # or replace an interrupt with a second error during stack unwinding.
  for (node in cluster) {
    single_node <- structure(list(node), class = class(cluster))
    tryCatch(suppressWarnings(parallel::stopCluster(single_node)),
             error = function(condition) NULL,
             interrupt = function(condition) NULL)
    if (!is.null(node$con)) {
      tryCatch(suppressWarnings(close(node$con)),
               error = function(condition) NULL,
               interrupt = function(condition) NULL)
    }
  }
  invisible(NULL)
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
    if (!initialized) stop_simulation_cluster(cluster)
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
                                        arguments, cluster = NULL,
                                        on_result = NULL) {
  # Each scenario restores its run-specific random stream inside
  # replication_function. Scheduling therefore does not change the stream
  # received by a dataset.
  run_one <- function(replicate, replication_function, arguments) {
    do.call(replication_function, c(arguments, list(replicate = replicate)))
  }
  environment(run_one) <- baseenv()
  if (!is.null(on_result)) {
    if (!is.function(on_result)) {
      stop("on_result must be a function or NULL.", call. = FALSE)
    }
    results <- vector("list", replications)
    record_result <- function(replicate, value) {
      results[replicate] <<- list(value)
      on_result(replicate, value)
    }
    if (is.null(cluster)) {
      for (replicate in seq_len(replications)) {
        value <- run_one(replicate, replication_function, arguments)
        record_result(replicate, value)
      }
    } else {
      run_simulation_socket_completions(
        cluster, replications, run_one, replication_function, arguments,
        record_result
      )
    }
    return(results)
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

run_simulation_socket_completions <- function(cluster, replications, run_one,
                                               replication_function, arguments,
                                               on_result) {
  # Isolate the two parallel internals needed for completion callbacks. The
  # public parLapplyLB() only exposes results after the entire batch finishes.
  send_call <- utils::getFromNamespace("sendCall", "parallel")
  receive_result <- utils::getFromNamespace("recvOneResult", "parallel")
  submit <- function(node, replicate) {
    send_call(cluster[[node]], run_one,
              list(replicate, replication_function, arguments), tag = replicate)
  }
  active <- min(length(cluster), replications)
  for (node in seq_len(active)) submit(node, node)
  for (completed in seq_len(replications)) {
    received <- receive_result(cluster)
    if (inherits(received$value, "try-error")) {
      stop(sprintf("Replication %d failed on worker %d: %s",
                   received$tag, received$node, as.character(received$value)),
           call. = FALSE)
    }
    # Persist the completed result before another task is submitted. If the
    # callback is interrupted, earlier results remain available to the driver.
    on_result(received$tag, received$value)
    next_replicate <- active + completed
    if (next_replicate <= replications) {
      submit(received$node, next_replicate)
    }
  }
  invisible(NULL)
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
  execution <- results$metadata$execution
  partial <- !is.null(execution$status) && execution$status != "completed"
  if (partial) title <- paste(title, "(partial results)")
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
      if (partial) sprintf(
        "The run was %s; %d of %d planned datasets are included. ",
        execution$status, execution$completed_replications,
        execution$planned_replications
      ) else "The complete aggregated results are shown below. ",
      "Individual simulation results remain available in the R session."
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

# Publication plotting adapters for the main-paper simulations.
# The plot design follows the original revised publication scripts.
check_publication_plot_packages <- function() {
  packages <- c("ggplot2", "patchwork", "labeling")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install the plotting packages: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
}

build_publication_simulation_figure <- function(data, method_order, colors,
                                                shapes, line_types,
                                                legend_rows = 1L) {
  check_publication_plot_packages()
  theme_pub <- ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      legend.key.width = grid::unit(1.5, "cm"),
      strip.text = ggplot2::element_text(size = 13, face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0)
    )
  data$Method_Label <- factor(data$Method_Label, levels = method_order)
  data <- data[!is.na(data$Method_Label), , drop = FALSE]
  data$Method_Label <- droplevels(data$Method_Label)
  colors <- stats::setNames(colors, method_order)
  shapes <- stats::setNames(shapes, method_order)
  line_types <- stats::setNames(line_types, method_order)
  make_panel <- function(data, metric, ylab, title, hline = NULL) {
    p <- ggplot2::ggplot(data, ggplot2::aes(
      x = n1, y = .data[[metric]], color = Method_Label,
      shape = Method_Label, linetype = Method_Label
    )) +
      ggplot2::geom_line(linewidth = 0.6) +
      ggplot2::geom_point(size = 2.2) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_shape_manual(values = shapes) +
      ggplot2::scale_linetype_manual(values = line_types) +
      ggplot2::labs(x = expression(n[1]), y = ylab, title = title) + theme_pub
    if (!is.null(hline)) {
      p <- p + ggplot2::geom_hline(
        yintercept = hline, linetype = "dashed", color = "gray30", linewidth = 0.5
      ) + ggplot2::scale_y_continuous(breaks = function(lims) {
        b <- labeling::extended(lims[1], lims[2], m = 5)
        sort(unique(c(b, hline)))
      })
    }
    p
  }
  settings <- levels(data$.setting)
  panels <- vector("list", length(settings) * 3L)
  for (i in seq_along(settings)) {
    d <- data[data$.setting == settings[i], , drop = FALSE]
    first <- (i - 1L) * 3L + 1L
    panels[[first]] <- make_panel(
      d, "coverage", "Coverage",
      paste0("(", letters[first], ") Coverage: ", settings[i]), 0.95
    )
    panels[[first + 1L]] <- make_panel(
      d, "bias", "Bias",
      paste0("(", letters[first + 1L], ") Bias: ", settings[i]), 0
    )
    panels[[first + 2L]] <- make_panel(
      d, "rmse", "RMSE",
      paste0("(", letters[first + 2L], ") RMSE: ", settings[i])
    )
  }
  rows <- lapply(seq_along(settings), function(i) {
    first <- (i - 1L) * 3L + 1L
    panels[[first]] + panels[[first + 1L]] + panels[[first + 2L]]
  })
  figure <- Reduce(`/`, rows) + patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")
  if (legend_rows > 1L) {
    figure <- figure & ggplot2::guides(
      color = ggplot2::guide_legend(nrow = legend_rows, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = legend_rows, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = legend_rows, byrow = TRUE)
    )
  }
  figure
}

build_scenario1_figure <- function(summary) {
  if (!nrow(summary)) return(NULL)
  methods <- c("Naive", "Robust sandwich", "Corrected sandwich",
               "MEC-Cox (KL-GLM(PS)-Cox(OR))")
  labels <- stats::setNames(methods, c("Naive", "Robust sandwich",
                                      "Corrected sandwich", "MEC-Cox"))
  data <- summary
  data$Method_Label <- unname(labels[as.character(data$method)])
  ratio_order <- unique(as.character(data$ratio))
  data$.setting <- factor(paste0("n1:n0 = ", data$ratio),
                          levels = paste0("n1:n0 = ", ratio_order))
  build_publication_simulation_figure(
    data, methods, c("gray25", "#0072B2", "#009E73", "#D55E00"),
    c(4, 1, 2, 16), c("dotted", "dashed", "longdash", "solid")
  )
}

build_scenario2_figure <- function(summary, settings) {
  settings <- settings[settings$setting %in% summary$setting, , drop = FALSE]
  if (!nrow(settings) || !nrow(summary)) return(NULL)
  methods <- c("Naive", "Robust sandwich", "Corrected sandwich",
               "MEC-Cox (KL-BART(PS)-Cox(OR))", "MEC-Cox (KL-BART(PS)-RSF(OR))")
  labels <- stats::setNames(methods, c("Naive", "Robust sandwich",
                                      "Corrected sandwich", "MEC-Cox (BART/Cox)",
                                      "MEC-Cox (BART/RSF)"))
  data <- summary
  data$Method_Label <- unname(labels[as.character(data$method)])
  data$.setting <- factor(paste0("Nonlinearity: ", tolower(data$setting)),
                          levels = paste0("Nonlinearity: ", tolower(settings$setting)))
  build_publication_simulation_figure(
    data, methods, c("gray25", "#0072B2", "#009E73", "#D55E00", "blue"),
    c(4, 1, 2, 16, 17), c("dotted", "dashed", "longdash", "solid", "solid"),
    legend_rows = 2L
  )
}

# Draw only on an already available graphics device; explicit exports are below.
plot_scenario1_results <- function(summary) {
  if (!simulation_graphics_available()) return(invisible(FALSE))
  figure <- build_scenario1_figure(summary)
  if (is.null(figure)) return(invisible(FALSE))
  print(figure)
  invisible(figure)
}

plot_scenario2_results <- function(summary, settings) {
  if (!simulation_graphics_available()) return(invisible(FALSE))
  figure <- build_scenario2_figure(summary, settings)
  if (is.null(figure)) return(invisible(FALSE))
  print(figure)
  invisible(figure)
}

save_simulation_figure <- function(figure, filename) {
  check_publication_plot_packages()
  ggplot2::ggsave(filename, figure, width = 12, height = 9, units = "in", dpi = 300)
  invisible(filename)
}
