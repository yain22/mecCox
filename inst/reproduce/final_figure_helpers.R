# Complete-only final-study comparison and explicit figure export.
# Sourcing this file starts no fits, reads no outcomes, and writes no files.
.final_figure_directory <- local({
  paths <- vapply(sys.frames(), function(frame)
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L], character(1))
  paths <- rev(paths[nzchar(paths)])
  paths <- paths[basename(paths) == "final_figure_helpers.R"]
  candidates <- c(dirname(paths), file.path("inst", "reproduce"), ".",
                   system.file("reproduce", package = "mecCox"))
  candidates <- candidates[file.exists(file.path(candidates, "final_study_helpers.R")) &
                             dir.exists(file.path(candidates, "final_plot_reference"))]
  if (!length(candidates)) stop("Cannot locate the final reproduction files.")
  normalizePath(candidates[1L], winslash = "/", mustWork = TRUE)
})

final_figure_bundle <- function(directory = .final_figure_directory) {
  api <- new.env(parent = globalenv())
  source(file.path(directory, "final_study_helpers.R"), local = api)
  api$load_final_study_bundle(file.path(directory, "final_reference"))
}

final_figure_roles <- function(labels) {
  roles <- as.character(labels)
  mec <- grepl("^MEC-Cox:", roles)
  roles[mec & grepl("OR=cox", labels, fixed = TRUE)] <- "MEC_Cox"
  roles[mec & grepl("OR=rsf", labels, fixed = TRUE)] <- "MEC_RSF"
  allowed <- c("ATT-IPW Cox: Naive model-based", "ATT-IPW Cox: Lin-Wei",
                "ATT-IPW Cox: Shu", "MEC_Cox", "MEC_RSF")
  if (anyNA(roles) || any(!roles %in% allowed)) stop("Unrecognized actual method label.")
  roles
}

final_figure_check_rows <- function(rows, bundle, figures = c("figure5", "figure6")) {
  required <- c("id", "figure", "setting", "config_id", "rep", "n1", "n0",
                 "Method", "role", "theta_true", "Estimate", "SE", "CI_L", "CI_U",
                 "Error", "target_MCSE", "legacy_target")
  if (!is.data.frame(rows) || !all(required %in% names(rows)))
    stop("Incomplete final result schema.")
  grid <- bundle$grid[bundle$grid$figure %in% figures, , drop = FALSE]
  if (!length(figures) || any(!figures %in% c("figure5", "figure6")) ||
      anyDuplicated(figures) || nrow(grid) != 15000L * length(figures))
    stop("Unexpected final-study figure selection.")
  expected_roles <- lapply(bundle$method_labels, final_figure_roles)
  expected_count <- sum(vapply(grid$figure, function(f) length(expected_roles[[f]]), integer(1)))
  if (nrow(rows) != expected_count || !setequal(rows$id, grid$id) ||
      anyNA(rows$id) || anyDuplicated(paste(rows$id, rows$role)))
    stop("Incomplete or duplicated final-study rows.")
  position <- match(rows$id, grid$id)
  for (column in c("figure", "setting", "config_id", "n1", "n0")) {
    if (!identical(unname(rows[[column]]), unname(grid[[column]][position])))
      stop("Misaligned frozen job metadata: ", column)
  }
  if (!identical(unname(rows$rep), unname(grid$rep_id[position])) ||
      !identical(unname(rows$role), unname(final_figure_roles(rows$Method))))
    stop("Misaligned replication or method role.")
  groups <- split(rows$role, rows$id)
  group_figures <- grid$figure[match(as.integer(names(groups)), grid$id)]
  if (!all(vapply(seq_along(groups), function(i) {
    setequal(groups[[i]], expected_roles[[group_figures[i]]])
  }, logical(1)))) stop("An attempted dataset is missing a method row.")
  target_index <- match(paste(rows$figure, rows$setting),
                         paste(tolower(bundle$targets$figure), bundle$targets$setting))
  if (anyNA(target_index)) stop("Unknown target mapping.")
  for (pair in list(c("theta_true", "numerical_target"), c("target_MCSE", "mean_MCSE"),
                    c("legacy_target", "legacy_target"))) {
    if (!identical(unname(rows[[pair[1L]]]),
                    unname(bundle$targets[[pair[2L]]][target_index])))
      stop("Wrong final target policy: ", pair[1L])
  }
  counts <- table(paste(rows$config_id, rows$n1, rows$role))
  if (any(counts != 1000L)) stop("Every final method/cell requires all 1000 attempts.")
  invisible(TRUE)
}

final_figure_primary_summary <- function(rows, bundle = final_figure_bundle()) {
  figures <- unique(rows$figure)
  final_figure_check_rows(rows, bundle, figures)
  groups <- split(seq_len(nrow(rows)), paste(rows$config_id, rows$n1, rows$role))
  out <- lapply(groups, function(index) {
    x <- rows[index, , drop = FALSE]
    finite <- is.finite(x$Estimate) & is.finite(x$SE)
    interval <- finite & is.finite(x$CI_L) & is.finite(x$CI_U) & x$CI_L <= x$CI_U
    delta <- x$Estimate[finite] - x$theta_true[finite]
    covered <- sum(interval & x$CI_L <= x$theta_true & x$CI_U >= x$theta_true)
    denominator <- sum(interval)
    data.frame(figure = x$figure[1L], setting = x$setting[1L],
      config_id = x$config_id[1L], n1 = x$n1[1L], n0 = x$n0[1L], role = x$role[1L],
      target_type = "independent_primary", eligibility = "legacy_finite",
      theta = x$theta_true[1L], target_MCSE = x$target_MCSE[1L],
      attempted = nrow(x), valid = sum(finite), legacy_nonfinite = sum(!finite),
      coverage_numerator = covered, coverage_denominator = denominator,
      coverage = if (denominator) covered / denominator else NA_real_,
      valid_and_covers_per_attempt = covered / nrow(x),
      bias = if (length(delta)) mean(delta) else NA_real_,
      RMSE = if (length(delta)) sqrt(mean(delta^2)) else NA_real_,
      errors = sum(!is.na(x$Error) & nzchar(x$Error)), stringsAsFactors = FALSE)
  })
  answer <- do.call(rbind, out)
  rownames(answer) <- NULL
  answer
}

final_figure_compare_tables <- function(actual, expected, key_columns, value_columns,
                                         tolerance = 0) {
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0)
    stop("Use an explicit nonnegative absolute comparison tolerance.")
  key <- function(x) do.call(paste, c(x[key_columns], sep = "|"))
  akey <- key(actual); ekey <- key(expected)
  if (anyDuplicated(akey) || anyDuplicated(ekey) || !setequal(akey, ekey))
    stop("Comparison keys are missing, duplicated or misaligned.")
  expected <- expected[match(akey, ekey), , drop = FALSE]
  differences <- list()
  for (column in value_columns) {
    a <- actual[[column]]; e <- expected[[column]]
    if (is.null(a) || is.null(e)) stop("Missing comparison column: ", column)
    equal <- is.na(a) & is.na(e)
    if (is.numeric(a) && is.numeric(e)) {
      # Distinguish missing NA/NaN, finite values, and signed infinities.
      equal <- equal & (is.nan(a) == is.nan(e))
      finite <- is.finite(a) & is.finite(e)
      equal[finite] <- abs(a[finite] - e[finite]) <= tolerance
      nonfinite <- !is.na(a) & !is.na(e) & !finite
      equal[nonfinite] <- a[nonfinite] == e[nonfinite]
    } else {
      observed <- !is.na(a) & !is.na(e)
      equal[observed] <- as.character(a[observed]) == as.character(e[observed])
    }
    bad <- which(!equal)
    if (length(bad)) differences[[length(differences) + 1L]] <- data.frame(
      key = akey[bad], column = column, actual = as.character(a[bad]),
      expected = as.character(e[bad]), stringsAsFactors = FALSE)
  }
  if (length(differences)) do.call(rbind, differences) else
    data.frame(key = character(), column = character(), actual = character(),
               expected = character(), stringsAsFactors = FALSE)
}

final_figure_bind_source_hashes <- function(gate_sources, manifest, bundle_sources) {
  # Private audit keys and public basenames intentionally differ. Bind the exact
  # ordered mapping through the audited manifest rather than discarding names.
  paths <- manifest$source_files
  if (!is.character(paths) || is.null(names(paths)) || anyDuplicated(names(paths)) ||
      !identical(names(paths), names(gate_sources)) ||
      !identical(gate_sources, manifest$source_md5))
    stop("Private source identities do not match the audited manifest.")
  canonical <- stats::setNames(unname(gate_sources), basename(gsub("\\\\", "/", paths)))
  if (!identical(canonical, bundle_sources))
    stop("Private-to-public source identity mapping differs from the frozen bundle.")
  canonical
}

read_audited_final_collector <- function(results_path, gate_path, target_binding_path,
                                          bundle = final_figure_bundle()) {
  md5 <- function(path) unname(tools::md5sum(path))
  gate <- readRDS(gate_path)
  if (!isTRUE(gate$PASS) || !identical(gate$results_md5, md5(results_path)) ||
      !identical(gate$target_binding_md5, md5(target_binding_path)) ||
      !identical(gate$manifest_md5, bundle$provenance$manifest_md5))
    stop("The complete independent audit does not bind these final inputs.")
  input <- readRDS(results_path)
  bound_sources <- final_figure_bind_source_hashes(gate$source_md5, input$manifest,
                                                  bundle$provenance$source_md5)
  if (!isTRUE(input$metadata$complete) || input$metadata$reps_per_cell != 1000L ||
      !identical(input$metadata$manifest_md5, gate$manifest_md5) ||
      !identical(input$metadata$collector_md5, gate$collector_md5) ||
      !identical(input$metadata$target_binding_md5, gate$target_binding_md5) ||
      !identical(input$metadata$primary_target, "independent_primary") ||
      !identical(input$target_binding, readRDS(target_binding_path)) ||
      !identical(input$target_binding$policy_md5, bundle$provenance$target_policy_md5))
    stop("Final collection, numerical policy or target binding is incomplete.")
  if (!"theta_true_raw" %in% names(input$results) ||
      !all(is.na(input$results$theta_true_raw))) stop("Raw target provenance changed.")
  final_figure_check_rows(input$results, bundle)
  summary <- final_figure_primary_summary(input$results, bundle)
  expected <- input$performance[input$performance$target_type == "independent_primary" &
                                 input$performance$eligibility == "legacy_finite", , drop = FALSE]
  if (nrow(expected) != 135L || nrow(summary) != 135L)
    stop("Missing final primary summary cells.")
  delta <- final_figure_compare_tables(summary, expected,
    c("figure", "setting", "n1", "role"), setdiff(names(summary),
      c("figure", "setting", "n1", "role")), tolerance = 1e-13)
  if (nrow(delta)) stop("Raw final values disagree with the audited primary summary.")
  list(input = input, summary = summary, audit = gate, bound_source_md5 = bound_sources,
       input_hashes = c(results = md5(results_path), audit = md5(gate_path),
                        target_binding = md5(target_binding_path)))
}

final_figure_new_directory <- function(path) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)))
    stop("Choose an empty output directory; existing artifacts are never replaced.")
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) stop("Cannot create output directory.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

final_figure_hash_files <- function(directory, files) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Install {digest}.")
  data.frame(file = files, sha256 = vapply(file.path(directory, files), function(path)
    digest::digest(file = path, algo = "sha256"), character(1)), row.names = NULL)
}

final_figure_portable_source_hashes <- function(hashes) {
  if (is.null(names(hashes)) || any(!nzchar(names(hashes))))
    stop("Source hash identities are missing.")
  identities <- basename(gsub("\\\\", "/", names(hashes)))
  # The immutable source vector can contain files with the same basename from
  # separate folders. Its frozen position supplies an unambiguous portable ID.
  identities <- sprintf("source_%02d__%s", seq_along(identities), identities)
  stats::setNames(unname(hashes), identities)
}

export_audited_final_reference <- function(results_path, gate_path, target_binding_path,
                                            output_directory, bundle = final_figure_bundle()) {
  if (dir.exists(output_directory) &&
      length(list.files(output_directory, all.files = TRUE, no.. = TRUE)))
    stop("Choose an empty output directory; existing artifacts are never replaced.")
  checked <- read_audited_final_collector(results_path, gate_path, target_binding_path, bundle)
  input <- checked$input
  out <- final_figure_new_directory(output_directory)
  # Carry all model rows and reporting flags, without private local filesystem paths.
  raw_columns <- c("id", "figure", "setting", "config_id", "rep", "n1", "n0", "Method",
    "role", "theta_true", "theta_true_raw", "target_MCSE", "legacy_target", "Estimate", "SE",
    "CI_L", "CI_U", "ESS", "Rel_ESS", "W_CV", "W_Min", "W_Max", "Cal_Grad", "Cal_Converged",
    "Error", "fatal_error", "raw_sandwich_SE", "SE_fallback_used", "KL_positive_finite_weights",
    "KL_exponential_weights_identical")
  if (!all(raw_columns %in% names(input$results))) stop("Missing retained reporting fields.")
  audit_fields <- c("PASS", "no_fits", "manifest_md5", "results_md5", "target_binding_md5",
    "collector_md5", "attempted_datasets", "attempted_per_cell", "method_rows", "summary_rows",
    "all_checkpoint_rows_and_diagnostics_verified", "PS_duplicates_not_pooled",
    "all_failures_extremes_preserved", "strict_domain_fallback_and_legacy_denominators_verified")
  portable_audit <- checked$audit[intersect(audit_fields, names(checked$audit))]
  portable_audit$source_md5 <- final_figure_portable_source_hashes(checked$bound_source_md5)
  reference <- list(results = input$results[raw_columns], summary = checked$summary,
    targets = bundle$targets,
    metadata = list(complete = TRUE, replications_per_cell = 1000L,
      manifest_md5 = bundle$provenance$manifest_md5,
      target_policy_md5 = bundle$provenance$target_policy_md5,
      target_type = "independent_primary", eligibility = "legacy_finite",
      input_hashes = checked$input_hashes, independent_complete_audit = portable_audit,
      raw_summary_absolute_tolerance = 1e-13,
      source_stage = bundle$provenance$stage))
  saveRDS(reference, file.path(out, "final_figure_reference.rds"))
  utils::write.csv(reference$summary, file.path(out, "primary_summary.csv"), row.names = FALSE)
  utils::write.csv(input$performance, file.path(out, "all_target_and_eligibility_summaries.csv"), row.names = FALSE)
  utils::write.csv(reference$results, file.path(out, "all_attempted_model_rows.csv"), row.names = FALSE)
  utils::write.csv(input$calibration_status, file.path(out, "raw_calibration_status.csv"), row.names = FALSE)
  utils::write.csv(input$targets[intersect(names(bundle$targets), names(input$targets))],
                   file.path(out, "numerical_targets.csv"), row.names = FALSE)
  files <- list.files(out)
  utils::write.csv(final_figure_hash_files(out, files), file.path(out, "SHA256.csv"), row.names = FALSE)
  invisible(reference)
}

read_final_figure_reference <- function(path, bundle = final_figure_bundle()) {
  hashes <- utils::read.csv(file.path(dirname(path), "SHA256.csv"), stringsAsFactors = FALSE)
  selected <- hashes[hashes$file == basename(path), , drop = FALSE]
  actual <- final_figure_hash_files(dirname(path), basename(path))
  if (nrow(selected) != 1L || !identical(selected$sha256, actual$sha256))
    stop("Portable final reference hash mismatch.")
  x <- readRDS(path)
  if (!isTRUE(x$metadata$complete) || x$metadata$replications_per_cell != 1000L ||
      !identical(x$metadata$manifest_md5, bundle$provenance$manifest_md5) ||
      !identical(x$metadata$target_policy_md5, bundle$provenance$target_policy_md5) ||
      !identical(x$metadata$target_type, "independent_primary") ||
      !identical(x$metadata$eligibility, "legacy_finite") ||
      !identical(x$targets, bundle$targets) || !isTRUE(x$metadata$independent_complete_audit$PASS) ||
      !identical(x$metadata$independent_complete_audit$manifest_md5, x$metadata$manifest_md5) ||
      !identical(x$metadata$independent_complete_audit$source_md5,
                  final_figure_portable_source_hashes(bundle$provenance$source_md5)) ||
      !identical(x$metadata$independent_complete_audit$results_md5,
                  unname(x$metadata$input_hashes[["results"]])) ||
      !identical(x$metadata$independent_complete_audit$target_binding_md5,
                  unname(x$metadata$input_hashes[["target_binding"]])))
    stop("Wrong or incomplete final reference provenance.")
  final_figure_check_rows(x$results, bundle)
  rebuilt <- final_figure_primary_summary(x$results, bundle)
  delta <- final_figure_compare_tables(rebuilt, x$summary,
    c("figure", "setting", "n1", "role"), setdiff(names(rebuilt),
      c("figure", "setting", "n1", "role")), tolerance = 1e-13)
  if (nrow(delta)) stop("Portable primary summary does not match all raw rows.")
  x
}

final_driver_rows_for_comparison <- function(result, figure, bundle) {
  if (!identical(result$status, "completed") ||
      result$metadata$execution$completed_replications != 15000L ||
      result$metadata$execution$planned_replications != 15000L ||
      result$metadata$design$replications != 1000L ||
      !identical(result$metadata$study_engine$manifest_md5, bundle$provenance$manifest_md5) ||
      !identical(result$metadata$target_policy, bundle$targets))
    stop("Comparison requires the complete current final-driver run.")
  d <- result$replications
  if (!"final_study_job_id" %in% names(d) || !"original_method" %in% names(d))
    stop("Missing frozen job IDs or actual method labels.")
  grid <- bundle$grid[match(d$final_study_job_id, bundle$grid$id), , drop = FALSE]
  if (anyNA(grid$id) || any(grid$figure != figure)) stop("Wrong driver grid.")
  rows <- data.frame(id = d$final_study_job_id, figure = grid$figure, setting = grid$setting,
    config_id = grid$config_id, rep = d$replicate, n1 = d$n1, n0 = d$n0,
    Method = d$original_method, role = final_figure_roles(d$original_method),
    theta_true = d$target, theta_true_raw = NA_real_, Estimate = d$estimate,
    SE = d$standard_error, CI_L = d$ci_lower, CI_U = d$ci_upper, Error = d$error,
    stringsAsFactors = FALSE)
  extra <- c("target_MCSE", "legacy_target", "ESS", "Rel_ESS", "W_CV", "W_Min", "W_Max",
    "Cal_Grad", "Cal_Converged", "fatal_error", "raw_sandwich_SE", "SE_fallback_used",
    "KL_positive_finite_weights", "KL_exponential_weights_identical")
  if (!all(extra %in% names(d))) stop("Missing final-driver diagnostic fields.")
  rows[extra] <- d[extra]
  final_figure_check_rows(rows, bundle, figure)
  # Require each attempted raw capture, including failed runs. Do not infer raw
  # target-NA provenance merely from a synthetic column in the display table.
  traces <- result$diagnostics
  if (length(traces) != 15000L) stop("Complete driver diagnostics are missing.")
  ids <- vapply(traces, function(x) x$job$id, integer(1))
  if (anyDuplicated(ids) || !setequal(ids, grid$id)) stop("Misaligned driver diagnostics.")
  row_index <- split(seq_len(nrow(rows)), rows$id)
  for (trace in traces) {
    x <- rows[row_index[[as.character(trace$job$id)]], , drop = FALSE]
    raw <- trace$results
    if (!all(is.na(raw$theta_true)) || anyNA(match(x$Method, raw$Method)))
      stop("Changed raw target or actual method provenance.")
    raw <- raw[match(x$Method, raw$Method), , drop = FALSE]
    for (column in c("rep", "n1", "n0", "Method", "Estimate", "SE", "CI_L", "CI_U", "Error")) {
      if (!identical(unname(x[[column]]), unname(raw[[column]])))
        stop("Driver table disagrees with raw capture: ", column)
    }
  }
  rows
}

compare_final_driver_results <- function(scenario1_results, scenario2_results,
                                           reference_path, report_directory,
                                           tolerance = 0, bundle = final_figure_bundle()) {
  reference <- read_final_figure_reference(reference_path, bundle)
  actual <- rbind(final_driver_rows_for_comparison(scenario1_results, "figure5", bundle),
                  final_driver_rows_for_comparison(scenario2_results, "figure6", bundle))
  final_figure_check_rows(actual, bundle)
  summary <- final_figure_primary_summary(actual, bundle)
  row_differences <- final_figure_compare_tables(actual, reference$results,
    c("id", "role"), setdiff(names(reference$results), c("id", "role")), tolerance)
  summary_differences <- final_figure_compare_tables(summary, reference$summary,
    c("figure", "setting", "n1", "role"), setdiff(names(summary),
      c("figure", "setting", "n1", "role")), tolerance)
  out <- final_figure_new_directory(report_directory)
  utils::write.csv(row_differences, file.path(out, "model_row_mismatches.csv"), row.names = FALSE)
  utils::write.csv(summary_differences, file.path(out, "primary_summary_mismatches.csv"), row.names = FALSE)
  verdict <- list(PASS = nrow(row_differences) == 0L && nrow(summary_differences) == 0L,
    datasets = 30000L, model_rows = 135000L, summary_cells = 135L,
    absolute_tolerance = tolerance, target_type = "independent_primary", eligibility = "legacy_finite",
    reference_sha256 = final_figure_hash_files(dirname(reference_path), basename(reference_path))$sha256,
    manifest_md5 = bundle$provenance$manifest_md5,
    row_mismatches = nrow(row_differences), summary_mismatches = nrow(summary_differences))
  saveRDS(verdict, file.path(out, "complete_numeric_comparison.rds"))
  if (!verdict$PASS) stop("Final reproduction differs; all mismatches saved in ", out)
  invisible(verdict)
}

final_figure_check_render_versions <- function() {
  expected <- c(ggplot2 = "4.0.0", patchwork = "1.3.2", gtable = "0.3.6")
  needed <- c(names(expected), "dplyr", "labeling", "digest")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the rendering packages: ", paste(missing, collapse = ", "))
  actual <- vapply(names(expected), function(p) as.character(utils::packageVersion(p)), character(1))
  if (!identical(actual, expected)) stop(
    "Use a separate rendering library with ggplot2 4.0.0, patchwork 1.3.2 and gtable 0.3.6.")
  invisible(actual)
}

load_original_final_plot_style <- function(figure, directory = .final_figure_directory) {
  figure <- match.arg(figure, c("figure5", "figure6"))
  final_figure_check_render_versions()
  reference <- file.path(directory, "final_plot_reference")
  name <- paste0(figure, "_original_revised.R")
  hashes <- utils::read.csv(file.path(reference, "SHA256.csv"), stringsAsFactors = FALSE)
  wanted <- hashes[hashes$file == name, , drop = FALSE]
  if (nrow(wanted) != 1L || !identical(wanted$sha256,
      final_figure_hash_files(reference, name)$sha256)) stop("Original plotting source hash mismatch.")
  keep <- c("method_labels", "method_order", "method_colors", "method_shapes",
             "method_linetypes", "theme_pub", "make_panel", "build_Nx3")
  e <- new.env(parent = asNamespace("ggplot2"))
  e$unit <- grid::unit
  e$filter <- dplyr::filter
  e$plot_layout <- patchwork::plot_layout
  found <- character()
  for (expression in parse(file.path(reference, name), keep.source = FALSE)) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in% keep) {
      found <- c(found, as.character(expression[[2L]]))
      eval(expression, envir = e)
    }
  }
  if (!setequal(found, keep) || anyDuplicated(found)) stop("Incomplete original style assignments.")
  e
}

build_exact_final_figures <- function(reference, bundle = final_figure_bundle(),
                                        directory = .final_figure_directory) {
  final_figure_check_rows(reference$results, bundle)
  data <- final_figure_primary_summary(reference$results, bundle)
  plots <- plot_data <- list()
  for (figure in c("figure5", "figure6")) {
    e <- load_original_final_plot_style(figure, directory)
    d <- data[data$figure == figure, , drop = FALSE]
    role_order <- c("ATT-IPW Cox: Naive model-based", "ATT-IPW Cox: Lin-Wei",
                    "ATT-IPW Cox: Shu", "MEC_Cox", if (figure == "figure6") "MEC_RSF")
    labels <- stats::setNames(e$method_order, role_order)
    settings <- if (figure == "figure5") c("ratio2", "ratio3", "ratio4") else c("none", "mild", "severe")
    titles <- if (figure == "figure5") paste0("n1:n0 = 1:", 2:4) else paste0("Nonlinearity: ", settings)
    d$Method_Label <- factor(unname(labels[d$role]), levels = e$method_order)
    d$.setting <- factor(d$setting, levels = settings)
    if (anyNA(d$Method_Label) || anyNA(d$.setting)) stop("Unmapped final plotting label.")
    d$Coverage_95 <- d$coverage
    d$Bias <- d$bias
    d <- d[order(d$.setting, d$Method_Label, d$n1), , drop = FALSE]
    plots[[figure]] <- e$build_Nx3(d, "setting", setting_titles = titles,
                                  include_methods = e$method_order)
    plot_data[[figure]] <- d
  }
  list(plots = plots, plot_data = plot_data)
}

redraw_final_figures <- function(reference_path, output_directory,
                                   bundle = final_figure_bundle(), directory = .final_figure_directory) {
  reference <- read_final_figure_reference(reference_path, bundle)
  rendered <- build_exact_final_figures(reference, bundle, directory)
  out <- final_figure_new_directory(output_directory)
  for (figure in names(rendered$plots)) {
    stem <- sub("figure", "Figure", figure)
    for (format in c("png", "pdf")) ggplot2::ggsave(
      file.path(out, paste0(stem, "_final_fresh1000.", format)), rendered$plots[[figure]],
      width = 12, height = 9, units = "in", dpi = 300, bg = "white")
    utils::write.csv(rendered$plot_data[[figure]],
                     file.path(out, paste0(stem, "_plot_data.csv")), row.names = FALSE)
  }
  metadata <- list(complete = TRUE, datasets = 30000L, attempts_per_cell = 1000L,
    target_type = "independent_primary", eligibility = "legacy_finite",
    reference_sha256 = final_figure_hash_files(dirname(reference_path), basename(reference_path))$sha256,
    source_manifest_md5 = bundle$provenance$manifest_md5, width_in = 12, height_in = 9,
    dpi = 300, background = "white", session = capture.output(utils::sessionInfo()),
    render_versions = final_figure_check_render_versions(),
    graphics_capabilities = capabilities(), pdf_device_defaults = grDevices::pdf.options(),
    ragg_version = if (requireNamespace("ragg", quietly = TRUE)) as.character(utils::packageVersion("ragg")) else NA_character_,
    style_hashes = utils::read.csv(file.path(directory, "final_plot_reference", "SHA256.csv")),
    helper_sha256 = final_figure_hash_files(directory, "final_figure_helpers.R")$sha256,
    note = "Exact original style assignments; all failures and finite extremes retained. Visual QA remains required.")
  saveRDS(metadata, file.path(out, "render_metadata.rds"))
  utils::write.csv(final_figure_hash_files(out, list.files(out)),
                   file.path(out, "SHA256.csv"), row.names = FALSE)
  invisible(rendered)
}
