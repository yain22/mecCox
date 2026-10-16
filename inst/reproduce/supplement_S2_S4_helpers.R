# Supplemental Figures S2-S4 reproduction helpers. Sourcing defines functions only.
# Keep the supplemental schema independent of the main-study adapter.

supp234_md5 <- function(path) unname(tools::md5sum(path))
supp234_sha <- function(path) digest::digest(file = path, algo = "sha256")

supp234_check_hashes <- function(directory) {
  h <- utils::read.csv(file.path(directory, "SHA256.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(h) > 0L, !anyDuplicated(h$file), !anyNA(h$file),
            all(nzchar(h$file)), !any(grepl("(^[/\\\\]|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$))", h$file)))
  actual <- vapply(file.path(directory, h$file), supp234_sha, character(1))
  stopifnot(identical(unname(actual), h$sha256))
  invisible(h)
}

load_supp234_bundle <- function(reference_dir) {
  h <- supp234_check_hashes(reference_dir)
  required <- c("runtime_configs.rds", "job_grid.rds", "rng_states.rds", "method_labels.rds",
    "targets.rds", "provenance.rds", "S2_functions.R", "S34_functions.R", "study_instrumentation.R",
    "S2_original_plot.R", "S3_original_plot.R", "S4_original_plot.R",
    "S2_original_plot_dat.R", "S34_original_plot_dat.R")
  stopifnot(all(required %in% h$file))
  read <- function(name) readRDS(file.path(reference_dir, paste0(name, ".rds")))
  b <- list(directory = normalizePath(reference_dir, winslash = "/", mustWork = TRUE),
    configs = read("runtime_configs"), grid = read("job_grid"), streams = read("rng_states"),
    method_labels = read("method_labels"), targets = read("targets"), provenance = read("provenance"))
  stopifnot(identical(b$provenance$schema, "supp_S2_S4_fresh1000_v1"),
    isTRUE(b$provenance$accepted_complete), identical(b$grid$id, seq_len(30000L)),
    length(b$streams) == 30000L, identical(names(b$streams), as.character(b$grid$id)),
    identical(names(b$configs), c("S2_p50", "S2_p200", "S3_ratio2", "S3_ratio4", "S4_ratio2", "S4_ratio4")),
    identical(names(b$method_labels), names(b$configs)),
    all(vapply(b$method_labels, length, integer(1)) == 8L),
    all(table(b$grid$config_id, b$grid$n1) == 1000L),
    identical(b$targets$config_id, names(b$configs)), nrow(b$targets) == 6L,
    all(b$targets$attempted == 100L), all(is.finite(b$targets$numerical_target)))
  b
}

supp234_context <- function() {
  needed <- c("survival", "MASS", "dplyr", "digest")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install: ", paste(absent, collapse = ", "))
  e <- new.env(parent = asNamespace("stats"))
  # Original production workers attach survival, MASS, dplyr, then glmnet.
  # The selected fitting branches call glmnet explicitly via ::.
  for (p in c("survival", "MASS", "dplyr"))
    for (n in getNamespaceExports(p)) assign(n, getExportedValue(p, n), e)
  e$head <- utils::head
  e$tail <- utils::tail
  e$packageVersion <- utils::packageVersion
  e
}

supp234_plain_engine <- function(family, bundle) {
  stopifnot(family %in% c("S2", "S34"))
  e <- supp234_context()
  sys.source(file.path(bundle$directory, paste0(family, "_functions.R")), envir = e, keep.source = FALSE)
  e
}

supp234_check_fit_environment <- function(bundle) {
  d <- bundle$provenance$dependencies
  actual <- vapply(names(d$package_versions), function(p)
    as.character(utils::packageVersion(p)), character(1))
  if (!identical(actual, d$package_versions) || !identical(R.version, d$R) ||
      !identical(Sys.getlocale(), d$locale))
    stop("Exact replay requires the recorded R, package versions and locale; see dependency snapshot.")
  invisible(TRUE)
}

fit_supp234_job <- function(job_id, bundle, capture_full = FALSE) {
  stopifnot(length(job_id) == 1L, is.numeric(job_id), is.finite(job_id),
    job_id == as.integer(job_id), job_id %in% bundle$grid$id)
  supp234_check_fit_environment(bundle)
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install {glmnet}.")
  settings <- c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1", R_RANGER_NUM_THREADS = "1")
  old_env <- Sys.getenv(names(settings), unset = NA_character_)
  old_option <- getOption("ranger.num.threads")
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  on.exit({
    absent <- is.na(old_env)
    if (any(absent)) Sys.unsetenv(names(old_env)[absent])
    if (any(!absent)) do.call(Sys.setenv, as.list(old_env[!absent]))
    options(ranger.num.threads = old_option)
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  do.call(Sys.setenv, as.list(settings))
  options(ranger.num.threads = 1L)
  context <- supp234_context()
  sys.source(file.path(bundle$directory, "study_instrumentation.R"), context, keep.source = FALSE)
  loader <- function(directory, family, derivative = "analytic") {
    stopifnot(identical(derivative, "analytic"), family %in% c("S2", "S34"))
    e <- new.env(parent = environment(study_load_engine))
    sys.source(file.path(directory, paste0(family, "_functions.R")), e, keep.source = FALSE)
    e
  }
  environment(loader) <- context
  context$study_load_engine <- loader
  m <- list(stage = bundle$provenance$stage, configs = bundle$configs, grid = bundle$grid,
    streams = bundle$streams, method_labels = bundle$method_labels,
    source_md5 = bundle$provenance$source_md5)
  # Exact frozen study_one preserves original two-stage rep_id seed nesting,
  # fatal rows, warnings, row order, and data/PS/H/calibration/RNG capture.
  context$study_one(as.integer(job_id), m, bundle$directory,
                    derivative = "analytic", logging = TRUE, capture_full = capture_full)
}

supp234_bind_job <- function(run, bundle) {
  j <- bundle$grid[bundle$grid$id == run$job$id, , drop = FALSE]
  stopifnot(identical(j, run$job), nrow(run$results) == 8L,
    all(is.na(run$results$theta_true)), !anyDuplicated(run$results$Method),
    setequal(run$results$Method, bundle$method_labels[[j$config_id]]))
  r <- run$results
  r$id <- j$id; r$config_id <- j$config_id; r$figure <- j$figure
  r$family <- j$family; r$p <- j$p; r$target_law <- j$target_law
  r$theta_true_raw <- r$theta_true
  r$theta_true <- bundle$targets$numerical_target[match(j$config_id, bundle$targets$config_id)]
  r$fatal_error <- run$fatal_error
  r
}

supp234_check_rows <- function(rows, bundle) {
  needed <- c("id", "config_id", "figure", "family", "p", "target_law", "rep", "n1", "n0",
    "Method", "theta_true_raw", "theta_true", "Estimate", "SE", "CI_L", "CI_U", "Error")
  stopifnot(is.data.frame(rows), all(needed %in% names(rows)), nrow(rows) == 240000L,
    !anyNA(rows$id), !anyNA(rows$Method), !anyDuplicated(rows[c("id", "Method")]),
    setequal(rows$id, bundle$grid$id), all(is.na(rows$theta_true_raw)))
  pos <- match(rows$id, bundle$grid$id)
  for (n in c("config_id", "figure", "family", "p", "target_law", "n1", "n0"))
    stopifnot(identical(unname(rows[[n]]), unname(bundle$grid[[n]][pos])))
  stopifnot(identical(unname(rows$rep), unname(bundle$grid$rep_id[pos])))
  for (k in names(bundle$configs)) {
    r <- rows[rows$config_id == k, , drop = FALSE]
    stopifnot(nrow(r) == 40000L, setequal(r$Method, bundle$method_labels[[k]]),
      all(table(r$n1, r$Method) == 1000L),
      all(r$theta_true == bundle$targets$numerical_target[bundle$targets$config_id == k]))
  }
  invisible(TRUE)
}

supp234_original_group_counts <- function(rows) {
  keys <- c("config_id", "figure", "n1", "n0", "Method")
  stopifnot(all(c(keys, "Estimate", "SE") %in% names(rows)), !anyNA(rows[keys]))
  groups <- unique(rows[keys]); rownames(groups) <- NULL
  key <- function(x) do.call(paste, c(x[keys], sep = "|"))
  index <- match(key(rows), key(groups))
  groups$Attempted <- tabulate(index, nbins = nrow(groups))
  groups$Trials <- tabulate(index[is.finite(rows$Estimate) & is.finite(rows$SE)],
                           nbins = nrow(groups))
  groups
}

supp234_check_original_groups <- function(summary, rows, counts = NULL) {
  if (is.null(counts)) counts <- supp234_original_group_counts(rows)
  # The original summarizer filters before grouping: a zero-eligible group is
  # absent, not an invented row of metrics. All attempts remain in rows/counts.
  eligible <- counts[counts$Trials > 0L, , drop = FALSE]
  supp234_same_table(summary, eligible, c("config_id", "n1", "Method"),
    c("config_id", "figure", "n1", "n0", "Method", "Trials"))
  invisible(counts)
}

supp234_original_summary <- function(rows, bundle) {
  supp234_check_rows(rows, bundle)
  summaries <- lapply(names(bundle$configs), function(k) {
    family <- if (startsWith(k, "S2")) "S2" else "S34"
    e <- supp234_plain_engine(family, bundle)
    # Deliberately call the exact original summarizer, including its finite
    # Estimate/SE denominator and NA semantics. Ordered-CI is diagnostic only.
    s <- as.data.frame(e$summarize_simulation(rows[rows$config_id == k, , drop = FALSE]))
    data.frame(config_id = rep.int(k, nrow(s)),
      figure = rep.int(substr(k, 1L, 2L), nrow(s)), s, check.names = FALSE)
  })
  ans <- do.call(rbind, summaries); rownames(ans) <- NULL
  counts <- supp234_original_group_counts(rows)
  stopifnot(nrow(counts) == 240L, all(counts$Attempted == 1000L))
  supp234_check_original_groups(ans, rows, counts)
  ans
}

supp234_same_table <- function(actual, expected, keys, columns = names(expected)) {
  key <- function(x) do.call(paste, c(x[keys], sep = "|"))
  a <- key(actual); b <- key(expected)
  stopifnot(!anyDuplicated(a), !anyDuplicated(b), setequal(a, b),
    all(columns %in% names(actual)), all(columns %in% names(expected)))
  expected <- expected[match(a, b), , drop = FALSE]
  for (n in columns) {
    if (!identical(unname(actual[[n]]), unname(expected[[n]])))
      stop("Exact table comparison failed for column: ", n)
  }
  invisible(TRUE)
}

read_supp234_reference <- function(path, bundle) {
  h <- supp234_check_hashes(dirname(path))
  stopifnot(basename(path) %in% h$file)
  x <- readRDS(path)
  stopifnot(isTRUE(x$metadata$complete),
    identical(x$metadata$schema, "supp_S2_S4_fresh1000_v1"),
    identical(x$metadata$manifest_md5, bundle$provenance$manifest_md5),
    identical(x$metadata$target_binding_md5, bundle$provenance$target_binding_md5),
    identical(x$metadata$acceptance_md5, bundle$provenance$acceptance_md5),
    identical(x$targets, bundle$targets))
  supp234_check_rows(x$results, bundle)
  supp234_same_table(supp234_original_summary(x$results, bundle), x$summary,
                    c("config_id", "n1", "Method"))
  x
}

supp234_plot_locale <- function() {
  previous <- Sys.getlocale("LC_CTYPE")
  if (!isTRUE(l10n_info()[["UTF-8"]])) {
    # Windows locale verified on the study platform. No numerical/fitting locale
    # is changed; callers use this only for plotting and restore it on exit.
    requested <- if (.Platform$OS.type == "windows") "English_United States.utf8" else "C.UTF-8"
    applied <- Sys.setlocale("LC_CTYPE", requested)
    if (!nzchar(applied) || !isTRUE(l10n_info()[["UTF-8"]])) {
      Sys.setlocale("LC_CTYPE", previous)
      stop("Plotting requires an available UTF-8 LC_CTYPE; no locale fallback is used.")
    }
  }
  list(previous = previous, active = Sys.getlocale("LC_CTYPE"), UTF8 = TRUE)
}

supp234_restore_plot_locale <- function(context) {
  if (!identical(Sys.getlocale("LC_CTYPE"), context$previous))
    Sys.setlocale("LC_CTYPE", context$previous)
  stopifnot(identical(Sys.getlocale("LC_CTYPE"), context$previous))
  invisible(TRUE)
}

supp234_check_plot_labels <- function(style, methods) {
  stopifnot(is.character(methods), length(methods) > 0L, !anyNA(methods),
    !anyDuplicated(methods), all(methods %in% style$method_order),
    all(methods %in% unname(style$method_labels)))
  for (name in c("method_colors", "method_shapes", "method_linetypes")) {
    values <- style[[name]]
    position <- match(methods, names(values))
    stopifnot(!anyDuplicated(names(values)), !anyNA(position), !anyNA(values[position]))
  }
  invisible(TRUE)
}
supp234_plot_style <- function(figure, bundle) {
  plot_locale <- supp234_plot_locale()
  on.exit(supp234_restore_plot_locale(plot_locale), add = TRUE)
  stopifnot(figure %in% c("S2", "S3", "S4"))
  versions <- c(ggplot2 = "4.0.0", patchwork = "1.3.2", gtable = "0.3.6")
  actual <- vapply(names(versions), function(p) as.character(utils::packageVersion(p)), character(1))
  stopifnot(identical(actual, versions))
  for (p in c("ggplot2", "patchwork", "gtable", "dplyr", "labeling")) loadNamespace(p)
  e <- new.env(parent = asNamespace("ggplot2"))
  e$unit <- grid::unit; e$filter <- dplyr::filter; e$plot_layout <- patchwork::plot_layout
  keep <- c("method_labels", "method_order", "method_colors", "method_shapes",
    "method_linetypes", "theme_pub", "make_panel", "build_2x3")
  found <- character()
  for (x in parse(file.path(bundle$directory, paste0(figure, "_original_plot.R")), keep.source = FALSE, encoding = "UTF-8")) {
    if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) &&
        as.character(x[[2L]]) %in% keep) {
      found <- c(found, as.character(x[[2L]])); eval(x, e)
    }
  }
  stopifnot(setequal(found, keep), !anyDuplicated(found))
  supp234_check_plot_labels(e, e$method_order)
  e
}

supp234_original_plot_dat <- function(rows, summary, cfg, bundle) {
  versions <- c(dplyr = "1.1.4", tidyr = "1.3.1")
  actual <- vapply(names(versions), function(p) as.character(utils::packageVersion(p)), character(1))
  stopifnot(identical(actual, versions), identical(bundle$provenance$plot_data_dependencies, versions))
  family <- unique(rows$family)
  stopifnot(length(family) == 1L)
  e <- supp234_plain_engine(family, bundle)
  e$res <- rows
  e$summ <- summary[setdiff(names(summary), c("config_id", "figure"))]
  e$Lin_Wei_ps_method <- cfg$Lin_Wei_ps_method
  e$include_unweighted_cox <- FALSE
  e$include_mec_no_cf <- FALSE
  # Exact five original assignments retain NA method/size plotting rows while
  # deriving sample-size levels only from the original eligible summary.
  sys.source(file.path(bundle$directory, paste0(family, "_original_plot_dat.R")),
             e, keep.source = FALSE)
  list(data = e$plot_dat, method_map = e$method_map)
}

build_supp234_figures <- function(reference, bundle) {
  summary <- supp234_original_summary(reference$results, bundle)
  supp234_same_table(summary, reference$summary, c("config_id", "n1", "Method"))
  plot_locale <- supp234_plot_locale()
  on.exit(supp234_restore_plot_locale(plot_locale), add = TRUE)
  figures <- plot_data <- list()
  for (fig in c("S2", "S3", "S4")) {
    style <- supp234_plot_style(fig, bundle)
    keys <- names(bundle$configs)[startsWith(names(bundle$configs), fig)]
    labels <- if (fig == "S2") c("M = 50", "M = 200") else c("n1:n0 = 1:2", "n1:n0 = 1:4")
    methods <- c("Naive", "Robust sandwich", "Corrected sandwich", "MEC-Cox (KL)")
    if (fig == "S3") methods <- c(methods, "MEC-Cox (Hellinger)", "MEC-Cox (EL)", "MEC-Cox (R\u00e9nyi)")
    supp234_check_plot_labels(style, methods)
    pieces <- lapply(seq_along(keys), function(i) {
      k <- keys[i]; cfg <- bundle$configs[[k]]
      r <- reference$results[reference$results$config_id == k, , drop = FALSE]
      s <- summary[summary$config_id == k, , drop = FALSE]
      original <- supp234_original_plot_dat(r, s, cfg, bundle)
      map <- original$method_map
      d <- as.data.frame(original$data)
      d$config_id <- rep.int(k, nrow(d))
      d$figure <- rep.int(fig, nrow(d))
      d$Method_Label <- factor(unname(style$method_labels[as.character(d$Method_Short)]), levels = style$method_order)
      stopifnot(!anyNA(d$Method_Label))
      d <- d[d$Method_Label %in% methods, , drop = FALSE]
      d$Method_Label <- droplevels(d$Method_Label)
      d$.setting <- factor(rep.int(labels[i], nrow(d)), levels = labels)
      displayed_methods <- map$Method[
        unname(style$method_labels[as.character(map$Method_Short)]) %in% methods]
      counts <- supp234_original_group_counts(r)
      counts <- counts[counts$Method %in% displayed_methods, , drop = FALSE]
      stopifnot(nrow(counts) == (if (fig == "S3") 35L else 20L),
        all(counts$Attempted == 1000L))
      # Plot rows include all displayed methods at summary-present sizes. A
      # joined NA row is plotting structure, never a fabricated summary result.
      visible <- counts[paste(counts$n1, counts$n0) %in% paste(s$n1, s$n0), , drop = FALSE]
      checked <- d; checked$Method <- as.character(checked$Method)
      supp234_same_table(checked, visible, c("config_id", "n1", "Method"),
        c("config_id", "figure", "n1", "n0", "Method"))
      supp234_check_original_groups(checked[!is.na(checked$Trials), , drop = FALSE], r, visible)
      metrics <- setdiff(names(s), c("config_id", "figure", "n1", "n0", "Method"))
      stopifnot(all(vapply(d[is.na(d$Trials), metrics, drop = FALSE],
        function(x) all(is.na(x)), logical(1))))
      d[order(d$Method_Label, d$n1), , drop = FALSE]
    })
    d <- do.call(rbind, pieces); rownames(d) <- NULL
    figures[[fig]] <- style$build_2x3(d, ".setting", setting_titles = labels, include_methods = methods)
    plot_data[[fig]] <- d
  }
  list(figures = figures, plot_data = plot_data)
}

supp234_new_directory <- function(path) {
  if (file.exists(path)) stop("Choose a new output directory; existing paths are never replaced.")
  if (!dir.create(path, recursive = TRUE)) stop("Cannot create output directory.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

supp234_write_hashes <- function(directory) {
  files <- sort(setdiff(list.files(directory), "SHA256.csv"))
  stopifnot(!any(dir.exists(file.path(directory, files))))
  h <- data.frame(file = files, sha256 = vapply(file.path(directory, files), supp234_sha, character(1)), row.names = NULL)
  utils::write.csv(h, file.path(directory, "SHA256.csv"), row.names = FALSE)
  invisible(h)
}
