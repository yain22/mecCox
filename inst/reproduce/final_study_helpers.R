# Prospective final-study adapter. Sourcing this file never starts a fit.
# The final_reference bundle preserves the accepted full-precision inputs.

.final_reference_directory <- local({
  paths <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
  }, character(1))
  paths <- rev(paths[nzchar(paths)])
  paths <- paths[basename(paths) == "final_study_helpers.R"]
  candidates <- c(file.path(dirname(paths), "final_reference"),
                  file.path("inst", "reproduce", "final_reference"),
                  file.path("final_reference"),
                  system.file("reproduce", "final_reference", package = "mecCox"))
  candidates <- candidates[nzchar(candidates) & dir.exists(candidates)]
  if (!length(candidates)) stop("Cannot locate the final_reference bundle.")
  normalizePath(candidates[1L], winslash = "/", mustWork = TRUE)
})

load_final_study_bundle <- function(reference_dir = .final_reference_directory) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Install {digest}.")
  hashes <- utils::read.csv(file.path(reference_dir, "SHA256.csv"), stringsAsFactors = FALSE)
  actual <- vapply(file.path(reference_dir, hashes$file), function(path) {
    digest::digest(file = path, algo = "sha256")
  }, character(1))
  if (!identical(unname(actual), hashes$sha256)) stop("Final-study bundle hash mismatch.")
  read <- function(name) readRDS(file.path(reference_dir, paste0(name, ".rds")))
  list(directory = reference_dir, configs = read("runtime_configs"),
       streams = read("rng_states"), grid = read("job_grid"),
       targets = read("targets"), method_labels = read("method_labels"),
       provenance = read("provenance"))
}

configure_final_study_worker <- function() {
  settings <- c(R_RANGER_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
                OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
                RCPP_PARALLEL_NUM_THREADS = "1")
  previous <- list(environment = Sys.getenv(names(settings), unset = NA_character_),
                   ranger = getOption("ranger.num.threads"))
  do.call(Sys.setenv, as.list(settings))
  options(ranger.num.threads = 1L)
  invisible(previous)
}

restore_final_study_worker <- function(previous) {
  unset <- is.na(previous$environment)
  if (any(unset)) Sys.unsetenv(names(previous$environment)[unset])
  if (any(!unset)) do.call(Sys.setenv, as.list(previous$environment[!unset]))
  options(ranger.num.threads = previous$ranger)
  invisible(NULL)
}

get_final_study_engine <- function(figure, capture_full = FALSE,
                                   reference_dir = .final_reference_directory) {
  figure <- match.arg(figure, c("figure5", "figure6"))
  needed <- c("survival", "MASS", "digest")
  if (figure == "figure6") needed <- c(needed, "dbarts", "ranger")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install: ", paste(missing, collapse = ", "))
  context <- new.env(parent = asNamespace("stats"))
  context$Surv <- survival::Surv
  context$coxph <- survival::coxph
  context$basehaz <- survival::basehaz
  context$head <- utils::head
  context$packageVersion <- utils::packageVersion
  context$main_dir <- reference_dir
  loader <- function(figure, derivative = "analytic") {
    stopifnot(identical(derivative, "analytic"))
    engine <- new.env(parent = environment(main_load_engine))
    sys.source(file.path(main_dir, paste0(figure, "_functions.R")), engine,
               keep.source = FALSE)
    engine
  }
  environment(loader) <- context
  context$main_load_engine <- loader
  # These function bodies are exact extractions from the frozen final run.
  # Only filesystem lookup and the namespace parent are supplied by this adapter.
  sys.source(file.path(reference_dir, "main_instrumentation.R"), context,
             keep.source = FALSE)
  stopifnot(identical(body(context$main_make_engine)[[2L]],
                      quote(ns <- new.env(parent = globalenv()))))
  body(context$main_make_engine)[[2L]] <-
    quote(ns <- new.env(parent = environment(main_make_engine)))
  context$main_make_engine(figure, "analytic", capture_full, TRUE)
}

fit_final_study_job <- function(job_id, bundle = load_final_study_bundle(),
                                capture_full = FALSE) {
  if (length(job_id) != 1L || !is.numeric(job_id) || !is.finite(job_id) ||
      job_id != as.integer(job_id) || !job_id %in% bundle$grid$id) {
    stop("Select one integer job ID from the frozen final-study grid.")
  }
  previous_threads <- configure_final_study_worker()
  on.exit(restore_final_study_worker(previous_threads), add = TRUE)
  job <- bundle$grid[bundle$grid$id == job_id, , drop = FALSE]
  config <- bundle$configs[[job$config_id]]
  engine <- get_final_study_engine(job$figure, capture_full, bundle$directory)
  arguments <- config[intersect(names(config), names(formals(engine$fit_one_rep)))]
  arguments$rep_id <- job$rep_id
  arguments$n1 <- job$n1
  arguments$n0 <- job$n0
  arguments$trim <- config$ps_trim
  arguments$auto_tune_seed <- config$auto_tune_seed + job$rep_id
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  assign(".Random.seed", bundle$streams[[as.character(job_id)]], .GlobalEnv)
  warnings <- character()
  fatal <- NA_character_
  start <- Sys.time()
  result <- tryCatch(withCallingHandlers(do.call(engine$fit_one_rep, arguments),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) { fatal <<- conditionMessage(e); NULL })
  if (is.null(result)) {
    result <- do.call(rbind, lapply(bundle$method_labels[[job$figure]], function(label) {
      engine$empty_row(job$rep_id, job$n1, job$n0, label, NA_real_,
                       paste("Fatal dataset error:", fatal))
    }))
  }
  stopifnot(nrow(result) == length(bundle$method_labels[[job$figure]]),
            all(is.na(result$theta_true)))
  list(job = job, results = result, PS_trace = engine$.PS_trace,
       nuisance_trace = engine$.nuisance_trace, calibration_trace = engine$.calibration_trace,
       estimator_trace = engine$.estimator_trace, point_trace = engine$.point_trace,
       math_trace = engine$.math_trace, solve_trace = engine$.solve_trace,
       OR_tuning_trace = engine$.OR_tuning_trace, data_hash = engine$.data_hash,
       fatal_error = fatal, warnings = warnings,
       R_rng_after = get(".Random.seed", .GlobalEnv, inherits = FALSE),
       elapsed = as.numeric(difftime(Sys.time(), start, units = "secs")),
       source_manifest_md5 = if (!is.null(bundle$provenance$manifest_md5_by_figure))
         bundle$provenance$manifest_md5_by_figure[[job$figure]] else
         bundle$provenance$manifest_md5)
}

bind_final_study_target <- function(run, bundle = load_final_study_bundle()) {
  target <- bundle$targets[bundle$targets$figure == sub("figure", "Figure", run$job$figure) &
                             bundle$targets$setting == run$job$setting, , drop = FALSE]
  stopifnot(nrow(target) == 1L, is.finite(target$numerical_target),
            all(is.na(run$results$theta_true)))
  rows <- run$results
  rows$theta_true <- target$numerical_target
  rows$target_MCSE <- target$mean_MCSE
  rows$legacy_target <- target$legacy_target
  rows
}
