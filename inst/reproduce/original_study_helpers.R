# Adapters for the main-paper simulation studies.
# Source this file; it does not load data or launch simulation jobs.
# Study functions are isolated from the package API and the caller's globals.

.original_study_directory <- local({
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
  }, character(1))
  files <- rev(files[nzchar(files)])
  files <- files[basename(files) == "original_study_helpers.R"]
  candidates <- c(
    file.path(dirname(files), "study_reference"),
    file.path("inst", "reproduce", "study_reference"),
    file.path("study_reference"),
    system.file("reproduce", "study_reference", package = "mecCox")
  )
  candidates <- candidates[nzchar(candidates) & dir.exists(candidates)]
  if (!length(candidates)) {
    stop("Cannot find the study_reference directory beside this helper.",
         call. = FALSE)
  }
  normalizePath(candidates[1L], winslash = "/", mustWork = TRUE)
})

original_study_metadata <- function(scenario,
                                    reference_dir = .original_study_directory) {
  scenario <- match.arg(as.character(scenario), c("scenario1", "scenario2"))
  metadata <- dget(file.path(reference_dir, "metadata.R"))
  metadata[[scenario]]
}

get_original_study_engine <- function(scenario,
                                      reference_dir = .original_study_directory) {
  scenario <- match.arg(as.character(scenario), c("scenario1", "scenario2"))
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Install {survival} to run the study engine.", call. = FALSE)
  }
  # Explicit bindings resolve stats/survival functions independently of
  # user-defined functions and the mecCox API.
  engine <- new.env(parent = asNamespace("stats"))
  engine$Surv <- survival::Surv
  engine$coxph <- survival::coxph
  engine$basehaz <- survival::basehaz
  engine$head <- utils::head
  sys.source(file.path(reference_dir, paste0(scenario, "_functions.R")),
             envir = engine, keep.source = FALSE)
  attr(engine, "study_reference") <- original_study_metadata(scenario, reference_dir)
  engine
}

get_original_study_config <- function(scenario, setting,
                                      reference_dir = .original_study_directory) {
  scenario <- match.arg(as.character(scenario), c("scenario1", "scenario2"))
  setting <- as.character(setting)
  if (length(setting) != 1L || is.na(setting)) {
    stop("Supply one study setting.", call. = FALSE)
  }
  if (scenario == "scenario1" && setting %in% c("2", "3", "4")) {
    setting <- paste0("ratio", setting)
  }
  configurations <- dget(file.path(reference_dir, "saved_configs.R"))[[scenario]]
  if (!setting %in% names(configurations)) {
    stop("Unknown study setting; choose ",
         paste(names(configurations), collapse = ", "), ".", call. = FALSE)
  }
  configurations[[setting]]
}

original_study_streams <- function(config, replicates = config$R) {
  if (!requireNamespace("rngtools", quietly = TRUE)) {
    stop("Install {rngtools} to generate the study RNG streams.",
         call. = FALSE)
  }
  if (!is.numeric(replicates) || length(replicates) != 1L ||
      !is.finite(replicates) || replicates < 1 || replicates != floor(replicates)) {
    stop("`replicates` must be a positive integer.", call. = FALSE)
  }
  sizes <- config$sample_size_grid
  if (!is.data.frame(sizes) || !all(c("n1", "n0") %in% names(sizes)) ||
      !nrow(sizes) || any(!is.finite(as.matrix(sizes[, c("n1", "n0")])))) {
    stop("The saved configuration lacks a valid sample-size grid.", call. = FALSE)
  }
  if (!is.numeric(config$seed) || length(config$seed) != 1L ||
      !is.finite(config$seed)) {
    stop("The saved configuration lacks a valid seed.", call. = FALSE)
  }
  # expand.grid() varies sample-size index fastest. Generate streams
  # on the entire saved grid before selecting jobs for a smaller quick run.
  grid <- expand.grid(ss = seq_len(nrow(sizes)), rep_id = seq_len(replicates))
  grid$n1 <- sizes$n1[grid$ss]
  grid$n0 <- sizes$n0[grid$ss]
  grid$job_id <- seq_len(nrow(grid))
  grid <- grid[, c("job_id", "ss", "rep_id", "n1", "n0")]

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  streams <- rngtools::RNGseq(nrow(grid), seed = config$seed,
                             simplify = FALSE, version = 2)
  list(grid = grid, streams = streams)
}

# This wrapper is self-contained so it can be exported alone to PSOCK workers.
# fit_one_rep uses separate nuisance fits, fixed fold/forest seeds,
# and replicate-specific tuning seeds.
fit_original_study_rep <- function(engine, config, n1, rep_id, rng_stream) {
  if (!is.environment(engine) ||
      !exists("fit_one_rep", envir = engine, inherits = FALSE)) {
    stop("`engine` must come from get_original_study_engine().", call. = FALSE)
  }
  if (!is.numeric(n1) || length(n1) != 1L || !is.finite(n1) ||
      !n1 %in% config$sample_size_grid$n1) {
    stop("`n1` must be one of the original sample sizes.", call. = FALSE)
  }
  if (!is.numeric(rep_id) || length(rep_id) != 1L || !is.finite(rep_id) ||
      rep_id < 1 || rep_id != floor(rep_id)) {
    stop("`rep_id` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(rng_stream) || length(rng_stream) != 7L ||
      anyNA(rng_stream) || any(!is.finite(rng_stream)) ||
      any(rng_stream != as.integer(rng_stream)) || rng_stream[1L] %% 100L != 7L) {
    stop("`rng_stream` must be a seven-integer L'Ecuyer-CMRG state.",
         call. = FALSE)
  }
  arguments <- config[intersect(names(config), names(formals(engine$fit_one_rep)))]
  arguments$rep_id <- as.integer(rep_id)
  arguments$n1 <- n1
  arguments$n0 <- config$sample_size_grid$n0[match(n1, config$sample_size_grid$n1)]
  arguments$trim <- config$ps_trim
  # Add rep_id at the job level; fit_one_rep also increments the tuning seed
  # when passing it to MEC.
  arguments$auto_tune_seed <- config$auto_tune_seed + rep_id

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  assign(".Random.seed", as.integer(rng_stream), envir = .GlobalEnv)
  do.call(engine$fit_one_rep, arguments)
}
