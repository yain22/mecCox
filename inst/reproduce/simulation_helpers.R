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
