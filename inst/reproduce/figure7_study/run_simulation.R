# Run from this directory after the independent DGP and estimator pilots pass:
# Rscript run_simulation.R 1000 10
# This file orchestrates the established estimator; it does not modify estimation.
source("dgp.R")
source("estimator.R")

simulation_seed <- function(replication, master_seed = 20261106L) {
  stopifnot(replication >= 1L, replication <= 9999L)
  # HR is deliberately absent so null and alternative use paired random draws.
  as.integer(master_seed + 2000000L + replication)
}

verify_original_source <- function() {
  audit <- read.csv("project_code_verification.csv", stringsAsFactors = FALSE)
  stopifnot(nrow(audit) >= 1L, all(audit$unchanged),
    length(unique(audit$source_md5)) == 1L,
    length(unique(audit$extracted_md5)) == 1L,
    identical(unique(audit$extracted_md5),
              unname(tools::md5sum("original_project_functions.R"))))
  invisible(TRUE)
}

validate_replication_rows <- function(rows, target_hr) {
  required <- c("method", "PS", "OR", "theta", "HR", "true_HR", "n1", "n0",
                "converged", "reason")
  methods <- c("Naive", "IPW", "MEC1", "MEC2", "MEC5", "MEC10", "MEC20", "Oracle")
  stopifnot(is.data.frame(rows), all(required %in% names(rows)), nrow(rows) == 32L,
    setequal(rows$method, methods),
    all(rows$PS %in% c("Correct", "Misspecified")),
    all(rows$OR %in% c("Correct", "Misspecified")),
    all(table(rows$PS, rows$OR, rows$method) == 1L),
    all(rows$true_HR == target_hr), all(rows$n1 == 200L), all(rows$n0 == 400L))
  invisible(TRUE)
}

failed_replication_rows <- function(target_hr, error_message) {
  # Preserve a complete failed replicate if nuisance fitting raises an exception.
  rows <- expand.grid(PS = c("Correct", "Misspecified"),
    OR = c("Correct", "Misspecified"),
    method = c("Naive", "IPW", "MEC1", "MEC2", "MEC5", "MEC10", "MEC20", "Oracle"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  rows$theta <- rows$HR <- NA_real_
  rows$true_HR <- target_hr
  rows$n1 <- 200L
  rows$n0 <- 400L
  rows$converged <- FALSE
  rows$reason <- paste("Replication-level estimation exception:", error_message)
  rows$warning <- ""
  rows$calibration_converged <- NA
  rows$calibration_grad_norm <- NA_real_
  rows$original_mean_imbalance <- NA_real_
  rows$calibration_iterations <- NA_integer_
  rows$calibration_coordinates <- ""
  rows$calibration_initial_error <- ""
  rows$calibration_kkt_residual <- NA_real_
  rows
}

run_one_dataset <- function(job, master_seed, fold_seed) {
  seed <- simulation_seed(job$replication, master_seed)
  config <- dgp_config(hr = job$true_HR)
  set.seed(seed)
  data <- dgp_draw_data(n0 = 400L, n1 = 200L, config = config)
  stopifnot(nrow(data) == 600L, sum(data$A == 0L) == 400L,
            sum(data$A == 1L) == 200L)
  rows <- tryCatch(mec_fit_replication(data, config, fold_seed = fold_seed, Kfold = 10L),
    error = function(e) failed_replication_rows(job$true_HR, conditionMessage(e)))
  validate_replication_rows(rows, job$true_HR)
  rows$replication <- as.integer(job$replication)
  rows$seed <- seed
  rows$design <- "original_workflow_mixed5_5"
  rows$total_n <- 600L
  rows
}

save_new_rds <- function(object, path) {
  # A pre-existing checkpoint is never overwritten.
  stopifnot(!file.exists(path))
  temporary <- paste0(path, ".tmp")
  saveRDS(object, temporary)
  stopifnot(file.rename(temporary, path))
  invisible(path)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  B <- if (length(args)) as.integer(args[1]) else 1000L
  workers <- if (length(args) > 1L) as.integer(args[2]) else 10L
  stopifnot(B >= 1L, B <= 9999L, workers >= 1L, workers <= 10L)
  verify_original_source()
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  dir.create("results/checkpoints", recursive = TRUE, showWarnings = FALSE)
  cfg <- dgp_config(hr = 0.70)
  # Only numerical settings are serialized; interpolation closures are rebuilt.
  cfg$transform <- list(nodes_per_interval = cfg$transform$nodes_per_interval,
                        grid_size = cfg$transform$grid_size)
  source_files <- c("dgp.R", "estimator.R", "original_project_functions.R",
                    "run_simulation.R", "project_code_verification.csv")
  specification <- list(
    B = B, n1 = 200L, n0 = 400L, sizes = 200L, total_n = 600L,
    true_HR = c(1, 0.70), designs = "original_workflow_mixed5_5",
    master_seed = 20261106L, fold_seed = 20260427L, crossfit_folds = 10L,
    source_md5 = tools::md5sum(source_files), configuration = cfg,
    landmark_counts = c(1L, 2L, 5L, 10L, 20L),
    single_landmark_probability = 0.10, landmark_quantile_type = 8L,
    oracle_landmarks = 5L,
    oracle = "MEC with true source odds and true conditional control survival at five landmarks",
    misspecified_PS_variables = "X1", misspecified_OR_variables = "X1",
    calibration_fallback = "Full square invertible SVD; all original constraints retained",
    fallback_original_mean_tolerance = 1e-10, fallback_KL_KKT_tolerance = 1e-6,
    R_version = R.version.string,
    survival_version = as.character(utils::packageVersion("survival")))
  manifest <- "results/specification.rds"
  if (file.exists(manifest)) {
    stopifnot(identical(readRDS(manifest), specification))
  } else save_new_rds(specification, manifest)
  tasks <- expand.grid(replication = seq_len(B), true_HR = specification$true_HR,
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  groups <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / 25L))
  cl <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, {
    Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    source("dgp.R")
    source("estimator.R")
    NULL
  })
  parallel::clusterExport(cl, c("simulation_seed", "validate_replication_rows", "failed_replication_rows", "run_one_dataset"),
                          envir = environment(run_one_dataset))
  started <- proc.time()[3]
  for (g in seq_along(groups)) {
    checkpoint <- sprintf("results/checkpoints/batch_%03d.rds", g)
    if (!file.exists(checkpoint)) {
      jobs <- lapply(groups[[g]], function(i) as.list(tasks[i, ]))
      fits <- parallel::parLapplyLB(cl, jobs, run_one_dataset,
        master_seed = specification$master_seed, fold_seed = specification$fold_seed)
      batch <- do.call(rbind, fits)
      stopifnot(nrow(batch) == 32L * length(jobs))
      save_new_rds(batch, checkpoint)
    } else {
      batch <- readRDS(checkpoint)
      expected <- tasks[groups[[g]], , drop = FALSE]
      stopifnot(nrow(batch) == 32L * nrow(expected),
        setequal(unique(batch$replication), expected$replication),
        setequal(unique(batch$true_HR), expected$true_HR))
    }
    cat(sprintf("Completed %d/%d datasets; elapsed %.1f seconds.\n",
                max(groups[[g]]), nrow(tasks), proc.time()[3] - started))
    flush.console()
  }
  results <- do.call(rbind, lapply(seq_along(groups), function(g)
    readRDS(sprintf("results/checkpoints/batch_%03d.rds", g))))
  rownames(results) <- NULL
  stopifnot(nrow(results) == 32L * nrow(tasks),
    !anyDuplicated(results[, c("true_HR", "replication", "PS", "OR", "method")]),
    all(with(results, table(true_HR, PS, OR, method)) == B))
  output <- "results/replications.rds"
  if (file.exists(output)) stopifnot(identical(readRDS(output), results)) else
    save_new_rds(results, output)
  cat(sprintf("Finished: %d datasets, %d fits, %d unconverged fits.\n",
              nrow(tasks), nrow(results), sum(!results$converged)))
}

if (sys.nframe() == 0L) main()
