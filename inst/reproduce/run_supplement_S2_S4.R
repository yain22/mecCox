# Supplemental Figures S2-S4 explicit selected-job runner.
# Source supplement_S2_S4_helpers.R first, then call run_supplement_S2_S4
# with explicit job_ids and a new output directory. Sourcing runs no jobs.
# The runner is serial and does not recalculate target roots.

run_supplement_S2_S4 <- function(job_ids, reference_dir, output_directory) {
  bundle <- load_supp234_bundle(reference_dir)
  stopifnot(is.numeric(job_ids), length(job_ids) > 0L, !anyNA(job_ids),
    !anyDuplicated(job_ids), all(is.finite(job_ids)), all(job_ids == as.integer(job_ids)),
    all(job_ids %in% bundle$grid$id))
  # Fail before creating files if the exact replay environment is unavailable.
  supp234_check_fit_environment(bundle)
  out <- supp234_new_directory(output_directory)
  records <- vector("list", length(job_ids))
  for (i in seq_along(job_ids)) {
    records[[i]] <- fit_supp234_job(job_ids[i], bundle)
    saveRDS(records[[i]], file.path(out, sprintf("dataset_%05d.rds", job_ids[i])), version = 3)
  }
  rows <- do.call(rbind, lapply(records, supp234_bind_job, bundle = bundle))
  rownames(rows) <- NULL
  complete <- length(job_ids) == 30000L && setequal(job_ids, bundle$grid$id)
  summary <- if (complete) supp234_original_summary(rows, bundle) else NULL
  result <- list(status = if (complete) "completed" else "selected_jobs_only",
    job_ids = job_ids, records = records, results = rows, summary = summary,
    manifest_md5 = bundle$provenance$manifest_md5,
    targets = bundle$targets, target_binding_md5 = bundle$provenance$target_binding_md5,
    note = "All attempted method rows retained. Selected jobs do not reproduce full-study figures.")
  saveRDS(result, file.path(out, "simulation_results.rds"), version = 3)
  supp234_write_hashes(out)
  invisible(result)
}
