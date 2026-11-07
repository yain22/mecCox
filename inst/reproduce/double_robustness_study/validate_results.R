# Integrity and audit information; this is separate from the displayed boxplots.
source("run_simulation.R")
spec <- readRDS("results/specification.rds")
raw <- readRDS("results/replications.rds")
verify_original_source()
stopifnot(identical(spec$source_md5, tools::md5sum(names(spec$source_md5))),
  spec$n1 == 200L, spec$n0 == 400L, spec$total_n == 600L,
  spec$crossfit_folds == 10L, spec$oracle_landmarks == 5L,
  spec$configuration$tau == log(8), spec$configuration$transform$grid_size == 12289L,
  nrow(raw) == spec$B * length(spec$true_HR) * 32L,
  all(raw$n1 == spec$n1), all(raw$n0 == spec$n0), all(raw$total_n == spec$total_n),
  all(with(raw, table(true_HR, PS, OR, method)) == spec$B),
  !anyDuplicated(raw[, c("true_HR", "replication", "PS", "OR", "method")]))
for (method in c("Naive", "Oracle")) {
  z <- raw[raw$method == method, , drop = FALSE]
  stopifnot(all(vapply(split(z$theta, interaction(z$true_HR, z$replication)),
    function(v) length(unique(v)) == 1L, logical(1))))
}
ipw <- subset(raw, method == "IPW")
stopifnot(all(vapply(split(ipw$theta,
  interaction(ipw$true_HR, ipw$replication, ipw$PS)),
  function(v) length(unique(v)) == 1L, logical(1))))
valid <- raw$converged & is.finite(raw$theta) & is.finite(raw$HR) & raw$HR > 0
stopifnot(all(!raw$converged | valid),
  max(abs(exp(raw$theta[valid]) - raw$HR[valid])) < 1e-12)
cal <- subset(raw, method %in% c("MEC1", "MEC2", "MEC5", "MEC10", "MEC20", "Oracle"))
cal_valid <- cal$converged & is.finite(cal$HR)
stopifnot(all(cal$calibration_converged[cal_valid]),
          all(cal$original_mean_imbalance[cal_valid] <= 1e-10))
full_svd <- cal_valid & cal$calibration_coordinates == "Full invertible SVD"
if (any(full_svd)) stopifnot(all(cal$calibration_kkt_residual[full_svd] <= 1e-6))
paired_cal <- merge(subset(cal, true_HR == 1), subset(cal, true_HR == 0.70),
  by = c("replication", "PS", "OR", "method"), suffixes = c("_null", "_alternative"))
stopifnot(identical(paired_cal$calibration_converged_null,
                    paired_cal$calibration_converged_alternative),
  identical(paired_cal$original_mean_imbalance_null,
            paired_cal$original_mean_imbalance_alternative))

# Check that the HR-paired data generator uses the same covariates and controls.
seed <- simulation_seed(1L, spec$master_seed)
set.seed(seed)
null <- dgp_draw_data(spec$n0, spec$n1, dgp_config(hr = 1))
set.seed(seed)
alternative <- dgp_draw_data(spec$n0, spec$n1, dgp_config(hr = 0.70))
stopifnot(identical(null[, paste0("X", 1:10)], alternative[, paste0("X", 1:10)]),
  identical(null[null$A == 0, ], alternative[alternative$A == 0, ]),
  all(alternative$time[alternative$A == 1] >= null$time[null$A == 1] - 1e-10))

counts <- aggregate(as.integer(valid), raw[c("true_HR", "PS", "OR", "method")], sum)
names(counts)[ncol(counts)] <- "valid_fits"
counts$generated_replications <- spec$B
counts$failed_fits <- spec$B - counts$valid_fits
write.csv(counts, "results/fit_counts.csv", row.names = FALSE)
write.csv(subset(raw, !valid), "results/failed_fits.csv", row.names = FALSE)
diagnostics <- list(integrity_passed = TRUE, generated_datasets = spec$B * length(spec$true_HR),
  total_rows = nrow(raw), valid_rows = sum(valid), failed_rows = sum(!valid),
  failure_reasons = table(raw$reason[!valid]),
  max_original_mean_imbalance = max(cal$original_mean_imbalance[cal_valid]),
  calibration_coordinates = table(cal$calibration_coordinates, cal$converged),
  paired_draws_verified = TRUE, source_hashes_verified = TRUE)
saveRDS(diagnostics, "results/integrity_checks.rds")
print(diagnostics)
print(counts[counts$failed_fits > 0, ])

