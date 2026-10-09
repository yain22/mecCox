# mecCox 0.1.0

- Added unweighted and ATT-IPW marginal Cox comparisons. The ATT-IPW fit
  reports naive, robust, and corrected standard errors for one coefficient.
- Added KL-calibrated MEC-Cox fits with cross-fitted source propensity and
  control-survival predictions, plus working stacked-sandwich inference.
- Added baseline covariate and predicted mortality-risk balance diagnostics.
- Added an entirely simulated external-control example and regression tests.
- Added display-only scripts for both simulation experiments, with editable
  replication counts and parallel execution using up to 20 workers by default.
  Complete summary, reference-target, and configuration tables appear as
  `kableExtra` HTML tables in the RStudio Viewer; figures use the active
  graphics device, and replication-level results remain in memory.
