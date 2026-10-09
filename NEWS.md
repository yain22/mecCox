# mecCox 0.1.1

- Restored Scenario 2's original small-subset BART tuning and randomized
  survival-forest settings. The general fitting API retains its existing
  defaults; `nuisance_settings = "original_study"` selects the study preset.
- Reuse cross-fitted propensity predictions between the two Scenario 2 fits,
  and omit unused learner predictions and derivative influence calculations.
- Report Scenario 2 progress within cells and retain completed datasets in
  memory after interruption, with explicit partial-result status.
- Resolve breast-cancer helper files from saved scripts, RStudio documents,
  project folders, or the installed package, with an explicit directory option.

# mecCox 0.1.0

- Added a public GBSG/Rotterdam breast-cancer case-study script that generates
  the two supplementary tables in the RStudio Viewer without result exports.
- Added references for Shu et al., Binder, and Lin--Wei, and a DAG hex logo.
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
