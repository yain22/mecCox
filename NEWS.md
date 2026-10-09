# mecCox 0.1.2

- Provide the Figure 4/5 study configurations and fitting functions with
  deterministic random streams, fixed reference targets, and successful-fit
  summaries based on the study confidence limits.
- Report completed datasets and partial results for both simulation scripts,
  with implementation identifiers in run metadata.
- Quick checks reduce the number of datasets and retain study learner settings.

# mecCox 0.1.1

- Added a small-subset BART tuning and randomized survival-forest preset,
  selected with `nuisance_settings = "original_study"` in the fitting API.
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
