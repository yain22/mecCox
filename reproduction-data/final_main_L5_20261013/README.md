# Final main simulation data

This compact bundle contains the completed data for Figures 5 and 6: 30,000 simulation datasets, with 1,000 attempted replications in every design cell. The RDS retains all 135,000 model rows, including failures and finite extreme results.

From the repository root, redraw both figures without fitting any models:

```sh
Rscript inst/reproduce/redraw_final_figures.R --reference=reproduction-data/final_main_L5_20261013/final_figure_reference.rds --output=final-main-figures
```

Use a fresh R session with ggplot2 4.0.0, patchwork 1.3.2 and gtable 0.3.6; dplyr, labeling and digest are also required. Put that rendering library on `R_LIBS_USER` before starting R. The output directory must be new or empty. The command writes both 12-by-9-inch PNG/PDF figures and their plotted data. The fitting dependency snapshot is separate and remains in `inst/reproduce/final_reference`.

Files:

- `final_figure_reference.rds`: all attempted model rows, primary summaries, numerical targets and portable provenance.
- `primary_summary.csv`: the 135 plotted summary cells with explicit denominators.
- `all_target_and_eligibility_summaries.csv`: primary, legacy-target and numerical-target uncertainty summaries, with both original finite-fit and stricter diagnostic eligibility rules.
- `numerical_targets.csv`: the numerical targets and Monte Carlo uncertainty.
- `raw_calibration_status.csv`: calibration convergence and error counts.
- `SHA256.csv`: integrity hashes for the data files and this README.

The plotted policy is `independent_primary` with `legacy_finite` eligibility: finite estimates and standard errors for bias/RMSE, and finite ordered intervals additionally required for coverage. Failures are retained in the attempt counts. The full raw-row CSV is omitted because the same rows are already in the RDS.

The completed private and public summaries agreed exactly (absolute tolerance zero). Both PNGs and the rendered PDF pages agreed pixel for pixel in the verified rendering environment. PDF file metadata may differ; device and font details are documented in `inst/reproduce/final_plot_reference/README.txt`.

This directory is available in the repository clone and excluded from the installed R package. To rerun simulations rather than redraw the saved results, use the current `scenario1.R` and `scenario2.R` drivers under `inst/reproduce`. See `inst/reproduce/FINAL_FIGURE_WORKFLOW.txt` for the complete workflow.
