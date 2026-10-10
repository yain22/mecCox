# Run the toy example, simulations, and public breast-cancer case study

This directory contains four runnable analysis scripts:

| Analysis | Script | Output |
| --- | --- | --- |
| Toy example | [toy_example.R](toy_example.R) | Three-panel figures and numerical summaries for both designs |
| Simulation 1 | [scenario1.R](scenario1.R) | Simulation summary and figure |
| Simulation 2 | [scenario2.R](scenario2.R) | Simulation summary and figure |
| Public breast-cancer case study | [breast_cancer.R](breast_cancer.R) | Two analysis tables in the format of supplementary Tables S2 and S3 |

The breast-cancer data are public datasets included in `survival`; they are
not subject to the restrictions of the separate SQUIRE/MSK-CHORD application.
No restricted patient records are included in these examples. The scripts
display results and retain R objects. The toy script also supports optional
file export; the other three scripts do not export result files.
The [case-study instructions](#public-breast-cancer-case-study) appear below
the simulation instructions.

## Toy example

The toy script reproduces both prognostic-balance designs in the paper.
It can run directly from a checkout without installing `mecCox`. Install
its dependencies once:

```r
install.packages(c("survival", "ggplot2", "patchwork"))
```

Each design uses 10,000 simulated datasets, 200 treated patients, 400 external
controls, and seed `20261007`. The binary covariate `X1` has probability 0.7
in the treated cohort and 0.3 in the external-control cohort. `X2` is an
independent standard normal variable in both cohorts. Event times under
treatment and control have the same exponential hazard, so the true ATT
log hazard ratio is zero:

| Design | Event hazard under either treatment | Oracle score at time 5 |
| --- | --- | --- |
| Precision gain | `0.08 * exp(0.5 * X1 + X2)` | `exp(-0.4 * exp(0.5 * X1 + X2))` |
| No additional precision gain | `0.08 * exp(0.5 * X1)` | `exp(-0.4 * exp(0.5 * X1))` |

Censoring is independently exponential with rate 0.03, with administrative
censoring at time 5. Both designs fit ATT odds from the logistic source model
with an intercept and `X1`, and use KL calibration with an intercept, `X1`,
and the oracle prognostic score. The two weighted Cox regressions are fitted
separately in every dataset. The second design still generates `X2`, keeping
the original random-draw order; only its coefficient in the hazard changes.
Its oracle score is a function of binary `X1`, which is already balanced by
the baseline weights. The calibration constraint is therefore redundant.

From the repository root, run both designs with one command:

```sh
Rscript inst/reproduce/toy_example.R
```

To check the code quickly or choose a replication count:

```sh
Rscript inst/reproduce/toy_example.R --quick
Rscript inst/reproduce/toy_example.R --replications=100
```

`--quick` uses 20 datasets per design when the configured replication count
is unchanged. An explicit `--replications=N` sets the count per design and
must be an integer of at least two. Quick runs check execution; they do not
reproduce the paper's empirical SDs. The script runs the precision-gain
design first, then the second design, resetting the seed before each.
Serial simulation preserves the original sequence of random draws.

In RStudio, open the script and click **Source**, or use `source()`:

```r
source("inst/reproduce/toy_example.R")
```

After installing the current GitHub package, the script is also available
through `system.file()`:

```r
source(system.file("reproduce", "toy_example.R", package = "mecCox"))
```

The script has editable `quick_run`, `replications`, `seed`, and `output_dir`
settings near the top. Set these inside the script before sourcing it.
Absolute script paths work from other working directories. When downloading
the files separately, keep `toy_helpers.R` in the same directory.

The results remain in `toy_example_results`, with `precision_gain` and
`no_precision_gain` components. Each component contains `replications`,
`summary`, `diagnostics`, `figure`, and `metadata`. For example:

```r
toy_example_results$precision_gain$summary
toy_example_results$no_precision_gain$summary
print(toy_example_results$precision_gain$figure)
print(toy_example_results$no_precision_gain$figure)
```

The script prints both numerical summaries. In an interactive graphics
session, the figures appear in sequence in the active graphics device,
normally RStudio's **Plots** pane. It creates no result files by default.
For explicit export:

```sh
Rscript inst/reproduce/toy_example.R --output-dir=toy-example-output
```

The output directory receives both three-panel figures as PDF and PNG,
replication-level estimates, summary and diagnostic CSVs, and R session
information. The filenames begin with `toy_precision_gain` or
`toy_no_precision_gain` to distinguish the designs. Exporting is optional
and does not change the simulation.

Expected results from the full manuscript runs are:

| Design | ATT-IPW Cox empirical SD | MEC-Cox empirical SD | MEC-Cox / ATT-IPW empirical variance |
| --- | ---: | ---: | ---: |
| Precision gain | 0.145398 | 0.123542 | 0.721960 |
| No additional precision gain | 0.152738 | 0.152738 | 1.000000 |

The figures use the same single-line gap labels as the manuscript. Panel
(a) compares empirical prognostic-score gaps; panel (b) plots the change in
the Cox estimate against the initial ATT-IPW gap; panel (c) compares the
sampling distributions. In the second design, panels (a) and (b) show a point
mass at zero and coincident points. Residuals below `1e-10` are displayed as
zero there; raw replication values are retained. Fixed seeds reproduce the
numerical results in the same R and package environment; exported session
information records that environment.

## Simulation dependencies

`scenario1.R` and `scenario2.R` run the two main-paper simulation designs
with the study fitting functions and fixed configurations. Install the package first, for
example with `R CMD INSTALL mecCox` from the checkout's parent directory.
Both scripts use `kableExtra`, `htmltools`, and `rstudioapi` to present HTML
tables, and `rngtools` for deterministic random-stream allocation.
Scenario 2 also requires the optional learners `dbarts` and `ranger`:

```r
install.packages(c("kableExtra", "htmltools", "rstudioapi", "rngtools"))
# Also install these learners for Scenario 2:
install.packages(c("dbarts", "ranger"))
```

## Run simulations from R or RStudio

Both simulation scripts support `source()` and RStudio's **Source** button. After
installing the package, run Scenario 1 with:

```r
script <- system.file("reproduce", "scenario1.R", package = "mecCox")
source(script)
```

For Scenario 2, use `"scenario2.R"` in `system.file()` instead. You can also
source a script from your local checkout:

```r
source("path/to/mecCox/inst/reproduce/scenario1.R")
```

The configuration block near the top of each script contains ordinary R
settings:

```r
quick_run <- FALSE
cores <- 20L
replications <- 1000L
```

Edit `replications` to choose the number of Monte Carlo runs per
design cell. The default is 1,000, matching the full simulation design.
Edit `cores` to choose the maximum worker count.

For Scenario 2, `replications <- 100L` with `quick_run <- FALSE` means
100 runs in each of 15 design cells: **1,500 simulated datasets**, each with
two MEC-Cox fits. With `quick_run <- TRUE`, the same explicit count gives
100 runs in each of three cells, or 300 datasets. For a small installation
check, set all three values as in the example below.

**The scripts keep results in memory and display them; they do not export
result files.** In RStudio, the complete summary, reference targets, and run
settings appear as formatted `kableExtra` HTML tables in the **Viewer** pane.
The figure is drawn on the active graphics device, normally the **Plots**
pane. The underlying tables, HTML objects, and metadata remain available in
the session; see [Results in the R session](#results-in-the-r-session) below.

To check the installation first, open the script with `file.edit(script)` and
edit its configuration block to:

```r
quick_run <- TRUE
cores <- 2L
replications <- 2L
```

Run the edited whole file with **Source** or `source(script)`.
For Scenario 2, edit the same settings inside `scenario2.R` and source that
file instead. Edit the settings **inside the script** before sourcing it:
values assigned only in the console are replaced by its configuration block.

An absolute script path works from any working directory; use forward slashes
in R paths on Windows. Helpers are found from the sourced file's location
or the installed package. When downloading files, keep the complete
`inst/reproduce` directory together, including `simulation_helpers.R`,
`original_study_helpers.R`, and the `study_reference` subdirectory.

## Run simulations from a terminal

From the checkout's parent directory:

```sh
Rscript mecCox/inst/reproduce/scenario1.R
Rscript mecCox/inst/reproduce/scenario2.R
```

For short installation and code-path checks:

```sh
Rscript mecCox/inst/reproduce/scenario1.R --quick --cores=2
Rscript mecCox/inst/reproduce/scenario2.R --quick --cores=2
```

To choose the number of runs explicitly:

```sh
Rscript mecCox/inst/reproduce/scenario1.R --cores=20 --replications=100
Rscript mecCox/inst/reproduce/scenario2.R --cores=20 --replications=100
```

`--replications=N` accepts a positive integer and overrides the setting in
the script. With `--quick`, the count defaults to two if the configured count
is still 1,000 and no command-line count is supplied. A nondefault configured
count or an explicit command-line count is retained.

Without an interactive graphics window, the scripts print their tables to
the console; use RStudio for the HTML report and plots. Command-line options
override the corresponding script settings. When running with
`source()` or interactively, the script uses its configuration block rather
than unrelated R session command-line arguments.

## Parallel simulation execution

The scripts request 20 workers by default and cap the actual worker count
at the number of detected logical cores and runs per design cell.
Choose another limit with `--cores=N`, or use `--cores=1` for serial execution
from a terminal. In R or RStudio, edit `cores` in the script instead:

```sh
Rscript mecCox/inst/reproduce/scenario1.R --cores=8
```

The PSOCK backend works on Windows, macOS, and Linux. A worker pool is reused
across design cells. Monte Carlo runs execute in parallel; summaries, HTML
tables, and plots run in the main R process. Reference targets are read from
the study configurations. Each dataset receives the L'Ecuyer random
stream allocated by the sample-size-first job grid.
Changing worker count or scheduling therefore preserves its data and fits
in the same R and package environment.

## Study implementation

The `study_reference` directory contains the fitting functions and
configurations for Figures 4 and 5. The scripts use fixed cross-fitting and
forest seeds, study-specific tuning, and a `1.96` confidence-limit multiplier.
Summary metrics include finite estimates with finite standard errors. Thus
1,000 runs means 1,000 attempted datasets per cell; successful counts are
reported separately. Implementation identifiers are available in the result
metadata.

**BART score extraction.** The study functions average latent predictions
and use a range-based probability check; otherwise they apply
`pnorm(raw + binary_offset)`. In `dbarts` 0.9-32, `yhat.test` already
includes `binaryOffset`, so that branch adds the offset twice. This is a
known limitation of the study calculation, not a general-purpose posterior
probability estimator. The documented posterior probability mean is
`colMeans(pnorm(fit$yhat.test))`.

The study environment used R 4.5.1, `dbarts` 0.9-32, `ranger` 0.17.0,
and `survival` 3.8-3. Numerical results may vary with software versions.

## Scenario 1: linear source selection and prognosis

At the default setting, the full study has 1,000 Monte Carlo runs at
each of 15 combinations of
`n1 = 200, 250, 300, 350, 400` and `n0/n1 = 2, 3, 4`. The source and outcome
models have 50 independent standard-normal covariates. The first five affect
source membership, the first ten affect the event hazard, and the remaining
40 are outcome noise. The source logit intercept is `-0.2`; its five slopes
are `0.75, 0.75, 0.65, 0.65, 0.55`. Source probabilities are clipped to
`[0.02, 0.98]` during generation. Cohort sizes are fixed by sampling from
the resulting source-conditional covariate distributions, without changing
the logit intercept. Consequently this is a linear-logit design with overlap
clipping, rather than an exactly logistic propensity throughout the tails.

The conditional event hazard is Weibull with scale `0.00008`, shape `2`, and
conditional log-hazard ratio `log(0.70)`. The ten prognostic coefficients are
the logarithms of `1.75, 1.75, 1.60, 1.60, 1.50` and five copies of `1.25`.
Independent exponential censoring has rate `0.0008`. The common reference
target is the ATT-weighted marginal Cox projection, approximated using 30,000
treated and 60,000 external-control superpopulation observations, **not**
the conditional coefficient `log(0.70)`. The saved reference target is
approximately `-0.2237787` on the log-hazard-ratio scale. The scripts use that
saved target in both full and quick runs.

The three ATT-IPW rows share the same weighted Cox coefficient and differ
only in their naive, robust sandwich, or corrected sandwich standard error.
MEC-Cox uses a ten-fold cross-fitted logistic source model and control-only
Cox survival learner. The five landmark times are the type-8 empirical
quantiles at probabilities `0.10, 0.30, 0.50, 0.70, 0.90` of observed
external-control event times. The survival predictions, along with an
intercept, form the KL calibration basis. Estimated source probabilities are
clipped to `[0.01, 0.99]` when constructing analysis weights. Both methods
use a Breslow weighted Cox fit.

The short `--quick` run retains `n1 = 200` and `n0/n1 = 2`, with two
runs unless a different count is selected. It retains the study's target,
all 50 covariates, ten folds, and five landmarks.

## Scenario 2: increasing nonlinearity

The full run uses 10 independent standard-normal covariates, the same five
treated sample sizes, and `n0/n1 = 4`. At the default setting, it
has 1,000 Monte Carlo runs at each of 15 sample-size and nonlinearity combinations:

| Setting | Source-selection multiplier `kappa_pi` | Prognostic multiplier `kappa_m` |
| --- | ---: | ---: |
| None | 0 | 0 |
| Mild | 1 | 2 |
| Severe | 2 | 5 |

The source logit adds `kappa_pi * r_pi(X)` to Scenario 1's linear component;
the control log-hazard adds `kappa_m * r_m(X)` to its linear component. The
nonlinear terms are:

```r
r_pi <- 0.70 * sin(1.25 * X1) + 0.45 * (X2^2 - 1) -
  0.55 * (as.numeric(X3 > 0) - 0.5) + 0.35 * X4 * X5 +
  0.25 * (cos(X1 + X2) - exp(-1))
r_m <- 0.45 * sin(X2) + 0.35 * (X3^2 - 1) +
  0.30 * (as.numeric(X4 > 0) - 0.5) + 0.25 * X1 * X5 +
  0.20 * (cos(X2 + X5) - exp(-1))
```

The remaining event-time, censoring, probability-clipping, and landmark
settings follow Scenario 1. Saved reference targets were computed separately
for each nonlinearity setting using 30,000 treated and 60,000 external controls.
The three logistic ATT-IPW comparators still share one point estimate. The
two MEC-Cox variants use cross-fitted BART source
probabilities and either Cox or RSF control-survival predictions, with ten
folds and five landmarks. The script uses lightweight tuning. Within each
training fold, BART selects 25, 50, or 100 trees using a source-balanced tuning
subset of at most 100 patients. Candidate and selected fits use 100 posterior
draws, 50 burn-in iterations, and shrinkage parameter 2. If tuning cannot be
performed, the fallback is 50 trees, 200 posterior draws, and 100 burn-in
iterations.

RSF tuning uses at most 100 external controls, split approximately equally
into training and validation sets. It selects `mtry` and minimum node size
from a nine-candidate grid, with 100 trees. Forests use randomized
splits (`extratrees`), one candidate split per variable, and a 63.2% sample
without replacement. The fallback is 300 trees and minimum node size 15.
The two MEC-Cox variants use fixed folds and tuning seeds. RSF tuning uses
orientation-free `max(C, 1 - C)` concordance criterion. Final forest fits use
the fixed study seed. The script and metadata record the
implementation and settings. More workers also require more memory; reduce
`cores` if concurrent fits exhaust available RAM.

Scenario 2 reports completed datasets while it runs. If interrupted, it
returns the results already received from workers in `scenario2_results`,
with `status = "interrupted"`, completion counts, and a partial summary.
Unfinished calculations are excluded. Results remain in memory; save any
objects you need before restarting R. Worker shutdown does not forcibly
terminate a native learner that is still computing.

The short `--quick` run keeps all three nonlinearity settings, but uses
`n1 = 200` and two runs per setting unless a different count is selected.
It retains BART/RSF tuning, ten-fold cross-fitting, five
landmarks, and saved reference targets. Both MEC-Cox variants are fitted.
Quick-run results cannot establish simulation performance.

## Public breast-cancer case study

`breast_cancer.R` uses the public `gbsg` and `rotterdam` datasets shipped
with `survival`. It loads them directly from the installed package; no
private data file or access request is needed. The target cohort comprises
GBSG patients who received hormonal therapy, and the external-control
cohort comprises Rotterdam patients who did not receive hormonal therapy.
The endpoint is recurrence-free survival administratively censored at five
years. Rotterdam recurrence and death records are combined to harmonize
the endpoint with the GBSG data.

The analysis uses age, menopausal status, tumor-size category, tumor grade,
and `log(1 + x)` transformations of positive-node count, progesterone
receptor level, and estrogen receptor level. It compares unweighted Cox,
ATT-IPW Cox with three standard errors for the same point estimate, and
MEC-Cox with GLM/Cox and DL/RSF learners. The MEC-Cox fits use ten-fold
cross-fitting and 20 survival landmarks.

Install the display packages and the optional DL/RSF dependencies:

```r
install.packages(c("kableExtra", "htmltools", "rstudioapi", "brulee", "torch", "ranger"))
# If the Torch runtime is not already installed:
torch::install_torch()
```

Run the script with RStudio's **Source** button or from the console:

```r
script <- system.file("reproduce", "breast_cancer.R", package = "mecCox")
source(script)
```

For a local checkout, use
`source("path/to/mecCox/inst/reproduce/breast_cancer.R")`. Keep
`breast_cancer_helpers.R` and `simulation_helpers.R` beside the script when
downloading individual files.
The script also checks the saved RStudio document location and project
ancestors. When using pasted code, an unsaved editor document has no directory;
install the current package or set `reproduce_dir` at the top of the script
to the folder containing both helpers. You can also set the persistent
option `options(mecCox.reproduce_dir = "path/to/mecCox/inst/reproduce")`.
An older installed package may lack the case-study helpers. Update it once:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("yain22/mecCox", upgrade = "never", force = TRUE)
```

Restart R after installation, then rerun the script. Helpers are loaded from
one directory so downloaded and installed versions are not mixed.

The two tables open together as `kableExtra` HTML tables in the RStudio
**Viewer**. They present the supplementary breast-cancer analyses:

- **Table S2 format:** log-hazard-ratio estimates, standard errors, hazard
  ratios, and 95% confidence intervals for all six method/variance combinations.
- **Table S3 format:** covariate differences before and after weighting,
  mean and maximum absolute balance differences, counts exceeding 0.10,
  effective sample size, and the coefficient of variation of external-control weights.

Continuous covariates are standardized using the target-cohort standard
deviation. Binary covariates use differences in proportions, matching the
supplement's reporting convention. The unweighted estimate is a
descriptive comparison; interpretation of weighted estimates depends on the
study's identifying assumptions.

The sidecar `breast_cancer_helpers.R` preserves the original case-study
model conventions. Grade 1 occurs in 33 treated patients and no external
controls, so weighting cannot remove that category's imbalance. The control
Cox fit retains the original factor coding and permits an aliased coefficient
(set to zero for prediction); prognosis for the unsupported grade therefore
depends on an extrapolation convention. The RSF uses the same factor levels.
The general `fit_mec_cox()` interface keeps its stricter support checks.
The Rotterdam endpoint also retains the original convention of counting the
first recorded recurrence or death and using the later follow-up time when
neither occurred; see the [dataset documentation](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/rotterdam.html)
for the issue of deaths recorded after recurrence follow-up ended.

The script creates no result exports. Its objects remain available in R:

```r
breast_cancer_estimates   # Table S2 format: estimates and uncertainty
breast_cancer_balance     # Table S3 format: balance and weight diagnostics
breast_cancer_results     # analysis results and metadata
breast_cancer_tables      # formatted HTML tables
breast_cancer_report      # report displayed in the Viewer
```

From a terminal, use `Rscript mecCox/inst/reproduce/breast_cancer.R`.
Tables print to the console when no RStudio Viewer is available. Numerical
results may vary from the manuscript tables with random seeds, learner
settings, and software versions.

## Results in the R session

When sourced, Scenario 1 creates:

- `scenario1_replications`: one row per run and method, with a visible
  error message for any failed fit.
- `scenario1_summary`: nominal 95% Wald coverage, Monte Carlo bias and RMSE on
  the **log-hazard-ratio** scale, and counts of successful and failed fits.
- `scenario1_results`: a list containing the table of individual runs, summary table,
  and run metadata, including the design, reference target, quick/full flag,
  worker counts, and R session.
- `scenario1_tables`: a named list of formatted HTML summary, reference-target,
  and configuration tables.
- `scenario1_report`: the HTML report displayed in the RStudio Viewer.

Scenario 2 creates the corresponding `scenario2_replications`,
`scenario2_summary`, `scenario2_results`, `scenario2_tables`, and
`scenario2_report`, plus `scenario2_targets`, which contains the reference
log-hazard ratio for each nonlinearity setting.

The HTML report includes all design cells, methods, performance metrics, and
failure counts in the aggregated summary. Individual simulation rows remain
in memory rather than being added to a very large Viewer report. You can
inspect the underlying objects without rerunning the simulation:

```r
scenario1_tables$summary
scenario1_summary
head(scenario1_replications)
scenario1_results$metadata
```

To redraw the figure without refitting, use `plot_scenario1_results(scenario1_summary)`.
For Scenario 2, use:

```r
plot_scenario2_results(scenario2_summary, scenario2_results$metadata$settings)
```

Before interpreting a full run, check the summary's failure counts and inspect
their messages in the table of individual runs. These scripts do not infer or replace
missing results when a fit fails.

## Numerical results

The study scripts specify the simulation designs, fitting procedures, and
random streams. Numerical results may vary with software versions. Use the
full design to assess Monte Carlo performance; quick runs check installation
and execution only.
