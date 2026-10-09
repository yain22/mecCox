# Run the simulations and public breast-cancer case study

This directory contains runnable code for three analyses:

| Analysis | Script | Output |
| --- | --- | --- |
| Simulation 1 | [scenario1.R](scenario1.R) | Simulation summary and figure |
| Simulation 2 | [scenario2.R](scenario2.R) | Simulation summary and figure |
| Public breast-cancer case study | [breast_cancer.R](breast_cancer.R) | Two analysis tables in the format of supplementary Tables S2 and S3 |

The breast-cancer data are public datasets included in `survival`; they are
not subject to the restrictions of the separate SQUIRE/MSK-CHORD application.
No restricted patient records are included in these examples. All three
scripts display results and retain R objects without exporting result files.
The [case-study instructions](#public-breast-cancer-case-study) appear below
the simulation instructions.

## Simulation dependencies

`scenario1.R` and `scenario2.R` rerun the two main-paper simulation designs
with the exported `mecCox` fitting functions. Install the package first, for
example with `R CMD INSTALL mecCox` from the checkout's parent directory.
Both scripts use `kableExtra`, `htmltools`, and `rstudioapi` to present HTML
tables. Scenario 2 also requires the optional learners `dbarts` and `ranger`:

```r
install.packages(c("kableExtra", "htmltools", "rstudioapi"))
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
in R paths on Windows. The helper is found from the sourced file's location,
the current project, or the installed package. If using downloaded copies,
keep `simulation_helpers.R` beside both scenario scripts.

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
across design cells. Only Monte Carlo runs execute in parallel; reference
target computation, summaries, HTML tables, and plots run in the main R
process. Each run receives its own deterministic seed, so changing
worker count or scheduling preserves its data and fits in the same R and
package environment.

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
the conditional coefficient `log(0.70)`. With the stated seed, the full-size
reference target is approximately `-0.2237787` on the log-hazard-ratio scale;
the short `--quick` run deliberately uses a smaller reference population.

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
runs unless a different count is selected, and a reference population
of 2,000 treated and 4,000 external controls. It retains all 50 covariates,
ten folds, and five landmarks.

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
settings follow Scenario 1. Reference targets are computed separately for
each nonlinearity setting using 30,000 treated and 60,000 external controls.
The three logistic ATT-IPW comparators still share one point estimate. The
two MEC-Cox variants use cross-fitted BART source
probabilities and either Cox or RSF control-survival predictions, with ten
folds and five landmarks. The full script uses 100 BART trees, 1,000 posterior
draws, 500 burn-in iterations, and shrinkage parameter 2. The RSF uses 500
trees, with `mtry` and minimum node size tuned within each training fold;
the fallback minimum node size is 15.
The two MEC-Cox variants use the same BART settings, folds, and seeds, so
their source-propensity predictions agree. The script and metadata record
the learner settings.
The study script fixes these BART settings, while the original study
code tuned BART within training folds. It runs the study design with
the public API; its fitting sequence and numerical results can differ from
those underlying the paper's figure.

The short `--quick` run keeps all three nonlinearity settings, but uses
`n1 = 200`, two runs per setting unless a different count is selected,
and reference populations of 2,000 treated and 4,000 external controls.
It uses 25 BART trees, 50 posterior
draws, 25 burn-in iterations, and 100 RSF trees without tuning. Ten-fold
cross-fitting and five landmarks are retained. It exercises both MEC-Cox
variants. Quick-run results cannot establish simulation performance.

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

The deterministic seeds make each script rerunnable, but the package API,
learner tuning, and random-number streams are not a bit-for-bit replay of the
earlier private parallel code used to produce the paper's figures. The
package's corrected ATT-IPW sandwich also differentiates the clipped analysis
weights locally; its values can differ from an earlier correction when fitted
probabilities cross a clipping boundary. Results should be compared at the
Monte Carlo level; a quick run cannot establish performance.
