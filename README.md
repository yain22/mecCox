<img src="man/figures/logo.svg" align="right" alt="mecCox hex logo" width="120" />

# mecCox

[![R package](https://img.shields.io/badge/R-%3E%3D%204.2.0-276DC3?logo=r&logoColor=white)](DESCRIPTION)
[![R package check](https://github.com/yain22/mecCox/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/yain22/mecCox/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)

`mecCox` implements the Cox comparisons and prognostic weight calibration used
in **Balancing Machine-Learned Prognostic Scores to Improve Efficiency in ATT
Marginal Hazard-Ratio Estimation with Inverse-Probability-Weighted Cox
Regression**. It is a research package prepared for GitHub; it is not a CRAN
submission.

The target cohort is coded `source = 1` and retains weight one. External
controls are coded `source = 0`. The package includes an unweighted Cox
reference, a logistic-propensity ATT-IPW Cox fit, and MEC-Cox fits that calibrate
external-control odds weights to cross-fitted control-survival predictions.

## Installation

Install the package from GitHub:

```r
install.packages("remotes")
remotes::install_github("yain22/mecCox", dependencies = NA)
```

The Cox/GLM implementation needs only `survival`. An MLP propensity learner
also needs `brulee` and its `torch` runtime; a random survival forest needs
`ranger`, and a BART propensity learner needs `dbarts`. These learners are
optional so the core package remains usable without them. In a noninteractive
R installation, `torch` may need a separate `torch::install_torch()` call to
install its runtime.

## First analysis

```r
library(mecCox)
data(example_external_controls)

covariates <- c("age", "sex", "marker")

unweighted <- fit_unweighted_cox(
  example_external_controls,
  time = "time", event = "event", source = "source"
)

ipw <- fit_att_ipw_cox(
  example_external_controls,
  time = "time", event = "event", source = "source",
  covariates = covariates
)

mec <- fit_mec_cox(
  example_external_controls,
  time = "time", event = "event", source = "source",
  covariates = covariates,
  ps_learner = "glm", survival_learner = "cox",
  n_folds = 5, n_landmarks = 5, seed = 2026
)

unweighted$hr
ipw$se                 # naive, robust, corrected: one coefficient
mec$hr
mec$calibration$max_mean_balance_residual

covariate_balance(mec)
predicted_risk_balance(mec, times = c(3, 6, 9, 12, 15))
```

The included `example_external_controls` data are entirely simulated. They
contain no SQUIRE or MSK-CHORD patient records.

## Run the simulations and public case study

The repository provides runnable code for two simulation experiments and
one public-data case study:

| Study | Script | Results |
| --- | --- | --- |
| Simulation 1: linear source selection and prognosis | [scenario1.R](inst/reproduce/scenario1.R) | Simulation summary and nine-panel figure |
| Simulation 2: increasing nonlinearity | [scenario2.R](inst/reproduce/scenario2.R) | Simulation summary and nine-panel figure |
| Breast-cancer case study | [breast_cancer.R](inst/reproduce/breast_cancer.R) | Hazard-ratio estimates and covariate/weight diagnostics in the format of supplementary Tables S2 and S3 |

The breast-cancer example uses the public `gbsg` and `rotterdam` datasets
distributed with the R package `survival`. These data are available to anyone
who installs that package. The separate SQUIRE/MSK-CHORD application uses
restricted patient data, which are not distributed here.

### Simulation experiments

The [Scenario 1 script](inst/reproduce/scenario1.R) runs the paper's
linear source-selection and outcome experiment. It compares the three
ATT-IPW Cox standard errors with GLM/Cox MEC-Cox for each treated sample size
`n1 = 200, 250, 300, 350, 400` and external-control ratio `n1:n0 = 1:2, 1:3,
1:4`. It uses 50 covariates, five event-time-quantile survival landmarks,
ten-fold cross-fitting, and 1,000 Monte Carlo runs per sample-size
combination. The results include estimates from individual runs, a table of
coverage, bias, and RMSE, and a nine-panel figure.

The [Scenario 2 script](inst/reproduce/scenario2.R) evaluates the three
nonlinearity settings with 10 covariates, a fixed external-control ratio of
`1:4`, and the same treated sample sizes and number of runs. It compares
the ATT-IPW methods with BART/Cox and BART/RSF MEC-Cox. Install its optional
learners before running this script:

```r
install.packages(c("dbarts", "ranger"))
```

Both scripts use HTML tables from `kableExtra` for their results display:

```r
install.packages(c("kableExtra", "htmltools", "rstudioapi"))
```

Run either script with RStudio's **Source** button, or from the R console:

```r
script <- system.file("reproduce", "scenario1.R", package = "mecCox")
source(script)
```

Use `"scenario2.R"` in the same command for Scenario 2. Each script has a
configuration block near the top:

```r
quick_run <- FALSE
cores <- 20L
replications <- 1000L
```

Edit `replications` to choose the number of Monte Carlo runs per
design cell, and edit `cores` to choose the maximum number of workers.
**The scripts display results and keep them in memory; they do not export
result files.** In RStudio, a formatted HTML report opens in the **Viewer**
pane, with the complete simulation summary, reference targets, and run
settings. The figure appears on the active graphics device, normally the
**Plots** pane. The results also remain in the R session:

```r
scenario1_summary       # coverage, bias, RMSE, and fit counts
scenario1_replications  # individual estimates and standard errors
scenario1_results       # individual runs, summary, and run metadata
scenario1_tables        # named list of formatted HTML tables
scenario1_report        # HTML report shown in the Viewer
```

Scenario 2 creates `scenario2_summary`, `scenario2_replications`,
`scenario2_results`, `scenario2_targets`, `scenario2_tables`, and
`scenario2_report`. Individual simulation rows remain available in R; the
Viewer report presents the aggregated results rather than thousands of raw
rows. To run a short check first, open the script with `file.edit(script)`,
set `quick_run <- TRUE`, `cores <- 2L`, and `replications <- 2L` in its
configuration block, then source the edited script. When using downloaded
copies, keep `simulation_helpers.R` in the same folder.

At the default 1,000 runs per design cell, each full study fits 15,000 simulated
datasets and can take substantial time.
Simulation runs use up to 20 parallel workers by default, capped by the
number of detected logical cores and runs. Edit `cores` to choose
a smaller worker count, or set `cores <- 1L` for a serial run. The socket-based
parallel backend works on Windows, macOS, and Linux. Fixed simulation seeds
preserve results across worker counts in the same R and package environment.

You can also run the scripts from a terminal:

```sh
Rscript mecCox/inst/reproduce/scenario1.R --quick --cores=2
Rscript mecCox/inst/reproduce/scenario2.R --quick --cores=2
Rscript mecCox/inst/reproduce/scenario1.R --cores=20 --replications=100
```

Without a graphical R session, the scripts print their tables to the console;
run them in RStudio to display the HTML report and plots. Command-line
`--replications=N` must be a positive integer and overrides the number of runs
configured in the script. The quick check defaults to two runs per
retained design cell when the configured count remains at 1,000 and no
count is supplied on the command line. A nondefault configured count or an
explicit `--replications` count is retained. Quick mode uses one treated
sample size and a smaller superpopulation reference; Scenario 2 also reduces
the BART and RSF settings and retains all three nonlinearity settings. These
runs check the code path and are not estimates from the paper's simulation.

The [running instructions](inst/reproduce/README.md) describe the
generators, in-memory objects, and display options.

### Public breast-cancer case study

The [breast-cancer script](inst/reproduce/breast_cancer.R) constructs the
GBSG hormonal-therapy target cohort and the Rotterdam external-control
cohort without hormonal therapy. It harmonizes recurrence-free survival
with administrative censoring at five years, then compares unweighted Cox,
ATT-IPW Cox, and MEC-Cox with GLM/Cox and DL/RSF learners. The ATT-IPW rows
report naive, robust sandwich, and corrected sandwich uncertainty estimates
for the same weighted Cox coefficient.

Install the display packages and optional learners, then run the script in
RStudio:

```r
install.packages(c("kableExtra", "htmltools", "rstudioapi", "brulee", "torch", "ranger"))
# If the Torch runtime is not already installed:
torch::install_torch()

script <- system.file("reproduce", "breast_cancer.R", package = "mecCox")
source(script)
```

The two `kableExtra` tables open together in the **Viewer** pane; results
remain in memory, with no result-file exports. The first table reports
log-hazard ratios, standard errors, hazard ratios, and 95% confidence
intervals. The second reports covariate balance and external-control weight
diagnostics. The script generates these analyses in the format of
supplementary Tables S2 and S3; numerical results may vary with random seeds,
learner settings, and software versions.

The R session retains `breast_cancer_estimates`, `breast_cancer_balance`,
`breast_cancer_results`, `breast_cancer_tables`, and `breast_cancer_report`.
See the [case-study instructions](inst/reproduce/README.md#public-breast-cancer-case-study)
for cohort construction and interpretation.

## Methods and interpretation

`fit_att_ipw_cox()` fits a logistic source-propensity model, clips fitted
probabilities to the chosen interval, normalizes external-control odds to sum
to the number of target patients, and fits one Breslow weighted Cox model. It
returns **three standard errors for that one coefficient**: naive model-based,
Lin–Wei/Binder robust, and a corrected sandwich that stacks the Cox and
logistic propensity scores. When clipping is active, the corrected sandwich
uses the fitted logistic probabilities for the logistic score and locally
differentiates the clipped weights used by Cox; interpretation close to a
clipping boundary needs care.

The robust sandwich builds on [Lin and Wei (1989)](https://doi.org/10.1080/01621459.1989.10478874)
and the weighted Cox extension of [Binder (1992)](https://doi.org/10.1093/biomet/79.1.139).
The corrected sandwich follows the stacked-estimating-equation approach to
IPW Cox variance estimation of [Shu, Young, Toh, and Wang (2021)](https://doi.org/10.1111/biom.13332),
applied here to normalized ATT odds.

`fit_mec_cox()` cross-fits both nuisance learners. Its default GLM/Cox fit uses
out-of-fold logistic source probabilities and Cox-predicted control survival.
The alternative MLP/RSF fit uses a neural network for source propensity and a
random survival forest for control survival:

```r
ml_fit <- fit_mec_cox(
  example_external_controls, "time", "event", "source", covariates,
  ps_learner = "mlp", survival_learner = "rsf",
  n_folds = 5, n_landmarks = 5, seed = 2026
)
```

`ps_learner = "bart"` is also available for nonlinear source-propensity
estimation, including the fits in the [Scenario 2 study
script](inst/reproduce/scenario2.R). BART tree count, posterior draws,
burn-in, and shrinkage are explicit arguments. The default API uses fixed
BART settings. Scenario 2 selects `nuisance_settings = "original_study"`
to use the original study's small-subset BART tuning and randomized survival
forests. Its settings are recorded in the run metadata. Completed datasets
are reported during execution and retained in memory after interruption.

MEC-Cox uses positive Kullback–Leibler calibrated weights. The stored
`basis`, `ps_oof`, `survival_oof`, `base_weights`, and `weights` make its
calibration auditable. A failed calibration stops with an error rather than
returning an unbalanced fit. `predicted_risk_balance()` applies the same
out-of-fold survival models at user-specified times and reports the target
minus weighted-external mean predicted mortality. At selected calibration
landmarks the final gap is approximately zero by construction; at other times
it is an empirical diagnostic. Its risk is predicted from baseline covariates,
not an observed death indicator or a patient stratification group.
If weights from a second fit are supplied to `predicted_risk_balance()`, every
column is still evaluated against the *first* fit's stored survival predictions.
This keeps the prognostic score fixed while comparing weighting methods.
`covariate_balance()` provides both raw and standardized differences;
`reported_difference` follows the paper's convention of raw proportion
differences for binary variables and CT-standardized differences for continuous
variables.

The MEC-Cox standard error is a **working stacked sandwich**: it accounts for
estimating the calibration multiplier while treating cross-fitted nuisance
predictions as fixed. The fit records the stacked Jacobian's numerical rank
and reports when a generalized inverse was needed; it never silently replaces
this standard error with a different variance formula. The package does not
claim that calibration alone
establishes a causal effect, valid risk prediction, or improved efficiency in
every data-generating process. The user must assess cohort construction,
overlap, censoring, model fit, and weight stability.

## Development checks

From the parent directory, run:

```sh
R CMD build mecCox
R CMD check mecCox_0.1.1.tar.gz --no-manual
```

The tests compare ATT-IPW output with the manuscript implementation and with
`survival::coxph`, and verify MEC-Cox balance at its fitted landmarks.

## References

Lee, S. Y. (2026). *Balancing Machine-Learned Prognostic Scores to Improve
Efficiency in ATT Marginal Hazard-Ratio Estimation with
Inverse-Probability-Weighted Cox Regression*. **Under review.**

Shu, D., Young, J. G., Toh, S., & Wang, R. (2021).
[Variance estimation in inverse probability weighted Cox models](https://doi.org/10.1111/biom.13332).
*Biometrics*, **77**(3), 1101–1117.

Binder, D. A. (1992).
[Fitting Cox's proportional hazards models from survey data](https://doi.org/10.1093/biomet/79.1.139).
*Biometrika*, **79**(1), 139–147.

Lin, D. Y., & Wei, L. J. (1989).
[The robust inference for the Cox proportional hazards model](https://doi.org/10.1080/01621459.1989.10478874).
*Journal of the American Statistical Association*, **84**(408), 1074–1078.

To obtain the manuscript and methodological references from R:

```r
citation("mecCox")
```

## License and copyright

Copyright © 2026 Se Yoon Lee. Released under the [MIT License](LICENSE.md).
