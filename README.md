# mecCox

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

## Reproduce the first simulation experiment

The [Scenario 1 script](inst/reproduce/scenario1.R) regenerates the paper's
linear source-selection and outcome experiment. It compares the three
ATT-IPW Cox standard errors with GLM/Cox MEC-Cox for each treated sample size
`n1 = 200, 250, 300, 350, 400` and external-control ratio `n1:n0 = 1:2, 1:3,
1:4`. It uses 50 covariates, five event-time-quantile survival landmarks,
ten-fold cross-fitting, and 1,000 Monte Carlo replications per sample-size
combination. The output includes the replication-level estimates, a table of
coverage, bias, and RMSE, and a nine-panel PDF figure. No other simulation
scenario is included in this reproduction entry point.

After installing the package, run from the checkout's parent directory:

```sh
Rscript mecCox/inst/reproduce/scenario1.R --output=scenario1-output
```

The full run fits 15,000 simulated datasets and can take substantial time. To
check that the code and dependencies work before starting it, run:

```sh
Rscript mecCox/inst/reproduce/scenario1.R --quick --output=scenario1-quick
```

`--quick` uses two replications, one sample-size combination, and a smaller
superpopulation reference. Its numerical results are only a code-path check;
they are not estimates from the paper's simulation. See the
[reproduction notes](inst/reproduce/README.md) for the exact generator and
interpretation of the output.

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
estimation, including the class of fits studied in the paper's second
simulation scenario. BART tree count, posterior draws, burn-in, and shrinkage
are explicit arguments. The package does not automatically tune BART; users
who need exact simulation settings should consult the study's simulation
script and supply its selected hyperparameters.

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
R CMD check mecCox_0.1.0.tar.gz --no-manual
```

The tests compare ATT-IPW output with the manuscript implementation and with
`survival::coxph`, and verify MEC-Cox balance at its fitted landmarks.
