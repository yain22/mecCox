# First simulation experiment

`scenario1.R` reruns the first main-paper simulation with the exported
`mecCox` fitting functions. It is intentionally limited to Scenario 1. Run it
after installing the package, for example with `R CMD INSTALL mecCox` from
the checkout's parent directory.

```sh
Rscript mecCox/inst/reproduce/scenario1.R --output=scenario1-output
```

For a short installation/code-path check:

```sh
Rscript mecCox/inst/reproduce/scenario1.R --quick --output=scenario1-quick
```

The complete run has 1,000 replications at each of 15 combinations of
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

The output directory contains:

- `replications.csv`: one row per replication and method, including a visible
  error message for any failed fit.
- `checkpoint_n1-*_n0-*.csv`: one file written after each completed design
  cell. These preserve finished cells if a long run is interrupted; the script
  does not automatically resume from them.
- `summary.csv`: nominal 95% Wald coverage, Monte Carlo bias and RMSE on the
  **log-hazard-ratio** scale, with counts of successful and failed fits.
- `scenario1.pdf`: nine panels arranged by cohort-size ratio and metric.
- `run_metadata.rds`: design, reference target, quick/full flag, and R session.

The deterministic seeds make this script rerunnable, but the package API and
replication stream are not a bit-for-bit replay of the earlier private
parallel code used to produce the published figure. The package's corrected
ATT-IPW sandwich also differentiates the clipped analysis weights locally;
its values can differ from an earlier correction when fitted probabilities
cross a clipping boundary. Results should be compared at the Monte Carlo
level; a quick run cannot establish performance.
Before interpreting a full run, check `summary.csv` for failures and inspect
their messages in `replications.csv`. This script does not infer or replace
missing results when a fit fails.
