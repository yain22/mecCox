# Simulation experiments

`scenario1.R` and `scenario2.R` rerun the two main-paper simulation designs
with the exported `mecCox` fitting functions. Install the package first, for
example with `R CMD INSTALL mecCox` from the checkout's parent directory.
Scenario 2 also requires `dbarts` and `ranger`:

```r
install.packages(c("dbarts", "ranger"))
```

## Run from a terminal

```sh
Rscript mecCox/inst/reproduce/scenario1.R --output=scenario1-output
Rscript mecCox/inst/reproduce/scenario2.R --output=scenario2-output
```

For short installation and code-path checks:

```sh
Rscript mecCox/inst/reproduce/scenario1.R --quick --output=scenario1-quick
Rscript mecCox/inst/reproduce/scenario2.R --quick --output=scenario2-quick
```

## Run from R or RStudio

Both scripts also support `source()` and RStudio's **Source** button. The
configuration block near the top of each script contains ordinary R settings:

```r
quick_run <- FALSE
cores <- 20L
output_directory <- "scenario1-output"
```

**These defaults start the full simulation.** To check the installation first,
open `scenario1.R` and edit its configuration block to:

```r
quick_run <- TRUE
cores <- 2L
output_directory <- "scenario1-quick"
```

Save the script, then run the whole file with **Source** or this R command:

```r
source("path/to/mecCox/inst/reproduce/scenario1.R")
```

For Scenario 2, edit the same settings inside `scenario2.R`, choose an output
directory such as `"scenario2-quick"`, and source that file instead. An absolute
path works from any working directory; use forward slashes in R paths on
Windows. Relative output directories are created under the R session's
working directory.

Edit the settings **inside the script** before sourcing it: values assigned
only in the console are replaced by its configuration block. Command-line
options remain available with `Rscript` and override the corresponding script
settings. When running with `source()` or interactively, the script uses its
configuration block rather than unrelated R session command-line arguments.
The helper is found from the sourced file's location, the current project, or
the installed package. If using downloaded copies, keep
`simulation_helpers.R` beside both scenario scripts.

## Parallel execution

The scripts request 20 workers by default and cap the actual worker count
at the number of detected logical cores and replications per design cell.
Choose another limit with `--cores=N`, or use `--cores=1` for serial execution
from a terminal. In R or RStudio, edit `cores` in the script instead:

```sh
Rscript mecCox/inst/reproduce/scenario1.R --cores=8 --output=scenario1-output
```

The PSOCK backend works on Windows, macOS, and Linux. A worker pool is reused
across design cells. Only Monte Carlo replications run in parallel; reference
target computation, checkpoint writing, summaries, and plots run in the main
R process. Each replication receives its own deterministic seed, so changing
worker count or scheduling preserves its data and fits in the same R and
package environment.

## Scenario 1: linear source selection and prognosis

The full run has 1,000 replications at each of 15 combinations of
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
replications and a reference population of 2,000 treated and 4,000 external
controls. It retains all 50 covariates, ten folds, and five landmarks.

## Scenario 2: increasing nonlinearity

The full run uses 10 independent standard-normal covariates, the same five
treated sample sizes, and `n0/n1 = 4`. It has 1,000 replications at each of 15
sample-size and nonlinearity combinations:

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
The reproduction script fixes these BART settings, while the original study
code tuned BART within training folds. It reproduces the study design with
the public API; its fitting sequence and numerical results can differ from
those underlying the paper's figure.

The short `--quick` run keeps all three nonlinearity settings, but uses
`n1 = 200`, two replications per setting, and reference populations of 2,000
treated and 4,000 external controls. It uses 25 BART trees, 50 posterior
draws, 25 burn-in iterations, and 100 RSF trees without tuning. Ten-fold
cross-fitting and five landmarks are retained. It exercises both MEC-Cox
variants. Quick-run results cannot establish simulation performance.

## Output and reproducibility

Each output directory contains:

- `replications.csv`: one row per replication and method, including a visible
  error message for any failed fit.
- `checkpoint_*.csv`: one file written after each completed design
  cell. These preserve finished cells if a long run is interrupted; the script
  does not automatically resume from them.
- `summary.csv`: nominal 95% Wald coverage, Monte Carlo bias and RMSE on the
  **log-hazard-ratio** scale, with counts of successful and failed fits.
- `scenario1.pdf` or `scenario2.pdf`: nine panels arranged by cohort-size
  ratio or nonlinearity setting, respectively, and performance metric.
- `run_metadata.rds`: design, reference target(s), quick/full flag, requested
  and actual worker counts, and R session.
- `reference_targets.csv` in Scenario 2: the reference log-hazard ratio for
  each nonlinearity setting.

The deterministic seeds make each script rerunnable, but the package API,
learner tuning, and replication streams are not a bit-for-bit replay of the
earlier private parallel code used to produce the paper's figures. The
package's corrected ATT-IPW sandwich also differentiates the clipped analysis
weights locally;
its values can differ from an earlier correction when fitted probabilities
cross a clipping boundary. Results should be compared at the Monte Carlo
level; a quick run cannot establish performance.
Before interpreting a full run, check `summary.csv` for failures and inspect
their messages in `replications.csv`. These scripts do not infer or replace
missing results when a fit fails.
