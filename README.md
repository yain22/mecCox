<img src="man/figures/logo.svg" align="right" alt="mecCox hex logo" width="120" />

# mecCox

[![R package](https://img.shields.io/badge/R-%3E%3D%204.2.0-276DC3?logo=r&logoColor=white)](DESCRIPTION)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)

`mecCox` provides Cox regression methods for analyses with external controls
and code to reproduce the paper's toy example, two main simulations, and
public case study. The treated cohort is coded `source = 1` and external
controls `source = 0`.

## Installation

```r
install.packages("remotes")
remotes::install_github("yain22/mecCox", dependencies = TRUE)
```

## MEC-Cox implementation

Following Supplement S2.4, MEC-Cox constructs cross-fitted survival features
and baseline ATT weights, formulates the calibration problem, and solves its
dual equation by Newton iteration. The calibrated weights then enter the
final weighted Cox regression. The [HTML manual](https://yain22.github.io/mecCox/mec-cox-implementation.html)
follows all six steps with the supplement's notation and executable R code;
the [R Markdown source](docs/mec-cox-implementation.Rmd) is also available.

## Toy analysis

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

The included `example_external_controls` data are simulated.

## Run the toy example, simulations, and public case study

The repository provides runnable code for a toy example with two designs,
two simulation experiments, and one public-data case study:

| Study | Script | Results |
| --- | --- | --- |
| Toy example | [toy_example.R](inst/reproduce/toy_example.R) | Three-panel figures, empirical SDs, and weight diagnostics for both designs |
| Simulation 1: linear source selection and prognosis | [scenario1.R](inst/reproduce/scenario1.R) | Simulation summary and nine-panel figure |
| Simulation 2: increasing nonlinearity | [scenario2.R](inst/reproduce/scenario2.R) | Simulation summary and nine-panel figure |
| Breast-cancer case study | [breast_cancer.R](inst/reproduce/breast_cancer.R) | Hazard-ratio estimates and covariate/weight diagnostics in the format of supplementary Tables S3 and S4 |

Run a script from R or RStudio:

```r
script <- system.file("reproduce", "scenario1.R", package = "mecCox")
source(script)
```

Replace `scenario1.R` with another filename from the table. Simulation run
settings are near the top of each script. Run `breast_cancer.R` in a fresh R
session; it installs the original learner versions and CPU Torch runtime in
a separate library on first use. Reproduction is verified with R 4.5.1 on Windows.

The breast-cancer example uses the public `gbsg` and `rotterdam` data from
`survival`. Restricted SQUIRE/MSK-CHORD patient data are not included.

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
