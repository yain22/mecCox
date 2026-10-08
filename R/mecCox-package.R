#' mecCox: prognostic balance for ATT-weighted Cox regression
#'
#' The package fits unweighted, ATT inverse-probability-weighted, and MEC-Cox
#' marginal source comparisons. The target cohort has `source = 1` and unit
#' weight; the external-control cohort has `source = 0`. MEC-Cox uses
#' cross-fitted estimates of source propensity and control survival, then
#' calibrates external-control odds weights to match selected prognostic
#' predictions in the target cohort.
#'
#' The returned Cox hazard ratio is a marginal source-comparison summary.
#' Its causal interpretation requires assumptions about cohort comparability,
#' consistency, positivity, and censoring that package software cannot verify.
#' The MEC-Cox stacked sandwich treats fitted nuisance learners as fixed;
#' cross-fitting alone does not establish unconditional Wald coverage.
#'
#' @keywords internal
"_PACKAGE"
