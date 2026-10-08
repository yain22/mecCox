#' Print an unweighted Cox comparison
#'
#' @param x An object returned by [fit_unweighted_cox()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.unweighted_cox_fit <- function(x, ...) {
  cat("Unweighted Cox source comparison\n")
  cat(sprintf("  Hazard ratio: %.4f\n", x$hr))
  cat(sprintf("  %.0f%% Wald interval: [%.4f, %.4f]\n",
              100 * x$conf_level,
              x$conf_int[["lower"]], x$conf_int[["upper"]]))
  cat(sprintf("  Log-HR standard error: %.4f\n", x$se))
  invisible(x)
}

#' Print an ATT-IPW Cox comparison
#'
#' @param x An object returned by [fit_att_ipw_cox()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.att_ipw_cox_fit <- function(x, ...) {
  cat("ATT-IPW Cox source comparison\n")
  cat(sprintf("  Hazard ratio: %.4f (one common estimate)\n", x$hr))
  for (method in names(x$se)) {
    cat(sprintf("  %-10s standard error: %.4f\n", method, x$se[[method]]))
  }
  cat(sprintf("  External-control ESS: %.1f of %d\n",
              x$diagnostics$external_effective_sample_size,
              x$diagnostics$external_count))
  invisible(x)
}

#' Print a MEC-Cox comparison
#'
#' @param x An object returned by [fit_mec_cox()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.mec_cox_fit <- function(x, ...) {
  cat(sprintf("MEC-Cox source comparison (%s/%s)\n",
              x$ps_learner, x$survival_learner))
  cat(sprintf("  Hazard ratio: %.4f\n", x$hr))
  cat(sprintf("  %.0f%% working Wald interval: [%.4f, %.4f]\n",
              100 * x$conf_level,
              x$conf_int[["lower"]], x$conf_int[["upper"]]))
  cat(sprintf("  Working stacked-sandwich SE: %.4f\n", x$se))
  cat(sprintf("  External-control ESS: %.1f\n",
              x$weight_diagnostics$ess))
  cat(sprintf("  Largest mean calibration residual: %.2g\n",
              x$calibration$max_mean_balance_residual))
  if (isTRUE(x$variance_diagnostics$used_pseudoinverse)) {
    cat("  Variance calculation used a generalized inverse of the stacked Jacobian.\n")
  }
  invisible(x)
}

#' Summarize a Cox comparison
#'
#' @param object An object returned by [fit_unweighted_cox()],
#'   [fit_att_ipw_cox()], or [fit_mec_cox()].
#' @param ... Unused.
#' @return A data frame containing the log hazard ratio, standard error,
#'   hazard ratio, and confidence limits. The ATT-IPW method has one row for
#'   each of its three variance estimates.
#' @name summary.mecCox
NULL

#' @rdname summary.mecCox
#' @export
summary.unweighted_cox_fit <- function(object, ...) {
  data.frame(
    method = "unweighted",
    log_hr = object$theta,
    se = object$se,
    hr = object$hr,
    hr_lower = object$conf_int[["lower"]],
    hr_upper = object$conf_int[["upper"]]
  )
}

#' @rdname summary.mecCox
#' @export
summary.att_ipw_cox_fit <- function(object, ...) {
  data.frame(
    method = paste0("att_ipw_", names(object$se)),
    log_hr = rep(object$theta, length(object$se)),
    se = as.numeric(object$se),
    hr = rep(object$hr, length(object$se)),
    hr_lower = object$conf_int$hazard_lower,
    hr_upper = object$conf_int$hazard_upper
  )
}

#' @rdname summary.mecCox
#' @export
summary.mec_cox_fit <- function(object, ...) {
  data.frame(
    method = paste0("mec_", object$ps_learner, "_", object$survival_learner),
    log_hr = object$theta,
    se = object$se,
    hr = object$hr,
    hr_lower = object$conf_int[["lower"]],
    hr_upper = object$conf_int[["upper"]]
  )
}
