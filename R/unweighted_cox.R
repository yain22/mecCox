#' Unweighted marginal Cox comparison
#'
#' Fits the two-source Cox reference model without inverse-probability or
#' calibration weights. Source 1 is the target cohort and source 0 is the
#' external-control cohort. The coefficient is the source-1 versus source-0
#' log hazard ratio.
#'
#' @param data A data frame containing the follow-up time, event indicator,
#'   and source indicator.
#' @param time,event,source Character names of columns in `data`. `event` and
#'   `source` must contain 0/1 values.
#' @param conf_level Confidence level for the Wald interval.
#'
#' @return An object of class `unweighted_cox_fit` with the log hazard ratio,
#'   hazard ratio, standard error, confidence interval, and fitted Cox model.
#' @export
#' @examples
#' data(example_external_controls)
#' fit_unweighted_cox(example_external_controls, "time", "event", "source")
fit_unweighted_cox <- function(data, time, event, source, conf_level = 0.95) {
  .pkg_validate_conf_level(conf_level)
  analysis_data <- .pkg_validate_data(
    data = data,
    time = time,
    event = event,
    source = source,
    covariates = character()
  )

  model <- survival::coxph(
    survival::Surv(time, delta) ~ A,
    data = analysis_data,
    ties = "breslow",
    robust = FALSE,
    control = survival::coxph.control(timefix = FALSE)
  )
  theta <- unname(stats::coef(model)[["A"]])
  se <- sqrt(unname(stats::vcov(model)["A", "A"]))
  if (!is.finite(theta) || !is.finite(se) || se <= 0) {
    stop("The unweighted Cox coefficient or standard error is not finite.",
         call. = FALSE)
  }

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  interval <- exp(theta + c(-1, 1) * z * se)
  names(interval) <- c("lower", "upper")

  result <- list(
    theta = theta,
    se = se,
    hr = exp(theta),
    conf_int = interval,
    conf_level = conf_level,
    weights = rep(1, nrow(analysis_data)),
    model = model,
    call = match.call()
  )
  class(result) <- "unweighted_cox_fit"
  result
}
