# The helpers in this file implement the Breslow weighted Cox score used by
# the ATT-IPW comparator.  They are deliberately kept separate from the
# MEC-Cox calibration code: all three ATT-IPW intervals use one point estimate.

.att_clip_probability <- function(probability, limits) {
  pmin(pmax(probability, limits[1]), limits[2])
}

.att_weights_for_gamma <- function(gamma, design_matrix, source, limits) {
  linear_predictor <- as.vector(design_matrix %*% gamma)
  fitted_probability <- stats::plogis(linear_predictor)
  clipped_probability <- .att_clip_probability(fitted_probability, limits)
  odds <- clipped_probability / (1 - clipped_probability)

  control <- source == 0L
  target_count <- sum(source == 1L)
  control_odds_total <- sum(odds[control])

  if (!is.finite(control_odds_total) || control_odds_total <= 0) {
    stop("The external-control odds cannot be normalized.", call. = FALSE)
  }

  normalization_factor <- target_count / control_odds_total
  weights <- rep(1, length(source))
  weights[control] <- normalization_factor * odds[control]

  list(
    weights = weights,
    raw_probability = fitted_probability,
    clipped_probability = clipped_probability,
    normalization_factor = normalization_factor
  )
}

.att_cox_score <- function(theta, prepared, weights) {
  source <- prepared$A
  observed_time <- prepared$time
  event <- prepared$delta
  sample_size <- nrow(prepared)

  relative_risk <- exp(theta * source)
  weighted_risk <- weights * relative_risk

  # Cumulative sums in descending time order give the Breslow risk-set totals
  # at every distinct observed time, including tied event and censoring times.
  descending_order <- order(observed_time, decreasing = TRUE)
  descending_time <- observed_time[descending_order]
  time_groups <- rle(descending_time)
  group_end <- cumsum(time_groups$lengths)

  risk_total <- cumsum(weighted_risk[descending_order])[group_end]
  treated_risk_total <- cumsum(
    (weighted_risk * source)[descending_order]
  )[group_end]

  time_position <- match(observed_time, time_groups$values)
  risk_mean_at_subject <- treated_risk_total[time_position] /
    risk_total[time_position]
  direct_contribution <- weights * event * (source - risk_mean_at_subject)

  event_times <- sort(unique(observed_time[event == 1L]))
  if (length(event_times) == 0L) {
    stop("At least one event is required to fit the Cox model.", call. = FALSE)
  }

  # The event mass is the sum of case weights at a tied event time.  This is
  # the Breslow convention used by survival::coxph(ties = "breslow").
  event_mass <- vapply(event_times, function(event_time) {
    sum(weights[event == 1L & observed_time == event_time])
  }, numeric(1))

  event_position <- match(event_times, time_groups$values)
  event_risk_total <- risk_total[event_position]
  event_risk_mean <- treated_risk_total[event_position] /
    event_risk_total

  cumulative_event_hazard <- cumsum(event_mass / event_risk_total)
  cumulative_centered_hazard <- cumsum(
    event_mass * event_risk_mean / event_risk_total
  )

  last_event_position <- findInterval(observed_time, event_times)
  hazard_at_subject <- numeric(sample_size)
  centered_hazard_at_subject <- numeric(sample_size)
  has_prior_event <- last_event_position > 0L
  hazard_at_subject[has_prior_event] <- cumulative_event_hazard[
    last_event_position[has_prior_event]
  ]
  centered_hazard_at_subject[has_prior_event] <-
    cumulative_centered_hazard[last_event_position[has_prior_event]]

  risk_set_contribution <- weighted_risk * (
    source * hazard_at_subject - centered_hazard_at_subject
  )
  individual_contribution <- direct_contribution - risk_set_contribution

  list(
    score = sum(direct_contribution),
    individual_contribution = individual_contribution
  )
}

.att_central_difference <- function(function_to_differentiate, parameter) {
  derivative <- numeric(length(parameter))

  for (column in seq_along(parameter)) {
    step_size <- 1e-5 * max(1, abs(parameter[column]))
    increased <- parameter
    decreased <- parameter
    increased[column] <- increased[column] + step_size
    decreased[column] <- decreased[column] - step_size

    derivative[column] <- (
      function_to_differentiate(increased) -
        function_to_differentiate(decreased)
    ) / (2 * step_size)
  }

  derivative
}

#' Fit an ATT-weighted Cox model with three variance estimates
#'
#' The target/source-1 cohort receives unit weight. External controls receive
#' propensity-score odds weights, normalized to sum to the number of treated
#' patients. One Breslow weighted Cox fit supplies the log hazard ratio; the
#' naive, Lin--Wei/Binder robust, and corrected sandwich calculations differ
#' only in their standard errors. The corrected calculation follows the
#' manuscript's joint Cox/logistic-score convention. When probabilities are
#' clipped, its logistic score and information use the *unclipped* fitted GLM
#' probabilities, while the Cox-to-logistic derivative follows the clipped,
#' normalized weights actually used by the estimator. Clipping introduces a
#' nondifferentiable boundary; the sandwich is a local approximation away from
#' that boundary.
#'
#' @param data A data frame with one row per patient.
#' @param time,event,source Names of the follow-up time, event indicator, and
#'   source indicator columns. `source` is one for the target cohort and zero for external
#'   control; `event` is one for an observed event.
#' @param covariates Character vector of baseline covariate column names for
#'   the logistic source-propensity model.
#' @param ps_clip Two probabilities defining the clipping interval.
#' @param conf_level Confidence level for two-sided Wald intervals.
#'
#' @return An `att_ipw_cox_fit` list with the common `theta` and `hr`, named
#'   `se` values and `conf_int` intervals for all three variance methods,
#'   patient weights, raw and clipped source-propensity scores, model fits,
#'   and weight diagnostics. The returned weights retain the input row order.
#' @export
#' @examples
#' data(example_external_controls)
#' att_fit <- fit_att_ipw_cox(
#'   example_external_controls,
#'   time = "time", event = "event", source = "source",
#'   covariates = c("age", "sex", "marker")
#' )
#' summary(att_fit)
fit_att_ipw_cox <- function(data, time, event, source, covariates,
                            ps_clip = c(0.01, 0.99), conf_level = 0.95) {
  if (!is.numeric(ps_clip) || length(ps_clip) != 2L ||
      anyNA(ps_clip) || any(!is.finite(ps_clip)) ||
      ps_clip[1] <= 0 || ps_clip[1] >= ps_clip[2] ||
      ps_clip[2] >= 1) {
    stop("ps_clip must satisfy 0 < lower < upper < 1.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must lie strictly between zero and one.", call. = FALSE)
  }

  prepared <- .pkg_validate_data(data, time, event, source, covariates)
  source_formula <- stats::reformulate(covariates, response = "A")
  propensity_fit <- stats::glm(
    source_formula,
    data = prepared,
    family = stats::binomial()
  )

  gamma <- stats::coef(propensity_fit)
  if (!isTRUE(propensity_fit$converged) || any(!is.finite(gamma))) {
    stop("The logistic source-propensity model did not converge with finite, identifiable coefficients.",
         call. = FALSE)
  }

  design_matrix <- stats::model.matrix(propensity_fit)
  weight_result <- .att_weights_for_gamma(
    gamma, design_matrix, prepared$A, ps_clip
  )
  weights <- weight_result$weights
  prepared$.att_weight <- weights

  cox_formula <- stats::as.formula("survival::Surv(time, delta) ~ A")
  cox_fit <- survival::coxph(
    cox_formula,
    data = prepared,
    weights = prepared$.att_weight,
    robust = FALSE,
    ties = "breslow",
    x = TRUE,
    control = survival::coxph.control(timefix = FALSE)
  )
  theta <- unname(stats::coef(cox_fit)["A"])
  if (!is.finite(theta)) {
    stop("The weighted Cox coefficient is not finite.", call. = FALSE)
  }

  naive_standard_error <- sqrt(as.numeric(stats::vcov(cox_fit)["A", "A"]))
  cox_score <- .att_cox_score(theta, prepared, weights)
  theta_derivative <- .att_central_difference(
    function(candidate) {
      .att_cox_score(candidate[1], prepared, weights)$score
    },
    theta
  )[1]

  if (!is.finite(theta_derivative) || abs(theta_derivative) < 1e-10) {
    stop("The weighted Cox score has insufficient information.", call. = FALSE)
  }

  robust_influence <- -cox_score$individual_contribution / theta_derivative
  robust_standard_error <- sqrt(sum(robust_influence^2))

  # Differentiate the *normalized* odds as gamma changes. This reproduces the
  # corrected-variance convention in the manuscript when clipping is inactive.
  cross_derivative <- .att_central_difference(
    function(candidate_gamma) {
      candidate_weights <- .att_weights_for_gamma(
        candidate_gamma, design_matrix, prepared$A, ps_clip
      )$weights
      .att_cox_score(theta, prepared, candidate_weights)$score
    },
    gamma
  )

  # The propensity model is fitted by glm before the clipping transformation.
  # Its score and information must therefore use the original GLM probabilities.
  raw_probability <- weight_result$raw_probability
  logistic_score <- sweep(
    design_matrix, 1L, prepared$A - raw_probability, "*"
  )
  logistic_information <- crossprod(
    design_matrix * (raw_probability * (1 - raw_probability)),
    design_matrix
  )

  parameter_count <- length(gamma)
  derivative_matrix <- matrix(
    0, nrow = parameter_count + 1L, ncol = parameter_count + 1L
  )
  derivative_matrix[1, 1] <- theta_derivative
  derivative_matrix[1, -1] <- cross_derivative
  derivative_matrix[-1, -1] <- -logistic_information

  if (qr(derivative_matrix, tol = 1e-10)$rank < ncol(derivative_matrix)) {
    stop("The corrected sandwich derivative matrix is singular.", call. = FALSE)
  }
  inverse_derivative <- solve(derivative_matrix)
  individual_scores <- cbind(
    cox_score$individual_contribution,
    logistic_score
  )
  corrected_influence <- -individual_scores %*% t(inverse_derivative)
  corrected_standard_error <- sqrt(sum(corrected_influence[, 1]^2))

  standard_errors <- c(
    naive = naive_standard_error,
    robust = robust_standard_error,
    corrected = corrected_standard_error
  )
  critical_value <- stats::qnorm((1 + conf_level) / 2)
  confidence_intervals <- data.frame(
    method = names(standard_errors),
    log_lower = theta - critical_value * standard_errors,
    log_upper = theta + critical_value * standard_errors,
    hazard_lower = exp(theta - critical_value * standard_errors),
    hazard_upper = exp(theta + critical_value * standard_errors),
    row.names = NULL
  )

  external_weights <- weights[prepared$A == 0L]
  diagnostics <- list(
    target_count = sum(prepared$A == 1L),
    external_count = length(external_weights),
    external_weight_sum = sum(external_weights),
    external_effective_sample_size =
      sum(external_weights)^2 / sum(external_weights^2),
    external_weight_cv = stats::sd(external_weights) /
      mean(external_weights),
    normalization_factor = weight_result$normalization_factor,
    clipped_probability_count = sum(
      weight_result$raw_probability != weight_result$clipped_probability
    )
  )

  result <- list(
    theta = theta,
    hr = exp(theta),
    se = standard_errors,
    conf_int = confidence_intervals,
    weights = weights,
    ps = weight_result$raw_probability,
    ps_clipped = weight_result$clipped_probability,
    diagnostics = diagnostics,
    analysis_data = prepared,
    covariates = covariates,
    ps_fit = propensity_fit,
    cox_fit = cox_fit,
    method = "ATT-IPW Cox",
    call = match.call()
  )
  class(result) <- "att_ipw_cox_fit"
  result
}
