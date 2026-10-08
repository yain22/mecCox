#' Inspect baseline covariate balance
#'
#' Compares the target cohort with external controls before weighting and under
#' each set of fitted control weights. Factors are expanded into indicator
#' columns. The `difference` column is always target minus weighted external.
#' The `standardized_difference` column divides that difference by the
#' unweighted target standard deviation. The `reported_difference` column
#' follows the accompanying paper's balance-table convention: binary
#' indicators use the raw proportion difference; continuous variables use
#' the target-standardized difference. A zero target standard deviation
#' produces `NA` for a standardized or reported continuous difference.
#'
#' @param object A result from [fit_att_ipw_cox()] or [fit_mec_cox()].
#'
#' @return A data frame with one row per expanded covariate and weight set,
#'   containing target and external means, their raw, standardized, and
#'   paper-convention differences.
#' @export
covariate_balance <- function(object) {
  if (!inherits(object, c("att_ipw_cox_fit", "mec_cox_fit"))) {
    stop("`object` must be an ATT-IPW or MEC-Cox fit.", call. = FALSE)
  }

  data <- object$analysis_data
  covariates <- object$covariates
  design <- .pkg_covariate_matrix(data, covariates)
  if (ncol(design) == 0L) {
    stop("No baseline covariates are available in this fit.", call. = FALSE)
  }

  target <- data$A == 1L
  external <- data$A == 0L
  target_means <- colMeans(design[target, , drop = FALSE])
  target_sd <- apply(design[target, , drop = FALSE], 2L, stats::sd)
  binary_column <- apply(design, 2L, function(values) {
    all(values %in% c(0, 1))
  })

  weights <- list(unweighted = rep(1, nrow(data)))
  if (inherits(object, "mec_cox_fit")) {
    weights$source_odds <- object$base_weights
    weights$mec <- object$weights
  } else {
    weights$att_ipw <- object$weights
  }

  result <- vector("list", length(weights))
  for (index in seq_along(weights)) {
    control_weights <- weights[[index]][external]
    external_means <- colSums(
      design[external, , drop = FALSE] * control_weights
    ) / sum(control_weights)
    difference <- target_means - external_means

    standardized <- rep(NA_real_, length(difference))
    estimable <- is.finite(target_sd) & target_sd > 0
    standardized[estimable] <- difference[estimable] / target_sd[estimable]
    reported <- standardized
    reported[binary_column] <- difference[binary_column]

    result[[index]] <- data.frame(
      covariate = colnames(design),
      weighting = names(weights)[index],
      variable_type = ifelse(binary_column, "binary", "continuous"),
      target_mean = as.numeric(target_means),
      external_mean = as.numeric(external_means),
      difference = as.numeric(difference),
      standardized_difference = standardized,
      reported_difference = reported,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, result)
}

#' Evaluate the predicted mortality-risk gap
#'
#' Applies the same out-of-fold survival models from a MEC-Cox fit to every
#' patient at the requested times. For each time, the reported gap is the
#' target mean predicted risk minus the weighted external-control mean
#' predicted risk (Equation 32 in the accompanying paper). It is a balance
#' diagnostic, not a Kaplan--Meier estimate or a stratification variable.
#'
#' The fitted calibration constraints imply an approximately zero MEC gap at
#' retained landmark times. At other times, including integer months between
#' landmarks, the gap is descriptive and need not be zero.
#'
#' @param object A result from [fit_mec_cox()].
#' @param times Positive follow-up times in the units used for fitting.
#' @param additional_weights Optional named list of full-length weight vectors
#'   aligned with the rows of the fitted data. All columns use the prediction
#'   models stored in `object`, including those from additional weight sets.
#'   This allows, for example, comparison with a Cox-basis MEC fit using the
#'   same patients while retaining RSF risk predictions as the diagnostic.
#'
#' @return A data frame with one row per time and weighting method. Risk and
#'   gap columns are on the probability scale; multiply by 100 for percentage
#'   points.
#' @export
predicted_risk_balance <- function(object, times,
                                   additional_weights = list()) {
  if (!inherits(object, "mec_cox_fit")) {
    stop("`object` must be a MEC-Cox fit.", call. = FALSE)
  }
  if (!is.numeric(times) || length(times) == 0L ||
      anyNA(times) || any(!is.finite(times)) || any(times <= 0)) {
    stop("`times` must contain positive, finite follow-up times.",
         call. = FALSE)
  }
  if (!is.list(additional_weights)) {
    stop("`additional_weights` must be a named list.", call. = FALSE)
  }
  if (length(additional_weights) > 0L &&
      (is.null(names(additional_weights)) ||
       anyNA(names(additional_weights)) ||
       any(names(additional_weights) == ""))) {
    stop("Every `additional_weights` vector must have a name.",
         call. = FALSE)
  }

  patient_count <- nrow(object$analysis_data)
  target <- object$analysis_data$A == 1L
  external <- object$analysis_data$A == 0L
  weight_sets <- c(
    list(unweighted = rep(1, patient_count),
         source_odds = object$base_weights,
         mec = object$weights),
    additional_weights
  )
  if (anyDuplicated(names(weight_sets))) {
    stop("Weighting-method names must be unique.", call. = FALSE)
  }
  for (name in names(weight_sets)) {
    vector <- weight_sets[[name]]
    if (!is.numeric(vector) || length(vector) != patient_count ||
        anyNA(vector) || any(!is.finite(vector[external])) ||
        any(vector[external] <= 0)) {
      stop("Weights for `", name,
           "` must be a full-length numeric vector with positive, finite external-control weights.",
           call. = FALSE)
    }
  }

  times <- as.numeric(times)
  survival <- .mec_oof_survival_at_times(object, times)
  landmark_index <- match(times, object$landmark_times)
  at_landmark <- which(!is.na(landmark_index))
  if (length(at_landmark) > 0L) {
    survival[, at_landmark] <- object$survival_oof[
      , landmark_index[at_landmark], drop = FALSE
    ]
  }
  risk <- 1 - survival

  result <- vector("list", length(times) * length(weight_sets))
  row <- 0L
  for (time_index in seq_along(times)) {
    target_mean <- mean(risk[target, time_index])
    for (name in names(weight_sets)) {
      control_weights <- weight_sets[[name]][external]
      external_mean <- stats::weighted.mean(
        risk[external, time_index], control_weights
      )
      gap <- target_mean - external_mean
      row <- row + 1L
      result[[row]] <- data.frame(
        time = times[time_index],
        weighting = name,
        target_risk = target_mean,
        external_risk = external_mean,
        gap = gap,
        absolute_gap = abs(gap),
        external_effective_sample_size =
          sum(control_weights)^2 / sum(control_weights^2),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, result)
}
