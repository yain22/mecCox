# Data checks shared by the three Cox estimators. The fitted models are always
# expressed using A = 1 for the target/trial cohort and A = 0 for the external
# control cohort.
.pkg_validate_data <- function(data, time, event, source, covariates) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (anyDuplicated(names(data))) {
    stop("`data` must have unique column names.", call. = FALSE)
  }

  column_arguments <- list(time = time, event = event, source = source)
  valid_column_arguments <- vapply(column_arguments, function(value) {
    is.character(value) && length(value) == 1L &&
      !is.na(value) && nzchar(value)
  }, logical(1))
  if (!all(valid_column_arguments)) {
    stop("`time`, `event`, and `source` must each name one column.",
         call. = FALSE)
  }
  if (anyDuplicated(c(time, event, source))) {
    stop("`time`, `event`, and `source` must name different columns.",
         call. = FALSE)
  }

  if (!is.character(covariates) || anyNA(covariates)) {
    stop("`covariates` must be a character vector of column names.",
         call. = FALSE)
  }
  if (anyDuplicated(covariates)) {
    stop("`covariates` contains duplicate names.", call. = FALSE)
  }
  if (any(make.names(covariates) != covariates)) {
    stop("Covariate names must be syntactic R names; rename columns containing spaces or punctuation.",
         call. = FALSE)
  }
  if (any(covariates %in% c(time, event, source, "time", "delta", "A"))) {
    stop("Covariates cannot be a time/event/source column or use reserved names `time`, `delta`, or `A`.",
         call. = FALSE)
  }

  required <- c(time, event, source, covariates)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0L) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }
  if (nrow(data) < 4L) {
    stop("`data` must contain at least four patients.", call. = FALSE)
  }

  observed_time <- data[[time]]
  if (!is.numeric(observed_time) ||
      anyNA(observed_time) ||
      any(!is.finite(observed_time)) ||
      any(observed_time <= 0)) {
    stop("Follow-up times must be finite, positive numbers.", call. = FALSE)
  }

  observed_event <- data[[event]]
  observed_source <- data[[source]]
  if (!.pkg_is_binary(observed_event)) {
    stop("`event` must contain only 0 and 1, with no missing values.",
         call. = FALSE)
  }
  if (!.pkg_is_binary(observed_source)) {
    stop("`source` must contain only 0 and 1, with no missing values.",
         call. = FALSE)
  }
  if (length(unique(observed_source)) != 2L) {
    stop("Both source groups (0 and 1) must be present.", call. = FALSE)
  }
  event_count_by_source <- tapply(
    as.integer(observed_event), as.integer(observed_source), sum
  )
  if (any(event_count_by_source == 0L)) {
    stop("At least one observed event is required in each source group.",
         call. = FALSE)
  }

  prepared <- data.frame(
    time = as.numeric(observed_time),
    delta = as.integer(observed_event),
    A = as.integer(observed_source)
  )

  for (name in covariates) {
    values <- data[[name]]
    if (is.character(values)) {
      values <- factor(values)
    }
    if (is.factor(values)) {
      if (anyNA(values) || nlevels(droplevels(values)) < 2L) {
        stop("Covariate `", name, "` must have at least two observed levels and no missing values.",
             call. = FALSE)
      }
      values <- droplevels(values)
    } else if (is.numeric(values) || is.logical(values)) {
      if (anyNA(values) || any(!is.finite(as.numeric(values)))) {
        stop("Covariate `", name, "` must have only finite, nonmissing values.",
             call. = FALSE)
      }
    } else {
      stop("Covariate `", name, "` must be numeric, logical, factor, or character.",
           call. = FALSE)
    }
    prepared[[name]] <- values
  }

  prepared
}

.pkg_is_binary <- function(values) {
  if (!(is.numeric(values) || is.logical(values)) || anyNA(values)) {
    return(FALSE)
  }
  all(values %in% c(0, 1))
}

.pkg_covariate_matrix <- function(prepared, covariates) {
  if (length(covariates) == 0L) {
    return(matrix(numeric(), nrow = nrow(prepared), ncol = 0L))
  }

  covariate_frame <- prepared[, covariates, drop = FALSE]
  design <- stats::model.matrix(~ ., data = covariate_frame)
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  storage.mode(design) <- "double"
  design
}

.pkg_validate_conf_level <- function(conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a number between 0 and 1.", call. = FALSE)
  }
  invisible(conf_level)
}
