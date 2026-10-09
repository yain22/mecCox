# MEC-Cox: cross-fitted prognostic calibration for an ATT-weighted Cox fit.
#
# A = 1 is the target cohort and retains weight one. A = 0 supplies the
# external controls. The control odds weights are calibrated so that selected
# out-of-fold predictions of control survival have the same weighted totals in
# the two cohorts.

.mec_with_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed))
  force(expression)
}

.mec_make_folds <- function(source, n_folds, seed) {
  .mec_with_seed(seed, {
    fold_id <- integer(length(source))
    for (group in c(0L, 1L)) {
      rows <- which(source == group)
      fold_id[rows] <- sample(rep(seq_len(n_folds), length.out = length(rows)))
    }
    fold_id
  })
}

.mec_landmarks <- function(data, n_landmarks, landmarks) {
  if (!is.null(landmarks)) {
    if (!is.numeric(landmarks) || !length(landmarks) ||
        any(!is.finite(landmarks)) || any(landmarks <= 0)) {
      stop("`landmarks` must contain positive, finite times.", call. = FALSE)
    }
    if (anyDuplicated(landmarks)) {
      stop("`landmarks` must not contain duplicates.", call. = FALSE)
    }
    return(sort(as.numeric(landmarks)))
  }

  if (!is.numeric(n_landmarks) || length(n_landmarks) != 1L ||
      !is.finite(n_landmarks) ||
      n_landmarks < 1L || n_landmarks != as.integer(n_landmarks)) {
    stop("`n_landmarks` must be a positive integer.", call. = FALSE)
  }

  event_times <- data$time[data$A == 0L & data$delta == 1L]
  if (!length(event_times)) {
    stop("External controls must have observed events.", call. = FALSE)
  }
  probabilities <- seq(0.10, 0.90, length.out = n_landmarks)
  selected <- as.numeric(stats::quantile(
    event_times, probabilities, names = FALSE, type = 8
  ))
  selected <- sort(unique(selected))
  if (length(selected) < n_landmarks) {
    warning("Repeated event-time quantiles reduced the number of landmarks.",
            call. = FALSE)
  }
  selected
}

.mec_check_factor_levels <- function(training, validation, covariates,
                                     learner) {
  for (variable in covariates) {
    if (!is.factor(training[[variable]])) {
      next
    }
    training_levels <- unique(as.character(training[[variable]]))
    validation_levels <- unique(as.character(validation[[variable]]))
    if (length(training_levels) < 2L) {
      stop(sprintf(
        "%s training data have fewer than two observed levels of `%s`.",
        learner, variable
      ), call. = FALSE)
    }
    if (length(setdiff(validation_levels, training_levels))) {
      stop(sprintf(
        "%s validation data contain a level of `%s` unseen in training.",
        learner, variable
      ), call. = FALSE)
    }
  }
  invisible(NULL)
}

.mec_propensity_glm <- function(training, validation, covariates) {
  .mec_check_factor_levels(training, validation, covariates,
                           "Propensity GLM")
  propensity_formula <- stats::reformulate(covariates, response = "A")
  fit <- stats::glm(propensity_formula, data = training,
                    family = stats::binomial())
  if (!isTRUE(fit$converged) ||
      any(!is.finite(stats::coef(fit)))) {
    stop("The propensity GLM did not converge to finite coefficients.",
         call. = FALSE)
  }
  as.numeric(stats::predict(fit, newdata = validation, type = "response"))
}

.mec_propensity_mlp <- function(training, validation, covariates, seed,
                                hidden_units, epochs) {
  if (!requireNamespace("brulee", quietly = TRUE)) {
    stop("Install {brulee} to use `ps_learner = 'mlp'`.", call. = FALSE)
  }
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("Install {torch} to use `ps_learner = 'mlp'`.", call. = FALSE)
  }
  .mec_check_factor_levels(training, validation, covariates,
                           "Propensity MLP")

  design_formula <- stats::reformulate(covariates)
  x_train <- stats::model.matrix(design_formula, training)[, -1L, drop = FALSE]
  x_valid <- stats::model.matrix(design_formula, validation)[, -1L, drop = FALSE]
  if (!identical(colnames(x_train), colnames(x_valid))) {
    stop("Covariate columns changed between an MLP training and validation fold.",
         call. = FALSE)
  }

  center <- colMeans(x_train)
  scale <- apply(x_train, 2L, stats::sd)
  scale[!is.finite(scale) | scale < 1e-12] <- 1
  train_frame <- as.data.frame(base::scale(x_train, center, scale))
  valid_frame <- as.data.frame(base::scale(x_valid, center, scale))
  names(train_frame) <- names(valid_frame) <- make.names(
    names(train_frame), unique = TRUE
  )
  train_frame$A <- factor(training$A, levels = c(0L, 1L))

  previous_threads <- torch::torch_get_num_threads()
  previous_torch_state <- torch::torch_get_rng_state()
  on.exit({
    torch::torch_set_rng_state(previous_torch_state)
    torch::torch_set_num_threads(previous_threads)
  }, add = TRUE)
  torch::torch_set_num_threads(2L)
  torch::torch_manual_seed(as.integer(seed))
  fit <- .mec_with_seed(seed, brulee::brulee_mlp(
    A ~ ., data = train_frame, hidden_units = hidden_units,
    activation = "relu", dropout = 0.1, penalty = 1e-3,
    epochs = as.integer(epochs), learn_rate = 0.01,
    batch_size = min(256L, nrow(train_frame)), stop_iter = 10L,
    verbose = FALSE
  ))
  probabilities <- stats::predict(fit, new_data = valid_frame, type = "prob")
  as.numeric(probabilities[[".pred_1"]])
}

.mec_balanced_tuning_rows <- function(source, tune_n) {
  if (length(source) <= tune_n) return(seq_along(source))
  target <- which(source == 1L)
  controls <- which(source == 0L)
  target_count <- min(length(target), floor(tune_n / 2))
  control_count <- min(length(controls), tune_n - target_count)
  remaining <- tune_n - target_count - control_count
  extra_target <- min(remaining, length(target) - target_count)
  target_count <- target_count + extra_target
  control_count <- control_count + min(
    remaining - extra_target, length(controls) - control_count
  )
  sort(c(utils::head(target, target_count),
         utils::head(controls, control_count)))
}

.mec_tune_bart <- function(training, covariates, seed, tune_n, shrinkage,
                           trim) {
  tuning <- training[.mec_balanced_tuning_rows(training$A, tune_n), ,
                     drop = FALSE]
  if (nrow(tuning) < 20L || length(unique(tuning$A)) < 2L) {
    warning("Too few source observations for BART tuning; using supplied settings.",
            call. = FALSE)
    return(NULL)
  }
  training_rows <- .mec_with_seed(seed, {
    sort(unlist(lapply(sort(unique(tuning$A)), function(source) {
      rows <- which(tuning$A == source)
      if (length(rows) <= 2L) return(rows)
      rows[sample.int(length(rows), max(1L, floor(0.70 * length(rows))))]
    }), use.names = FALSE))
  })
  validation_rows <- setdiff(seq_len(nrow(tuning)), training_rows)
  if (length(unique(tuning$A[training_rows])) < 2L ||
      length(unique(tuning$A[validation_rows])) < 2L) {
    warning("The BART tuning split lacks both source groups; using supplied settings.",
            call. = FALSE)
    return(NULL)
  }
  grid <- data.frame(num_trees = c(25L, 50L, 100L),
                     posterior_draws = 100L, burn_in = 50L)
  loss <- rep(Inf, nrow(grid))
  for (candidate in seq_len(nrow(grid))) {
    loss[candidate] <- tryCatch({
      probability <- .mec_propensity_bart(
        tuning[training_rows, , drop = FALSE],
        tuning[validation_rows, , drop = FALSE], covariates, seed,
        num_trees = grid$num_trees[candidate],
        posterior_draws = grid$posterior_draws[candidate],
        burn_in = grid$burn_in[candidate], shrinkage = shrinkage
      )
      probability <- pmin(pmax(probability, trim[1L]), trim[2L])
      observed <- tuning$A[validation_rows]
      -mean(observed * log(probability) +
              (1 - observed) * log1p(-probability))
    }, error = function(error) Inf)
  }
  if (all(!is.finite(loss))) {
    warning("All BART tuning candidates failed; using supplied settings.",
            call. = FALSE)
    return(NULL)
  }
  selected <- which.min(loss)
  list(num_trees = grid$num_trees[selected],
       posterior_draws = grid$posterior_draws[selected],
       burn_in = grid$burn_in[selected], shrinkage = shrinkage,
       log_loss = loss[selected], candidates = cbind(grid, log_loss = loss),
       tuning_rows = nrow(tuning), training_rows = length(training_rows),
       validation_rows = length(validation_rows), seed = seed)
}

.mec_propensity_bart <- function(training, validation, covariates, seed,
                                 num_trees, posterior_draws, burn_in,
                                 shrinkage, auto_tune = FALSE, tune_n = 100L,
                                 trim = c(0.01, 0.99), return_details = FALSE) {
  if (!requireNamespace("dbarts", quietly = TRUE)) {
    stop("Install {dbarts} to use `ps_learner = 'bart'`.", call. = FALSE)
  }
  .mec_check_factor_levels(training, validation, covariates,
                           "Propensity BART")
  tuning <- if (auto_tune) {
    .mec_tune_bart(training, covariates, seed, tune_n, shrinkage, trim)
  } else NULL
  if (!is.null(tuning)) {
    num_trees <- tuning$num_trees
    posterior_draws <- tuning$posterior_draws
    burn_in <- tuning$burn_in
  }

  design_formula <- stats::reformulate(covariates)
  x_train <- stats::model.matrix(design_formula, training)[, -1L, drop = FALSE]
  x_valid <- stats::model.matrix(design_formula, validation)[, -1L, drop = FALSE]
  if (!identical(colnames(x_train), colnames(x_valid))) {
    stop("Covariate columns changed between BART training and validation.",
         call. = FALSE)
  }
  observed_source <- as.integer(training$A)
  source_fraction <- min(max(mean(observed_source), 0.01), 0.99)
  binary_offset <- stats::qnorm(source_fraction)

  fit <- .mec_with_seed(seed, dbarts::bart(
    x.train = x_train,
    y.train = observed_source,
    x.test = x_valid,
    ntree = as.integer(num_trees),
    ndpost = as.integer(posterior_draws),
    nskip = as.integer(burn_in),
    k = shrinkage,
    binaryOffset = binary_offset,
    verbose = FALSE,
    keeptrees = FALSE,
    keeptrainfits = FALSE,
    nthread = 1L,
    seed = as.integer(seed)
  ))

  if (!is.null(fit$prob.test.mean)) {
    probabilities <- as.numeric(fit$prob.test.mean)
  } else if (!is.null(fit$prob.test)) {
    probabilities <- as.numeric(colMeans(fit$prob.test))
  } else if (!is.null(fit$yhat.test)) {
    probabilities <- as.numeric(colMeans(
      stats::pnorm(fit$yhat.test + binary_offset)
    ))
  } else if (!is.null(fit$yhat.test.mean)) {
    probabilities <- stats::pnorm(
      as.numeric(fit$yhat.test.mean) + binary_offset
    )
  } else {
    stop("BART did not return validation-set source probabilities.",
         call. = FALSE)
  }
  if (return_details) {
    return(list(probabilities = probabilities, tuning = tuning,
                parameters = list(num_trees = num_trees,
                                  posterior_draws = posterior_draws,
                                  burn_in = burn_in, shrinkage = shrinkage,
                                  seed = seed)))
  }
  probabilities
}

.mec_predict_cox_survival <- function(training, validation, covariates,
                                      landmarks) {
  controls <- training[training$A == 0L, , drop = FALSE]
  if (sum(controls$delta) < 2L) {
    stop("A Cox training fold has fewer than two external-control events.",
         call. = FALSE)
  }
  .mec_check_factor_levels(controls, validation, covariates,
                           "Control Cox learner")

  survival_formula <- stats::reformulate(
    covariates, response = "survival::Surv(time, delta)"
  )
  fit <- survival::coxph(survival_formula, data = controls,
                         ties = "breslow", x = TRUE, singular.ok = FALSE,
                         control = survival::coxph.control(timefix = FALSE))
  baseline <- survival::basehaz(fit, centered = FALSE)
  if (!nrow(baseline)) {
    stop("The Cox survival learner has no estimated baseline hazard.",
         call. = FALSE)
  }

  # Before the first observed event the cumulative baseline hazard is zero.
  event_index <- findInterval(landmarks, baseline$time)
  hazard_at_time <- numeric(length(landmarks))
  observed_index <- event_index > 0L
  hazard_at_time[observed_index] <- baseline$hazard[event_index[observed_index]]

  linear_predictor <- as.numeric(stats::predict(
    fit, newdata = validation, type = "lp", reference = "zero"
  ))
  survival <- exp(-outer(exp(linear_predictor), hazard_at_time))
  if (any(!is.finite(survival))) {
    stop("Cox survival predictions are non-finite.", call. = FALSE)
  }
  list(survival = survival, fit = fit, baseline = baseline)
}

.mec_tune_rsf <- function(training, covariates, landmarks, seed,
                           nuisance_settings = "standard", tune_n = 100L) {
  controls <- training[training$A == 0L, , drop = FALSE]
  p <- length(covariates)
  default <- list(mtry = max(1L, floor(sqrt(p))), min_node_size = 15L)
  original_study <- identical(nuisance_settings, "original_study")
  if (original_study) {
    controls <- utils::head(controls, tune_n)
    default$num_trees <- 300L
  }
  if (nrow(controls) < 30L || sum(controls$delta) < 8L) {
    warning("Too few external-control outcomes for RSF tuning; using defaults.",
            call. = FALSE)
    return(default)
  }

  training_rows <- .mec_with_seed(seed, {
    event_rows <- which(controls$delta == 1L)
    censored_rows <- which(controls$delta == 0L)
    if (original_study) {
      draw_half <- function(rows) {
        if (!length(rows)) return(integer())
        rows[sample.int(length(rows), max(1L, floor(0.5 * length(rows))))]
      }
      sort(c(draw_half(event_rows), draw_half(censored_rows)))
    } else {
      sort(c(
        sample(event_rows, floor(0.6 * length(event_rows))),
        sample(censored_rows, floor(0.6 * length(censored_rows)))
      ))
    }
  })
  validation_rows <- setdiff(seq_len(nrow(controls)), training_rows)
  if (length(validation_rows) < 10L ||
      sum(controls$delta[validation_rows]) < 3L) {
    warning("The RSF tuning split has too few validation events; using defaults.",
            call. = FALSE)
    return(default)
  }

  tuning_train <- controls[training_rows, , drop = FALSE]
  tuning_valid <- controls[validation_rows, , drop = FALSE]
  if (original_study) {
    valid_levels <- tryCatch({
      .mec_check_factor_levels(tuning_train, tuning_valid, covariates,
                               "RSF tuning")
      TRUE
    }, error = function(error) FALSE)
    if (!valid_levels) {
      warning("The RSF tuning subset has inadequate factor levels; using defaults.",
              call. = FALSE)
      return(default)
    }
  } else {
    .mec_check_factor_levels(tuning_train, tuning_valid, covariates,
                             "RSF tuning")
  }
  train_frame <- tuning_train[, c("time", "delta", covariates), drop = FALSE]
  validation_frame <- tuning_valid[, covariates, drop = FALSE]
  Surv <- survival::Surv
  survival_formula <- stats::as.formula("Surv(time, delta) ~ .")
  grid <- if (original_study) {
    expand.grid(
      mtry = sort(unique(pmax(1L, pmin(p, c(
        floor(sqrt(p)), ceiling(p / 3), ceiling(p / 2)
      ))))), min_node_size = c(15L, 30L, 50L)
    )
  } else {
    expand.grid(
      mtry = sort(unique(pmax(1L, c(floor(sqrt(p)), ceiling(p / 2))))),
      min_node_size = c(15L, 30L)
    )
  }
  candidate_trees <- if (original_study) 100L else 200L
  discrimination <- rep(-Inf, nrow(grid))
  evaluation_time <- stats::median(landmarks)

  for (candidate in seq_len(nrow(grid))) {
    discrimination[candidate] <- tryCatch({
      fit <- ranger::ranger(
        survival_formula, data = train_frame,
        num.trees = candidate_trees, mtry = grid$mtry[candidate],
        min.node.size = grid$min_node_size[candidate],
        splitrule = if (original_study) "extratrees" else "logrank",
        num.random.splits = 1L, replace = !original_study,
        sample.fraction = if (original_study) 0.632 else 1,
        write.forest = TRUE, oob.error = FALSE,
        num.threads = 1L, seed = as.integer(seed)
      )
      prediction <- stats::predict(fit, data = validation_frame,
                                   num.threads = 1L)
      event_times <- prediction$unique.death.times
      if (is.null(event_times)) {
        event_times <- fit$unique.death.times
      }
      event_index <- findInterval(evaluation_time, event_times)
      predicted_risk <- if (event_index == 0L) {
        rep(0, nrow(tuning_valid))
      } else {
        1 - prediction$survival[, event_index]
      }
      .mec_risk_concordance(
        tuning_valid$time, tuning_valid$delta, predicted_risk
      )
    }, error = function(error) -Inf)
  }
  if (all(!is.finite(discrimination))) {
    if (original_study) {
      warning("All RSF tuning candidates failed; using defaults.", call. = FALSE)
      return(default)
    }
    stop("All RSF tuning candidates failed.", call. = FALSE)
  }
  selected <- which.max(discrimination)
  answer <- list(mtry = as.integer(grid$mtry[selected]),
                 min_node_size = as.integer(grid$min_node_size[selected]),
                 c_index = discrimination[selected],
                 candidate_num_trees = candidate_trees,
                 candidates = cbind(grid, c_index = discrimination),
                 tuning_rows = nrow(controls), training_rows = nrow(tuning_train),
                 validation_rows = nrow(tuning_valid), seed = seed)
  if (original_study) answer$num_trees <- candidate_trees
  answer
}

.mec_risk_concordance <- function(time, event, predicted_risk) {
  survival::concordance(
    survival::Surv(time, event) ~ predicted_risk,
    reverse = TRUE
  )$concordance
}

.mec_predict_rsf_survival <- function(training, validation, covariates,
                                      landmarks, seed, num_trees,
                                      min_node_size, auto_tune,
                                      nuisance_settings = "standard",
                                      tune_n = 100L) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install {ranger} to use `survival_learner = 'rsf'`.", call. = FALSE)
  }
  controls <- training[training$A == 0L, , drop = FALSE]
  if (sum(controls$delta) < 5L) {
    stop("An RSF training fold has fewer than five external-control events.",
         call. = FALSE)
  }
  .mec_check_factor_levels(controls, validation, covariates,
                           "Control RSF learner")

  training_frame <- controls[, c("time", "delta", covariates), drop = FALSE]
  prediction_frame <- validation[, covariates, drop = FALSE]
  Surv <- survival::Surv
  survival_formula <- stats::as.formula("Surv(time, delta) ~ .")
  original_study <- identical(nuisance_settings, "original_study")
  if (auto_tune) {
    tuning <- .mec_tune_rsf(training, covariates, landmarks, seed,
                            nuisance_settings, tune_n)
    mtry <- tuning$mtry
    min_node_size <- tuning$min_node_size
    if (!is.null(tuning$num_trees)) num_trees <- tuning$num_trees
  } else {
    tuning <- NULL
    mtry <- max(1L, floor(sqrt(length(covariates))))
  }
  fit <- ranger::ranger(
    survival_formula, data = training_frame,
    num.trees = as.integer(num_trees),
    mtry = mtry,
    min.node.size = as.integer(min_node_size),
    splitrule = if (original_study) "extratrees" else "logrank",
    num.random.splits = 1L, replace = !original_study,
    sample.fraction = if (original_study) 0.632 else 1,
    write.forest = TRUE, oob.error = FALSE,
    num.threads = 1L, seed = as.integer(seed)
  )
  prediction <- stats::predict(fit, data = prediction_frame, num.threads = 1L)
  event_times <- prediction$unique.death.times
  if (is.null(event_times)) {
    event_times <- fit$unique.death.times
  }
  if (!length(event_times)) {
    stop("The RSF survival learner has no estimated event times.",
         call. = FALSE)
  }

  event_index <- findInterval(landmarks, event_times)
  survival <- matrix(1, nrow = nrow(validation), ncol = length(landmarks))
  observed_index <- which(event_index > 0L)
  if (length(observed_index)) {
    survival[, observed_index] <- prediction$survival[,
      event_index[observed_index], drop = FALSE
    ]
  }
  if (any(!is.finite(survival))) {
    stop("RSF survival predictions are non-finite.", call. = FALSE)
  }
  list(survival = survival, fit = fit, event_times = event_times,
       tuning = tuning,
       parameters = list(num_trees = num_trees, mtry = mtry,
                         min_node_size = min_node_size,
                         splitrule = if (original_study) "extratrees" else "logrank",
                         replace = !original_study,
                         sample_fraction = if (original_study) 0.632 else 1,
                         num_random_splits = 1L, seed = seed))
}

.mec_crossfit_nuisance <- function(data, covariates, ps_learner,
                                   survival_learner, landmarks, fold_id,
                                   ps_trim, seed, mlp_hidden_units, mlp_epochs,
                                   bart_num_trees, bart_posterior_draws,
                                   bart_burn_in, bart_shrinkage,
                                   rsf_num_trees, rsf_min_node_size,
                                   rsf_auto_tune, ps_predictions = NULL,
                                   nuisance_settings = "standard",
                                   bart_auto_tune = FALSE, ml_tune_n = 100L) {
  n <- nrow(data)
  n_landmarks <- length(landmarks)
  ps_oof <- if (is.null(ps_predictions)) rep(NA_real_, n) else ps_predictions
  survival_oof <- matrix(NA_real_, nrow = n, ncol = n_landmarks)
  outcome_fits <- vector("list", max(fold_id))
  propensity_tuning <- vector("list", max(fold_id))

  for (fold in sort(unique(fold_id))) {
    train_rows <- fold_id != fold
    valid_rows <- fold_id == fold
    training <- data[train_rows, , drop = FALSE]
    validation <- data[valid_rows, , drop = FALSE]

    if (is.null(ps_predictions)) {
      if (ps_learner == "glm") {
        ps_oof[valid_rows] <- .mec_propensity_glm(
          training, validation, covariates
        )
      } else if (ps_learner == "mlp") {
        ps_oof[valid_rows] <- .mec_propensity_mlp(
          training, validation, covariates,
          seed = seed + fold,
          hidden_units = mlp_hidden_units, epochs = mlp_epochs
        )
      } else {
        propensity_result <- .mec_propensity_bart(
          training, validation, covariates,
          seed = seed + fold,
          num_trees = bart_num_trees,
          posterior_draws = bart_posterior_draws,
          burn_in = bart_burn_in,
          shrinkage = bart_shrinkage, auto_tune = bart_auto_tune,
          tune_n = ml_tune_n, trim = ps_trim, return_details = TRUE
        )
        ps_oof[valid_rows] <- propensity_result$probabilities
        propensity_tuning[[fold]] <- propensity_result[c("parameters", "tuning")]
      }
    }

    if (survival_learner == "cox") {
      outcome_result <- .mec_predict_cox_survival(
        training, validation, covariates, landmarks
      )
    } else {
      outcome_result <- .mec_predict_rsf_survival(
        training, validation, covariates, landmarks,
        seed = seed + fold,
        num_trees = rsf_num_trees,
        min_node_size = rsf_min_node_size,
        auto_tune = rsf_auto_tune, nuisance_settings = nuisance_settings,
        tune_n = ml_tune_n
      )
    }
    survival_oof[valid_rows, ] <- outcome_result$survival
    outcome_result$survival <- NULL
    outcome_fits[[fold]] <- outcome_result
  }

  if (any(!is.finite(ps_oof)) || any(!is.finite(survival_oof))) {
    stop("Cross-fitting produced missing or non-finite predictions.",
         call. = FALSE)
  }
  if (any(ps_oof < 0 | ps_oof > 1) ||
      any(survival_oof < 0 | survival_oof > 1)) {
    stop("A learner returned predictions outside [0, 1].", call. = FALSE)
  }

  ps_oof <- pmin(pmax(ps_oof, ps_trim[1L]), ps_trim[2L])
  colnames(survival_oof) <- paste0("survival_", seq_len(n_landmarks))
  basis <- cbind(intercept = 1, survival_oof)

  varying <- apply(basis[, -1L, drop = FALSE], 2L, stats::sd) > 1e-10
  retained <- c(TRUE, varying)
  if (!all(varying)) {
    warning("Constant survival-prediction columns were removed from calibration.",
            call. = FALSE)
  }
  basis <- basis[, retained, drop = FALSE]

  list(ps_oof = ps_oof, survival_oof = survival_oof,
       basis = basis, retained_basis = retained,
       outcome_fits = outcome_fits, propensity_tuning = propensity_tuning)
}

.mec_kl_calibrate <- function(control_basis, target_basis, base_weights,
                              tolerance, max_iterations) {
  if (any(!is.finite(base_weights)) || any(base_weights <= 0)) {
    stop("Baseline control weights must be positive and finite.",
         call. = FALSE)
  }
  if (ncol(control_basis) >= nrow(control_basis)) {
    stop("There are too many calibration features for the external controls.",
         call. = FALSE)
  }

  target_total <- colSums(target_basis)
  lambda <- numeric(ncol(control_basis))
  converged <- FALSE

  dual_value <- function(candidate) {
    linear_term <- as.numeric(control_basis %*% candidate)
    weights <- base_weights * exp(linear_term)
    if (any(!is.finite(weights)) || any(weights <= 0)) {
      return(Inf)
    }
    sum(weights) - sum(target_total * candidate)
  }

  for (iteration in seq_len(max_iterations)) {
    weights <- base_weights * exp(as.numeric(control_basis %*% lambda))
    if (any(!is.finite(weights)) || any(weights <= 0)) {
      stop("KL calibration produced non-finite or nonpositive weights.",
           call. = FALSE)
    }
    gradient <- as.numeric(crossprod(control_basis, weights) - target_total)
    maximum_gap <- max(abs(gradient)) / nrow(target_basis)
    if (maximum_gap <= tolerance) {
      converged <- TRUE
      break
    }

    hessian <- crossprod(control_basis * weights, control_basis)
    ridge <- 1e-10 * max(1, max(diag(hessian)))
    direction <- tryCatch(
      solve(hessian + diag(ridge, ncol(hessian)), gradient),
      error = function(error) NULL
    )
    if (is.null(direction) || any(!is.finite(direction))) {
      stop("The KL calibration Hessian is singular.", call. = FALSE)
    }

    old_value <- dual_value(lambda)
    slope <- sum(gradient * direction)
    step_size <- 1
    accepted <- FALSE
    for (attempt in seq_len(40L)) {
      candidate <- lambda - step_size * as.numeric(direction)
      candidate_value <- dual_value(candidate)
      if (is.finite(candidate_value) &&
          candidate_value <= old_value - 1e-4 * step_size * slope) {
        lambda <- candidate
        accepted <- TRUE
        break
      }
      step_size <- step_size / 2
    }
    if (!accepted) {
      stop("KL calibration could not find a decreasing Newton step.",
           call. = FALSE)
    }
  }

  weights <- base_weights * exp(as.numeric(control_basis %*% lambda))
  balance_residual <- as.numeric(crossprod(control_basis, weights) -
                                   target_total)
  if (!converged) {
    stop(sprintf(
      "KL calibration did not converge; maximum mean balance residual is %.4g.",
      max(abs(balance_residual)) / nrow(target_basis)
    ), call. = FALSE)
  }

  list(
    weights = as.numeric(weights), lambda = lambda,
    balance_residual = balance_residual,
    max_mean_balance_residual = max(abs(balance_residual)) /
      nrow(target_basis),
    iterations = iteration, converged = TRUE,
    hessian = crossprod(control_basis * weights, control_basis)
  )
}

.mec_cox_score <- function(theta, data, weights, individual = TRUE) {
  source <- data$A
  time <- data$time
  event <- data$delta
  n <- nrow(data)
  relative_risk <- exp(theta * source)

  descending <- order(time, decreasing = TRUE)
  time_descending <- time[descending]
  weighted_risk <- weights[descending] * relative_risk[descending]
  weighted_source_risk <- weighted_risk * source[descending]
  end_of_time_group <- cumsum(rle(time_descending)$lengths)
  distinct_times <- rle(time_descending)$values
  risk_sum <- cumsum(weighted_risk)[end_of_time_group]
  source_risk_sum <- cumsum(weighted_source_risk)[end_of_time_group]

  group_index <- match(time, distinct_times)
  source_fraction <- source_risk_sum[group_index] / risk_sum[group_index]
  event_contribution <- weights * event * (source - source_fraction)
  score <- sum(event_contribution)
  if (!individual) return(list(score = score))

  event_times <- sort(unique(time[event == 1L]))
  event_weight <- as.numeric(tapply(
    weights * event, factor(time, levels = event_times), sum
  ))
  risk_index <- match(event_times, distinct_times)
  event_risk <- risk_sum[risk_index]
  event_source_fraction <- source_risk_sum[risk_index] / event_risk
  cum_risk_term <- cumsum(event_weight / event_risk)
  cum_source_term <- cumsum(
    event_weight * event_source_fraction / event_risk
  )

  last_event <- findInterval(time, event_times)
  risk_term <- source_term <- numeric(n)
  has_prior_event <- last_event > 0L
  risk_term[has_prior_event] <- cum_risk_term[last_event[has_prior_event]]
  source_term[has_prior_event] <- cum_source_term[last_event[has_prior_event]]
  risk_contribution <- weights * relative_risk *
    (source * risk_term - source_term)

  list(score = score,
       contribution = event_contribution - risk_contribution)
}

.mec_central_difference <- function(fun, x, step = 1e-5) {
  result <- numeric(length(x))
  for (column in seq_along(x)) {
    increment <- step * max(1, abs(x[column]))
    positive <- negative <- x
    positive[column] <- positive[column] + increment
    negative[column] <- negative[column] - increment
    result[column] <- (fun(positive) - fun(negative)) / (2 * increment)
  }
  result
}

.mec_stacked_variance <- function(theta, data, weights, base_weights,
                                  control_basis, target_basis,
                                  control_rows, target_rows, lambda,
                                  calibration_hessian) {
  n <- nrow(data)
  score_result <- .mec_cox_score(theta, data, weights)
  derivative_theta <- .mec_central_difference(
    function(candidate) .mec_cox_score(
      candidate[1L], data, weights, individual = FALSE
    )$score,
    x = theta
  )

  score_at_lambda <- function(candidate) {
    control_weights <- base_weights * exp(
      as.numeric(control_basis %*% candidate)
    )
    candidate_weights <- numeric(n)
    candidate_weights[target_rows] <- 1
    candidate_weights[control_rows] <- control_weights
    .mec_cox_score(theta, data, candidate_weights, individual = FALSE)$score
  }
  derivative_lambda <- .mec_central_difference(score_at_lambda, lambda)

  n_basis <- ncol(control_basis)
  calibration_contribution <- matrix(0, nrow = n, ncol = n_basis)
  calibration_contribution[control_rows, ] <-
    control_basis * weights[control_rows]
  calibration_contribution[target_rows, ] <- -target_basis
  subject_contribution <- cbind(
    score_result$contribution, calibration_contribution
  )

  jacobian <- rbind(
    c(derivative_theta, derivative_lambda),
    cbind(rep(0, n_basis), calibration_hessian)
  )
  if (any(!is.finite(jacobian))) {
    stop("The stacked-sandwich Jacobian contains non-finite values.",
         call. = FALSE)
  }
  decomposition <- svd(jacobian)
  singular_values <- decomposition$d
  retained <- singular_values > max(singular_values) * 1e-8
  if (!any(retained)) {
    stop("The stacked-sandwich Jacobian has no estimable directions.",
         call. = FALSE)
  }
  used_pseudoinverse <- !all(retained)
  if (used_pseudoinverse) {
    warning(
      "The stacked-sandwich Jacobian is rank deficient; using its generalized inverse.",
      call. = FALSE
    )
  }
  inverse_jacobian <- decomposition$v[, retained, drop = FALSE] %*%
    (t(decomposition$u[, retained, drop = FALSE]) /
       singular_values[retained])
  influence <- -subject_contribution %*% t(inverse_jacobian)
  standard_error <- sqrt(sum(influence[, 1L]^2))
  if (!is.finite(standard_error) || standard_error <= 0) {
    stop("The stacked-sandwich standard error is invalid.", call. = FALSE)
  }

  list(standard_error = standard_error,
       influence = influence[, 1L],
       jacobian_condition = max(singular_values) / min(singular_values),
       jacobian_rank = sum(retained),
       used_pseudoinverse = used_pseudoinverse)
}

.mec_weight_diagnostics <- function(weights) {
  total <- sum(weights)
  list(
    ess = total^2 / sum(weights^2),
    cv = stats::sd(weights) / mean(weights),
    minimum = min(weights), maximum = max(weights),
    sum = total
  )
}

#' Fit MEC-Cox with cross-fitted prognostic calibration
#'
#' External-control source-propensity odds are normalized to the target sample
#' size, then adjusted by Kullback-Leibler calibration. The calibrated weights
#' balance an intercept and cross-fitted predictions of control survival at
#' selected landmarks. Both the source-propensity and survival learners are
#' cross-fitted by default and use the same source-stratified folds.
#'
#' @param data A data frame containing one row per person.
#' @param time,event,source Names of the follow-up time, event indicator, and
#'   source indicator columns. `source = 1` identifies the target cohort;
#'   `source = 0` identifies external controls.
#' @param covariates Names of baseline covariates used by both learners.
#' @param ps_learner Source-propensity learner: logistic GLM, MLP, or BART.
#' @param survival_learner Control-survival learner: Cox or random survival
#'   forest. The survival learner is trained on external controls only.
#' @param n_folds Number of source-stratified cross-fitting folds.
#' @param n_landmarks Number of event-time quantile landmarks when `landmarks`
#'   is not supplied.
#' @param landmarks Optional numeric vector of times in the units of `time`.
#' @param seed Seed for fold assignment and stochastic learners.
#' @param conf_level Confidence level for a working Wald interval.
#' @param ps_trim Lower and upper limits for cross-fitted source propensity.
#' @param calibration_tol Maximum allowed balance discrepancy per target person.
#' @param max_calibration_iter Maximum number of KL Newton steps.
#' @param mlp_hidden_units,mlp_epochs MLP architecture and training epochs.
#' @param bart_num_trees,bart_posterior_draws,bart_burn_in,bart_shrinkage
#'   BART training settings, used directly when tuning is disabled or as
#'   fallbacks if tuning fails. The standard defaults are 100 trees, 1000
#'   posterior draws, 500 burn-in draws, and shrinkage 2. Under
#'   `nuisance_settings = "original_study"`, unspecified values instead use
#'   50 trees, 200 posterior draws, and 100 burn-in draws, with shrinkage 2.
#' @param rsf_num_trees,rsf_min_node_size Random survival forest controls.
#' @param rsf_auto_tune If `TRUE`, select the RSF `mtry` and minimum node size
#'   within each training fold using a source-control training/validation split.
#' @param ps_predictions Optional numeric vector of out-of-fold source
#'   probabilities, one per row of `data`, strictly inside (0, 1). Supplying this
#'   vector skips propensity fitting and applies `ps_trim` as usual. The caller
#'   must ensure the predictions use the same rows, row order, covariates, and
#'   source-stratified folds defined by `n_folds` and `seed`. To reuse a previous
#'   fit's `ps_oof` when changing the survival learner, keep these inputs and
#'   `ps_trim` unchanged; earlier clipping cannot be undone. This function checks
#'   the vector's length and values, not its training provenance.
#' @param nuisance_settings Learner configuration. `"standard"` uses the
#'   general-purpose settings. `"original_study"` uses the lightweight tuning
#'   and randomized survival-forest preset:
#'   source-balanced BART tuning on at most `ml_tune_n` rows with a 70/30
#'   stratified split; and RSF tuning on the first at most `ml_tune_n` controls
#'   with a 50/50 event-stratified split. The RSF uses extremely randomized
#'   splits, one candidate split, sampling without replacement at fraction
#'   0.632, and 100 trees after successful tuning (300 trees by default without
#'   tuning or if tuning fails). This API uses fold-specific seeds, averaging
#'   of transformed BART draws, and directional mortality concordance.
#'   The main-paper simulation scripts use their bundled study functions.
#' @param bart_auto_tune Whether to choose BART tree count from 25, 50, and 100
#'   using held-out log loss within each training fold. Candidates and the final
#'   selected fit use 100 posterior draws and 50 burn-in draws; shrinkage stays
#'   at `bart_shrinkage`. `NULL` enables tuning for `"original_study"` and
#'   disables it for `"standard"`. Set `FALSE` to use the supplied BART controls.
#' @param ml_tune_n Maximum subset size for BART tuning and for RSF tuning under
#'   `"original_study"`. Standard RSF tuning continues to use all training-fold
#'   controls. Tuning subsets never include the held-out cross-fitting fold.
#' @param ... Reserved; unsupported arguments cause an error.
#'
#' @return An object of class `mec_cox_fit` containing the log-hazard-ratio
#'   estimate, working stacked-sandwich standard error, final and baseline
#'   weights, out-of-fold predictions, calibration basis, fold assignments,
#'   model fits, and calibration and weight diagnostics. The working sandwich
#'   accounts for estimation of the calibration multiplier but treats the
#'   cross-fitted nuisance learners as fixed.
#'   `propensity_tuning` and `outcome_fits` record the fold-specific learner
#'   parameters and tuning results. Supplied propensity predictions have no
#'   recorded propensity-training metadata.
#' @examples
#' data(example_external_controls)
#' fit <- fit_mec_cox(
#'   example_external_controls,
#'   time = "time", event = "event", source = "source",
#'   covariates = c("age", "sex", "marker"),
#'   ps_learner = "glm", survival_learner = "cox",
#'   n_folds = 3, n_landmarks = 3, seed = 12
#' )
#' fit$hr
#' @export
fit_mec_cox <- function(data, time, event, source, covariates,
                        ps_learner = c("glm", "mlp", "bart"),
                        survival_learner = c("cox", "rsf"),
                        n_folds = 10L, n_landmarks = 20L,
                        landmarks = NULL, seed = 20260427L,
                        conf_level = 0.95,
                        ps_trim = c(0.01, 0.99),
                        calibration_tol = 1e-8,
                        max_calibration_iter = 100L,
                        mlp_hidden_units = c(32L, 16L),
                        mlp_epochs = 200L,
                        bart_num_trees = 100L,
                        bart_posterior_draws = 1000L,
                        bart_burn_in = 500L,
                        bart_shrinkage = 2,
                        rsf_num_trees = 500L,
                        rsf_min_node_size = 15L,
                        rsf_auto_tune = TRUE, ps_predictions = NULL,
                        nuisance_settings = c("standard", "original_study"),
                        bart_auto_tune = NULL, ml_tune_n = 100L, ...) {
  call <- match.call()
  if (length(list(...))) {
    stop("Unsupported argument in `...`.", call. = FALSE)
  }
  ps_learner <- match.arg(ps_learner)
  survival_learner <- match.arg(survival_learner)
  nuisance_settings <- match.arg(nuisance_settings)
  if (is.null(bart_auto_tune)) {
    bart_auto_tune <- identical(nuisance_settings, "original_study")
  }
  if (!is.logical(bart_auto_tune) || length(bart_auto_tune) != 1L ||
      is.na(bart_auto_tune)) {
    stop("`bart_auto_tune` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }
  if (!is.numeric(ml_tune_n) || length(ml_tune_n) != 1L ||
      !is.finite(ml_tune_n) || ml_tune_n < 1L ||
      ml_tune_n != as.integer(ml_tune_n)) {
    stop("`ml_tune_n` must be a positive integer.", call. = FALSE)
  }
  if (nuisance_settings == "original_study") {
    if (missing(bart_num_trees)) bart_num_trees <- 50L
    if (missing(bart_posterior_draws)) bart_posterior_draws <- 200L
    if (missing(bart_burn_in)) bart_burn_in <- 100L
    if (missing(rsf_num_trees)) rsf_num_trees <- 300L
  }
  if (!is.character(covariates) || !length(covariates)) {
    stop("`covariates` must name at least one baseline covariate.",
         call. = FALSE)
  }
  analysis_data <- .pkg_validate_data(
    data, time, event, source, covariates
  )
  if (!is.null(ps_predictions)) {
    if (!is.numeric(ps_predictions) || !is.null(dim(ps_predictions)) ||
        length(ps_predictions) != nrow(analysis_data) ||
        any(!is.finite(ps_predictions)) ||
        any(ps_predictions <= 0 | ps_predictions >= 1)) {
      stop("`ps_predictions` must be a finite numeric vector with one ",
           "probability strictly inside (0, 1) per row of `data`.",
           call. = FALSE)
    }
    ps_predictions <- as.numeric(ps_predictions)
  }

  if (!is.numeric(n_folds) || length(n_folds) != 1L ||
      !is.finite(n_folds) ||
      n_folds < 2L || n_folds != as.integer(n_folds)) {
    stop("`n_folds` must be an integer of at least two.", call. = FALSE)
  }
  group_size <- table(factor(analysis_data$A, levels = c(0L, 1L)))
  if (any(group_size < n_folds)) {
    stop("Each source must contain at least `n_folds` patients.",
         call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("`seed` must be a finite integer.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be between zero and one.", call. = FALSE)
  }
  if (!is.numeric(ps_trim) || length(ps_trim) != 2L ||
      any(!is.finite(ps_trim)) ||
      ps_trim[1L] <= 0 || ps_trim[1L] >= ps_trim[2L] ||
      ps_trim[2L] >= 1) {
    stop("`ps_trim` must contain two increasing values inside (0, 1).",
         call. = FALSE)
  }
  if (!is.numeric(calibration_tol) || length(calibration_tol) != 1L ||
      !is.finite(calibration_tol) ||
      calibration_tol <= 0) {
    stop("`calibration_tol` must be positive and finite.", call. = FALSE)
  }
  if (!is.numeric(max_calibration_iter) ||
      length(max_calibration_iter) != 1L ||
      !is.finite(max_calibration_iter) ||
      max_calibration_iter < 1L ||
      max_calibration_iter != as.integer(max_calibration_iter)) {
    stop("`max_calibration_iter` must be a positive integer.",
         call. = FALSE)
  }
  if (!is.logical(rsf_auto_tune) || length(rsf_auto_tune) != 1L ||
      is.na(rsf_auto_tune)) {
    stop("`rsf_auto_tune` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(rsf_num_trees) || length(rsf_num_trees) != 1L ||
      !is.finite(rsf_num_trees) || rsf_num_trees < 1L ||
      rsf_num_trees != as.integer(rsf_num_trees)) {
    stop("`rsf_num_trees` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(rsf_min_node_size) || length(rsf_min_node_size) != 1L ||
      !is.finite(rsf_min_node_size) || rsf_min_node_size < 1L ||
      rsf_min_node_size != as.integer(rsf_min_node_size)) {
    stop("`rsf_min_node_size` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(mlp_hidden_units) || !length(mlp_hidden_units) ||
      any(!is.finite(mlp_hidden_units)) || any(mlp_hidden_units < 1L) ||
      any(mlp_hidden_units != as.integer(mlp_hidden_units))) {
    stop("`mlp_hidden_units` must contain positive integers.",
         call. = FALSE)
  }
  if (!is.numeric(mlp_epochs) || length(mlp_epochs) != 1L ||
      !is.finite(mlp_epochs) || mlp_epochs < 1L ||
      mlp_epochs != as.integer(mlp_epochs)) {
    stop("`mlp_epochs` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(bart_num_trees) || length(bart_num_trees) != 1L ||
      !is.finite(bart_num_trees) || bart_num_trees < 1L ||
      bart_num_trees != as.integer(bart_num_trees)) {
    stop("`bart_num_trees` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(bart_posterior_draws) ||
      length(bart_posterior_draws) != 1L ||
      !is.finite(bart_posterior_draws) || bart_posterior_draws < 1L ||
      bart_posterior_draws != as.integer(bart_posterior_draws)) {
    stop("`bart_posterior_draws` must be a positive integer.",
         call. = FALSE)
  }
  if (!is.numeric(bart_burn_in) || length(bart_burn_in) != 1L ||
      !is.finite(bart_burn_in) || bart_burn_in < 1L ||
      bart_burn_in != as.integer(bart_burn_in)) {
    stop("`bart_burn_in` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(bart_shrinkage) || length(bart_shrinkage) != 1L ||
      !is.finite(bart_shrinkage) || bart_shrinkage <= 0) {
    stop("`bart_shrinkage` must be a positive number.", call. = FALSE)
  }

  landmark_times <- .mec_landmarks(
    analysis_data, n_landmarks, landmarks
  )
  fold_id <- .mec_make_folds(analysis_data$A, n_folds, seed)
  nuisance <- .mec_crossfit_nuisance(
    analysis_data, covariates, ps_learner, survival_learner,
    landmark_times, fold_id, ps_trim, seed,
    mlp_hidden_units, mlp_epochs,
    bart_num_trees, bart_posterior_draws,
    bart_burn_in, bart_shrinkage,
    rsf_num_trees, rsf_min_node_size,
    rsf_auto_tune, ps_predictions, nuisance_settings, bart_auto_tune, ml_tune_n
  )

  target_rows <- which(analysis_data$A == 1L)
  control_rows <- which(analysis_data$A == 0L)
  source_odds <- nuisance$ps_oof / (1 - nuisance$ps_oof)
  baseline_control <- source_odds[control_rows]
  baseline_control <- length(target_rows) * baseline_control /
    sum(baseline_control)

  control_basis <- nuisance$basis[control_rows, , drop = FALSE]
  target_basis <- nuisance$basis[target_rows, , drop = FALSE]
  calibration <- .mec_kl_calibrate(
    control_basis, target_basis, baseline_control,
    tolerance = calibration_tol,
    max_iterations = as.integer(max_calibration_iter)
  )

  baseline_weights <- weights <- numeric(nrow(analysis_data))
  baseline_weights[target_rows] <- weights[target_rows] <- 1
  baseline_weights[control_rows] <- baseline_control
  weights[control_rows] <- calibration$weights

  cox_fit <- survival::coxph(
    survival::Surv(time, delta) ~ A,
    data = analysis_data, weights = weights,
    robust = FALSE, ties = "breslow",
    control = survival::coxph.control(timefix = FALSE)
  )
  theta <- unname(stats::coef(cox_fit)["A"])
  if (!is.finite(theta)) {
    stop("The final weighted Cox coefficient is not finite.", call. = FALSE)
  }
  variance <- .mec_stacked_variance(
    theta, analysis_data, weights, baseline_control,
    control_basis, target_basis, control_rows, target_rows,
    calibration$lambda, calibration$hessian
  )
  standard_error <- variance$standard_error
  critical_value <- stats::qnorm((1 + conf_level) / 2)
  log_interval <- theta + c(-1, 1) * critical_value * standard_error

  calibration$weights <- NULL
  calibration$hessian <- NULL
  fit <- list(
    call = call,
    theta = theta,
    se = standard_error,
    hr = exp(theta),
    conf_int = stats::setNames(exp(log_interval), c("lower", "upper")),
    conf_level = conf_level,
    variance_method = "stacked_working",
    weights = weights,
    baseline_weights = baseline_weights,
    base_weights = baseline_weights,
    ps_oof = nuisance$ps_oof,
    survival_oof = nuisance$survival_oof,
    basis = nuisance$basis,
    landmark_times = landmark_times,
    retained_basis = nuisance$retained_basis,
    fold_id = fold_id,
    outcome_fits = nuisance$outcome_fits,
    propensity_tuning = nuisance$propensity_tuning,
    analysis_data = analysis_data,
    covariates = covariates,
    ps_learner = ps_learner,
    ps_predictions_supplied = !is.null(ps_predictions),
    nuisance_settings = nuisance_settings,
    bart_auto_tune = bart_auto_tune,
    ml_tune_n = ml_tune_n,
    survival_learner = survival_learner,
    rsf_auto_tune = rsf_auto_tune,
    calibration = calibration,
    weight_diagnostics = .mec_weight_diagnostics(weights[control_rows]),
    variance_diagnostics = list(
      jacobian_condition = variance$jacobian_condition,
      jacobian_rank = variance$jacobian_rank,
      used_pseudoinverse = variance$used_pseudoinverse
    ),
    cox_fit = cox_fit
  )
  class(fit) <- "mec_cox_fit"
  fit
}

# Reuse the fold-specific outcome fits when diagnostics request times that were
# not part of the calibration basis. The returned rows follow analysis_data.
.mec_oof_survival_at_times <- function(object, times) {
  if (!inherits(object, "mec_cox_fit")) {
    stop("`object` must be a MEC-Cox fit.", call. = FALSE)
  }
  if (!is.numeric(times) || !length(times) ||
      any(!is.finite(times)) || any(times <= 0)) {
    stop("`times` must contain positive, finite values.", call. = FALSE)
  }

  survival <- matrix(NA_real_, nrow = nrow(object$analysis_data),
                     ncol = length(times))
  for (fold in sort(unique(object$fold_id))) {
    rows <- object$fold_id == fold
    fold_fit <- object$outcome_fits[[fold]]
    validation <- object$analysis_data[rows, , drop = FALSE]
    if (object$survival_learner == "cox") {
      baseline <- fold_fit$baseline
      event_index <- findInterval(times, baseline$time)
      hazard <- numeric(length(times))
      observed <- event_index > 0L
      hazard[observed] <- baseline$hazard[event_index[observed]]
      linear_predictor <- as.numeric(stats::predict(
        fold_fit$fit, newdata = validation,
        type = "lp", reference = "zero"
      ))
      survival[rows, ] <- exp(-outer(exp(linear_predictor), hazard))
    } else {
      prediction <- stats::predict(
        fold_fit$fit,
        data = validation[, object$covariates, drop = FALSE],
        num.threads = 1L
      )
      event_times <- fold_fit$event_times
      event_index <- findInterval(times, event_times)
      fold_survival <- matrix(1, nrow = nrow(validation),
                              ncol = length(times))
      observed <- which(event_index > 0L)
      if (length(observed)) {
        fold_survival[, observed] <- prediction$survival[,
          event_index[observed], drop = FALSE
        ]
      }
      survival[rows, ] <- fold_survival
    }
  }
  colnames(survival) <- paste0("time_", times)
  survival
}
