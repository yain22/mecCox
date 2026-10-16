# Driver orchestration only. The verified final-study adapter owns model fitting.
load_final_driver_reference <- function(directory) {
  api <- new.env(parent = globalenv())
  source(file.path(directory, "final_study_helpers.R"), local = api)
  list(api = api, bundle = api$load_final_study_bundle(
    file.path(directory, "final_reference")))
}

check_final_driver_replications <- function(replications) {
  if (length(replications) != 1L || !is.numeric(replications) ||
      !is.finite(replications) || replications < 1L || replications > 1000L ||
      replications != as.integer(replications)) {
    stop("The frozen final study supports 1 to 1000 replications per cell.",
         call. = FALSE)
  }
  invisible(replications)
}

final_driver_target <- function(reference, figure, setting) {
  targets <- reference$bundle$targets
  target <- targets[tolower(targets$figure) == figure &
                      targets$setting == setting, , drop = FALSE]
  if (nrow(target) != 1L || !is.finite(target$numerical_target)) {
    stop("The requested final-study target is not uniquely defined.")
  }
  target
}

final_driver_job <- function(reference, figure, setting, n1, replicate) {
  grid <- reference$bundle$grid
  job <- grid[grid$figure == figure & grid$setting == setting &
                grid$n1 == n1 & grid$rep_id == replicate, , drop = FALSE]
  if (nrow(job) != 1L) stop("The requested run is not in the frozen final grid.")
  job
}

final_driver_methods <- function(raw, figure) {
  method <- rep(NA_character_, nrow(raw))
  method[raw$Method == "ATT-IPW Cox: Naive model-based"] <- "Naive"
  method[raw$Method == "ATT-IPW Cox: Lin-Wei"] <- "Robust sandwich"
  method[raw$Method == "ATT-IPW Cox: Shu"] <- "Corrected sandwich"
  cox <- grepl("^MEC-Cox:", raw$Method) & grepl("OR=cox", raw$Method, fixed = TRUE)
  rsf <- grepl("^MEC-Cox:", raw$Method) & grepl("OR=rsf", raw$Method, fixed = TRUE)
  method[cox] <- if (figure == "figure5") "MEC-Cox" else "MEC-Cox (BART/Cox)"
  if (figure == "figure6") method[rsf] <- "MEC-Cox (BART/RSF)"
  expected <- c("Naive", "Robust sandwich", "Corrected sandwich",
                if (figure == "figure5") "MEC-Cox" else
                  c("MEC-Cox (BART/Cox)", "MEC-Cox (BART/RSF)"))
  if (anyNA(method) || anyDuplicated(method) || !setequal(method, expected)) {
    stop("The final-study engine returned an unexpected set of method rows.")
  }
  method
}

final_driver_run <- function(reference, figure, setting, n1, replicate) {
  job <- final_driver_job(reference, figure, setting, n1, replicate)
  run <- reference$api$fit_final_study_job(job$id, reference$bundle,
                                          capture_full = FALSE)
  raw <- reference$api$bind_final_study_target(run, reference$bundle)
  list(job = job, raw = raw, run = run,
       methods = final_driver_methods(raw, figure))
}

final_driver_add_diagnostics <- function(answer, fitted) {
  raw <- fitted$raw
  answer$ci_lower <- raw$CI_L
  answer$ci_upper <- raw$CI_U
  answer$original_method <- raw$Method
  answer$final_study_job_id <- fitted$job$id
  # Retain the existing column as an alias of the actual frozen global job ID.
  answer$original_rng_job <- fitted$job$id
  for (name in c("ESS", "Rel_ESS", "W_CV", "W_Min", "W_Max", "Cal_Grad",
                 "Cal_Converged", "target_MCSE", "legacy_target")) {
    answer[[name]] <- raw[[name]]
  }
  answer$raw_sandwich_SE <- NA_real_
  answer$SE_fallback_used <- NA
  answer$KL_positive_finite_weights <- NA
  answer$KL_exponential_weights_identical <- NA
  answer$actual_solve_path <- NA_character_
  answer$fatal_error <- fitted$run$fatal_error
  for (trace in fitted$run$estimator_trace) {
    index <- grepl(paste0("OR=", trace$or_method), raw$Method, fixed = TRUE)
    answer$raw_sandwich_SE[index] <- trace$raw_sandwich_SE
    answer$SE_fallback_used[index] <- trace$fallback_used
  }
  for (trace in fitted$run$math_trace) {
    index <- grepl(paste0("OR=", trace$or_method), raw$Method, fixed = TRUE)
    answer$KL_positive_finite_weights[index] <- trace$positive_finite_weights
    answer$KL_exponential_weights_identical[index] <- trace$exponential_weights_identical
    answer$actual_solve_path[index] <- trace$solve_path
  }
  # The controller retains complete capture_full=FALSE records separately from
  # the ordinary table, including raw target-NA rows, errors and diagnostic traces.
  attr(answer, "final_diagnostics") <- fitted$run
  answer
}

final_driver_coverage <- function(successful) {
  valid <- is.finite(successful$ci_lower) & is.finite(successful$ci_upper) &
    successful$ci_lower <= successful$ci_upper
  covered <- valid & successful$ci_lower <= successful$target &
    successful$ci_upper >= successful$target
  list(numerator = sum(covered), denominator = sum(valid),
       coverage = if (any(valid)) sum(covered) / sum(valid) else NA_real_)
}

final_driver_export_names <- function() {
  c("final_driver_target", "final_driver_job", "final_driver_methods",
    "final_driver_run", "final_driver_add_diagnostics", "final_driver_coverage")
}
