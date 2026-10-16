study_hash <- function(x) digest::digest(x, algo = "sha256")

study_rng <- function() get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

study_md5 <- function(p) unname(tools::md5sum(p))

study_sources <- function(m) {
    stopifnot(identical(vapply(m$source_files, study_md5, character(1)), m$source_md5))
    invisible(TRUE)
}

study_load_engine <- function(directory, family, derivative = "analytic") {
    stopifnot(family %in% c("S2", "S34"), derivative %in% c("analytic", "fd", "original"))
    name <- if (derivative == "original") 
        paste0(family, "_original.R")
    else paste0(family, "_engine_", derivative, ".R")
    e <- new.env(parent = globalenv())
    for (x in parse(file.path(directory, name), keep.source = FALSE)) if (is.call(x) && identical(x[[1]], as.name("<-")) && is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function"))) 
        eval(x, e)
    e
}

study_make_engine <- function(directory, family, derivative = "analytic", logging = TRUE, capture_full = FALSE) {
    e <- study_load_engine(directory, family, derivative)
    e$.data_hash <- NULL
    e$.PS_trace <- e$.nuisance_trace <- e$.calibration_trace <- e$.estimator_trace <- e$.math_trace <- list()
    e$.current_generator <- NA_character_
    if (!logging) 
        return(e)
    e$generate_observed_data <- local({
        raw <- e$generate_observed_data
        eng <- e
        function(...) {
            ans <- raw(...)
            eng$.data_hash <- study_hash(ans)
            ans
        }
    })
    e$estimate_ps_ml <- local({
        raw <- e$estimate_ps_ml
        eng <- e
        function(...) {
            args <- list(...)
            ans <- do.call(raw, args)
            pred <- if (is.null(args$newdata)) 
                args$data
            else args$newdata
            p <- as.numeric(ans$ehat)
            y <- as.integer(pred$A)
            eng$.PS_trace[[length(eng$.PS_trace) + 1L]] <- list(generator = eng$.current_generator, ps_method = args$ps_method, tune_seed = args$tune_seed, n_train = nrow(args$data), train_ids_hash = study_hash(args$data$id), ids = pred$id, y = y, p = p, brier = mean((p - y)^2), logloss = eng$binary_logloss(y, p), tuned_params = ans$tuned_params, R_rng_after = study_rng())
            ans
        }
    })
    e$estimate_mec_nuisance <- local({
        raw <- e$estimate_mec_nuisance
        eng <- e
        full <- capture_full
        function(...) {
            args <- list(...)
            ans <- do.call(raw, args)
            eng$.nuisance_trace[[length(eng$.nuisance_trace) + 1L]] <- list(generator = eng$.current_generator, data_hash = study_hash(args$data), ehat = ans$ehat, H_hash = study_hash(ans$H), H = if (full) ans$H else NULL, fold_id = ans$fold_id, landmark_times = ans$landmark_times, basis_method = ans$basis_method, basis_dim = ans$basis_dim, lambda_selected = ans$lambda_selected, n_selected = ans$n_selected, R_rng_after = study_rng())
            ans
        }
    })
    e$calibrate_weights_bregman_from_design_return <- local({
        raw <- e$calibrate_weights_bregman_from_design_return
        eng <- e
        function(...) {
            args <- list(...)
            ans <- tryCatch(do.call(raw, args), error = function(ex) ex)
            if (inherits(ans, "error")) {
                eng$.calibration_trace[[length(eng$.calibration_trace) + 1L]] <- list(generator = eng$.current_generator, error = conditionMessage(ans))
                stop(ans)
            }
            eng$.calibration_trace[[length(eng$.calibration_trace) + 1L]] <- list(generator = eng$.current_generator, Xs_hash = study_hash(args$Xs), Xp_hash = study_hash(args$Xp), d_hash = study_hash(args$d), weights_hash = study_hash(ans$w), lambda = ans$lambda, C = ans$D_lambda_lambda, gradient_norm = ans$grad_norm, converged = ans$converged, finite_weights = all(is.finite(ans$w)), weight_min = min(ans$w), weight_max = max(ans$w), weight_sum = sum(ans$w), error = NA_character_)
            ans
        }
    })
    fit_name <- if (family == "S2") 
        "fit_mec_cox_from_nuisance"
    else "fit_mec_cox"
    f <- e[[fit_name]]
    b <- as.list(body(f))
    pos <- which(vapply(b, function(x) identical(x, quote(se_mec <- sqrt(sum(phi_theta^2)))), logical(1)))
    stopifnot(length(pos) == 1L)
    e$.capture_full <- capture_full
    inserted <- quote(.math_trace[[length(.math_trace) + 1L]] <<- list(generator = divergence, a = D_theta_theta, b = D_theta_lambda, C = D_lambda_lambda, lambda = cal$lambda, raw_sandwich_SE = se_mec, positive_finite_weights = all(is.finite(w0_mec) & w0_mec > 0), H_hash = study_hash(H), W_hash = study_hash(W_mec), eta_hash = study_hash(eta_theta), full = if (.capture_full) list(H0 = H0, eta_controls = eta_theta[idx0], D = D, IF = IF, rho = rho) else NULL))
    body(f) <- as.call(append(b, list(inserted), after = pos))
    environment(f) <- e
    e[[fit_name]] <- local({
        raw <- f
        eng <- e
        full <- capture_full
        function(...) {
            args <- list(...)
            eng$.current_generator <- args$divergence
            before_count <- length(eng$.math_trace)
            ans <- tryCatch(do.call(raw, args), error = function(ex) ex)
            if (inherits(ans, "error")) {
                eng$.estimator_trace[[length(eng$.estimator_trace) + 1L]] <- list(generator = eng$.current_generator, error = conditionMessage(ans))
                stop(ans)
            }
            point <- ans
            point$se <- NULL
            raw_se <- if (length(eng$.math_trace) > before_count) 
                tail(eng$.math_trace, 1)[[1]]$raw_sandwich_SE
            else NA_real_
            eng$.estimator_trace[[length(eng$.estimator_trace) + 1L]] <- list(generator = eng$.current_generator, point_hash = study_hash(point), point = if (full) point else NULL, theta = ans$theta, W_hash = study_hash(ans$W), H_hash = study_hash(ans$H), reported_SE = ans$se, raw_sandwich_SE = raw_se, fallback_used = !is.finite(raw_se), error = NA_character_)
            ans
        }
    })
    e
}

study_one <- function(id, m, directory, derivative = "analytic", logging = TRUE, capture_full = FALSE) {
    job <- m$grid[m$grid$id == id, , drop = FALSE]
    stopifnot(nrow(job) == 1L)
    cfg <- m$configs[[job$config_id]]
    e <- study_make_engine(directory, job$family, derivative, logging, capture_full)
    aa <- cfg[intersect(names(cfg), names(formals(e$fit_one_rep)))]
    aa$rep_id <- job$rep_id
    aa$n1 <- job$n1
    aa$n0 <- job$n0
    aa$trim <- cfg$ps_trim
    aa$auto_tune_seed <- cfg$auto_tune_seed + job$rep_id
    assign(".Random.seed", m$streams[[as.character(id)]], envir = .GlobalEnv)
    before <- study_rng()
    warnings <- character()
    fatal <- NA_character_
    start <- Sys.time()
    res <- tryCatch(withCallingHandlers(do.call(e$fit_one_rep, aa), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(ex) {
        fatal <<- conditionMessage(ex)
        NULL
    })
    if (is.null(res)) 
        res <- do.call(rbind, lapply(m$method_labels[[job$config_id]], function(label) e$empty_row(job$rep_id, job$n1, job$n0, label, NA_real_, paste("Fatal dataset error:", fatal))))
    stopifnot(nrow(res) == 8L, !anyDuplicated(res$Method), setequal(res$Method, m$method_labels[[job$config_id]]), all(is.na(res$theta_true)))
    list(job = job, stage = m$stage, source_md5 = m$source_md5, derivative = derivative, results = res, data_hash = e$.data_hash, PS_trace = e$.PS_trace, nuisance_trace = e$.nuisance_trace, calibration_trace = e$.calibration_trace, estimator_trace = e$.estimator_trace, math_trace = e$.math_trace, R_rng_before = before, R_rng_after = study_rng(), warnings = warnings, fatal_error = fatal, elapsed = as.numeric(difftime(Sys.time(), start, units = "secs")), finished = Sys.time())
}

study_validate <- function(z, id, m) {
    job <- m$grid[m$grid$id == id, , drop = FALSE]
    stopifnot(identical(z$stage, m$stage), identical(z$job, job), identical(z$source_md5, m$source_md5), identical(z$derivative, "analytic"), identical(z$R_rng_before, m$streams[[as.character(id)]]), nrow(z$results) == 8L, !anyDuplicated(z$results$Method), setequal(z$results$Method, m$method_labels[[job$config_id]]), all(is.na(z$results$theta_true)))
    invisible(TRUE)
}

