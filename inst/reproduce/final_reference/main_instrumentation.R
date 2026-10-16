main_hash <- function(x) digest::digest(x, algo = "sha256")

main_rng <- function() get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

main_install_math_logging <- function(e, capture_full = FALSE) {
    e$.math_trace <- list()
    e$.solve_trace <- list()
    e$.capture_math_full <- capture_full
    logged_solve <- function(A, b = NULL, ridge = 1e-08) {
        A <- as.matrix(A)
        path <- "direct"
        ans <- tryCatch({
            if (is.null(b)) 
                solve(A)
            else solve(A, b)
        }, error = function(ex) {
            path <<- "ridge"
            A2 <- A + diag(ridge, nrow(A))
            tryCatch({
                if (is.null(b)) 
                  solve(A2)
                else solve(A2, b)
            }, error = function(ex2) {
                path <<- "ginv"
                if (is.null(b)) 
                  MASS::ginv(A)
                else MASS::ginv(A) %*% b
            })
        })
        .solve_trace[[length(.solve_trace) + 1L]] <<- list(or_method = .current_or, path = path, ridge = ridge, dimension = dim(A), A = if (.capture_math_full) A else NULL, b = if (.capture_math_full) b else NULL, answer = if (.capture_math_full) ans else NULL)
        ans
    }
    environment(logged_solve) <- e
    e$solve_safe <- logged_solve
    wrapenv <- environment(e$fit_mec_cox)
    raw <- get("raw", envir = wrapenv, inherits = FALSE)
    bodylist <- as.list(body(raw))
    j <- which(vapply(bodylist, function(x) identical(x, quote(phi_theta <- IF[, 1])), logical(1)))
    stopifnot(length(j) == 1L)
    inserted <- quote(.math_trace[[length(.math_trace) + 1L]] <<- list(or_method = .current_or, a = D_theta_theta, b = D_theta_lambda, C = D_lambda_lambda, lambda = cal$lambda, positive_finite_weights = all(is.finite(w0_mec) & w0_mec > 0), exponential_weights_identical = identical(w0_mec, fam$ginv(as.numeric(gd0 + H0 %*% cal$lambda))), D_condition = tryCatch(suppressWarnings(kappa(D, exact = TRUE)), error = function(ex) NA_real_), solve_trace_index = length(.solve_trace), solve_path = .solve_trace[[length(.solve_trace)]]$path, 
        full = if (.capture_math_full) list(H0 = H0, eta_controls = eta_theta[idx0], W = W_mec, rho = rho, D = D, IF = IF, phi = phi_theta) else NULL))
    body(raw) <- as.call(append(bodylist, list(inserted), after = j))
    assign("raw", raw, envir = wrapenv)
    invisible(e)
}

main_make_engine <- function(figure, derivative = "analytic", capture_full = FALSE, extra_logging = TRUE) {
    ns <- new.env(parent = globalenv())
    sys.source(file.path(main_dir, "instrumentation_reference.R"), ns)
    ns$moderate_load_engine <- local({
        f <- figure
        k <- derivative
        function() main_load_engine(f, k)
    })
    e <- ns$moderate_make_engine(ns$moderate_candidates()$fixed100_k2, TRUE, FALSE, capture_full)
    e$.math_trace <- e$.solve_trace <- e$.point_trace <- e$.OR_tuning_trace <- list()
    if (!extra_logging) 
        return(e)
    main_install_math_logging(e, capture_full)
    if (figure == "figure5") 
        e$estimate_ps_ml <- local({
            raw <- e$estimate_ps_ml
            engine <- e
            function(...) {
                args <- list(...)
                ans <- do.call(raw, args)
                if (identical(args$ps_method, "glm") && !is.null(args$newdata)) {
                  p <- as.numeric(ans$ehat)
                  y <- as.integer(args$newdata$A)
                  engine$.PS_trace[[length(engine$.PS_trace) + 1L]] <- list(or_method = engine$.current_or, tune_seed = args$tune_seed, cache_hit = FALSE, n_train = nrow(args$data), n_pred = nrow(args$newdata), ids = args$newdata$id, y = y, p = p, tuned_params = ans$tuned_params, R_rng_after = main_rng())
                }
                ans
            }
        })
    e$tune_rsf_survival <- local({
        raw <- e$tune_rsf_survival
        engine <- e
        function(...) {
            args <- list(...)
            before <- main_rng()
            ans <- do.call(raw, args)
            engine$.OR_tuning_trace[[length(engine$.OR_tuning_trace) + 1L]] <- list(tune_seed = args$tune_seed, tune_n = args$tune_n, fast_rsf = args$fast_rsf, landmark_times = args$landmark_times, t0 = median(args$landmark_times), training_ids = args$train_data$id, xvars = args$xvars, selected_params = ans, R_rng_before = before, R_rng_after = main_rng())
            ans
        }
    })
    e$fit_mec_cox <- local({
        raw <- e$fit_mec_cox
        engine <- e
        full <- capture_full
        function(...) {
            args <- list(...)
            ans <- do.call(raw, args)
            point <- ans
            point$se <- NULL
            engine$.point_trace[[length(engine$.point_trace) + 1L]] <- list(or_method = args$or_method, point_hash = main_hash(point), W_hash = main_hash(ans$W), H_hash = main_hash(ans$H), theta = ans$theta, point = if (full) point else NULL)
            ans
        }
    })
    e
}

main_role <- function(method) ifelse(grepl("OR=cox", method, fixed = TRUE), "MEC_Cox", ifelse(grepl("OR=rsf", method, fixed = TRUE), "MEC_RSF", method))

