moderate_candidates <- function() list(reference = list(id = "reference", mode = "reference", k = 2, description = "Original25/50/100 tree tuning,100post/50burn; exact mean-first offset-only source"), fixed100_k2 = list(id = "fixed100_k2", mode = "fixed", ntree = 100L, k = 2, ndpost = 300L, nskip = 150L), fixed100_k1p5 = list(id = "fixed100_k1p5", mode = "fixed", ntree = 100L, k = 1.5, ndpost = 300L, nskip = 150L))

moderate_rng <- function() get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

moderate_hash <- function(x) digest::digest(x, algo = "sha256")

moderate_add_k <- function(f) {
    formals(f) <- c(formals(f), alist(bart_k = 2))
    walk <- function(x) {
        if (!is.call(x)) 
            return(x)
        if (identical(x[[1]], quote(dbarts::bart))) 
            x$k <- as.name("bart_k")
        for (j in seq_along(x)[-1L]) if (is.call(x[[j]])) 
            x[[j]] <- walk(x[[j]])
        x
    }
    body(f) <- walk(body(f))
    f
}

moderate_make_engine <- function(candidate, cache_enabled = TRUE, plain_reference = FALSE, capture_full = FALSE) {
    stopifnot(!plain_reference || candidate$id == "reference")
    e <- moderate_load_engine()
    e$.PS_cache <- new.env(parent = emptyenv())
    e$.PS_trace <- list()
    e$.nuisance_trace <- list()
    e$.calibration_trace <- list()
    e$.estimator_trace <- list()
    e$.current_or <- NA_character_
    e$.moderate_raw_se <- NA_real_
    e$.data_hash <- NULL
    generator <- e$generate_observed_data
    e$generate_observed_data <- local({
        raw <- generator
        engine <- e
        function(...) {
            ans <- raw(...)
            engine$.data_hash <- moderate_hash(ans)
            ans
        }
    })
    original_ps <- e$estimate_ps_ml
    fixed_ps <- moderate_add_k(original_ps)
    e$estimate_ps_ml <- local({
        raw <- original_ps
        fixed <- fixed_ps
        engine <- e
        cand <- candidate
        use_cache <- cache_enabled
        plain <- plain_reference
        function(...) {
            args <- list(...)
            outer <- identical(args$ps_method, "bart") && isTRUE(args$auto_tune) && !is.null(args$newdata)
            if (!outer) 
                return(do.call(raw, args))
            stopifnot(length(args$tune_seed) == 1L, is.finite(args$tune_seed))
            if (!plain) 
                set.seed(args$tune_seed)
            before <- moderate_rng()
            key_data <- list(candidate = cand, args = args, R_rng_before = before, dbarts = as.character(packageVersion("dbarts")))
            key_raw <- serialize(key_data, NULL, version = 3)
            key <- digest::digest(key_raw, serialize = FALSE, algo = "sha256")
            hit <- use_cache && !plain && exists(key, envir = engine$.PS_cache, inherits = FALSE)
            if (hit) {
                stored <- get(key, envir = engine$.PS_cache)
                stopifnot(identical(stored$key_raw, key_raw))
                ans <- stored$result
                assign(".Random.seed", stored$R_rng_after, envir = .GlobalEnv)
                for (w in stored$warnings) warning(w, call. = FALSE)
            }
            else {
                warnings <- character()
                ans <- withCallingHandlers({
                  if (cand$mode == "reference") 
                    do.call(raw, args)
                  else {
                    aa <- args
                    aa$auto_tune <- FALSE
                    aa$bart_ntree <- cand$ntree
                    aa$bart_ndpost <- cand$ndpost
                    aa$bart_nskip <- cand$nskip
                    aa$bart_k <- cand$k
                    fit <- do.call(fixed, aa)
                    fit$tuned_params <- list(bart_ntree = cand$ntree, bart_ndpost = cand$ndpost, bart_nskip = cand$nskip, bart_k = cand$k)
                    fit
                  }
                }, warning = function(w) warnings <<- c(warnings, conditionMessage(w)))
                stored <- list(key_raw = key_raw, result = ans, R_rng_after = moderate_rng(), warnings = warnings)
                if (use_cache && !plain) 
                  assign(key, stored, envir = engine$.PS_cache)
            }
            p <- as.numeric(ans$ehat)
            y <- as.integer(args$newdata$A)
            engine$.PS_trace[[length(engine$.PS_trace) + 1L]] <- list(or_method = engine$.current_or, tune_seed = args$tune_seed, cache_hit = hit, key = key, n_train = nrow(args$data), n_pred = nrow(args$newdata), ids = args$newdata$id, y = y, p = p, logloss = engine$binary_logloss(y, p), brier = mean((y - p)^2), predicted_prevalence = mean(p), observed_prevalence = mean(y), clipping_rate = mean(p <= args$trim[1L] | p >= args$trim[2L]), tuned_params = ans$tuned_params, R_rng_after = stored$R_rng_after)
            ans
        }
    })
    nuisance <- e$estimate_mec_nuisance
    e$estimate_mec_nuisance <- local({
        raw <- nuisance
        engine <- e
        full <- capture_full
        function(...) {
            args <- list(...)
            ans <- do.call(raw, args)
            engine$.nuisance_trace[[length(engine$.nuisance_trace) + 1L]] <- list(or_method = args$or_method, ehat = ans$ehat, H_hash = moderate_hash(ans$H), H = if (full) ans$H else NULL, fold_id = ans$fold_id, landmark_times = ans$landmark_times, actual_landmarks = ncol(ans$H) - 1L)
            ans
        }
    })
    calibration <- e$calibrate_weights_bregman_from_design_return
    e$calibrate_weights_bregman_from_design_return <- local({
        raw <- calibration
        engine <- e
        function(...) {
            args <- list(...)
            Xs <- args$Xs
            Xp <- args$Xp
            condition <- function(x) tryCatch(suppressWarnings(kappa(x, exact = TRUE)), error = function(e) NA_real_)
            row <- list(or_method = engine$.current_or, constraints = ncol(Xs), design_rank = tryCatch(qr(Xs)$rank, error = function(e) NA_integer_), design_condition = condition(Xs), D_lambda_condition = NA_real_, gradient_norm = NA_real_, converged = NA, finite_weights = NA, weight_min = NA_real_, weight_max = NA_real_, weight_sum = NA_real_, max_mean_balance_residual = NA_real_, error = NA_character_)
            ans <- tryCatch(do.call(raw, args), error = function(ex) ex)
            if (inherits(ans, "error")) {
                row$error <- conditionMessage(ans)
                engine$.calibration_trace[[length(engine$.calibration_trace) + 1L]] <- row
                stop(ans)
            }
            row$D_lambda_condition <- condition(ans$D_lambda_lambda)
            row$gradient_norm <- ans$grad_norm
            row$converged <- ans$converged
            row$finite_weights <- all(is.finite(ans$w))
            if (any(is.finite(ans$w))) {
                row$weight_min <- min(ans$w[is.finite(ans$w)])
                row$weight_max <- max(ans$w[is.finite(ans$w)])
            }
            row$weight_sum <- sum(ans$w)
            row$max_mean_balance_residual <- max(abs(colSums(Xs * ans$w) - colSums(Xp)))/nrow(Xp)
            engine$.calibration_trace[[length(engine$.calibration_trace) + 1L]] <- row
            ans
        }
    })
    fit <- e$fit_mec_cox
    b <- as.list(body(fit))
    j <- which(vapply(b, function(x) identical(x, quote(se_mec <- sqrt(sum(phi_theta^2)))), logical(1)))
    stopifnot(length(j) == 1L)
    body(fit) <- as.call(append(b, list(quote(.moderate_raw_se <<- se_mec)), after = j))
    environment(fit) <- e
    e$fit_mec_cox <- local({
        raw <- fit
        engine <- e
        function(...) {
            args <- list(...)
            engine$.current_or <- args$or_method
            engine$.moderate_raw_se <- NA_real_
            ans <- tryCatch(do.call(raw, args), error = function(ex) ex)
            row <- list(or_method = args$or_method, raw_sandwich_SE = engine$.moderate_raw_se, reported_SE = NA_real_, fallback_used = NA, error = NA_character_)
            if (inherits(ans, "error")) {
                row$error <- conditionMessage(ans)
                engine$.estimator_trace[[length(engine$.estimator_trace) + 1L]] <- row
                stop(ans)
            }
            row$reported_SE <- ans$se
            row$fallback_used <- !is.finite(engine$.moderate_raw_se)
            engine$.estimator_trace[[length(engine$.estimator_trace) + 1L]] <- row
            ans
        }
    })
    e
}

