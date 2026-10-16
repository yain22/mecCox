expit <- function(x) {
    1/(1 + exp(-x))
}

clip01 <- function(x, lower = 0.01, upper = 0.99) {
    pmin(pmax(as.numeric(x), lower), upper)
}

round_numeric_df <- function(x, digits = 4) {
    x <- as.data.frame(x)
    num_cols <- vapply(x, is.numeric, logical(1))
    x[num_cols] <- lapply(x[num_cols], round, digits = digits)
    x
}

solve_safe <- function(A, b = NULL, ridge = 1e-08) {
    A <- as.matrix(A)
    ans <- tryCatch({
        if (is.null(b)) {
            solve(A)
        }
        else {
            solve(A, b)
        }
    }, error = function(e) {
        A2 <- A + diag(ridge, nrow(A))
        tryCatch({
            if (is.null(b)) {
                solve(A2)
            }
            else {
                solve(A2, b)
            }
        }, error = function(e2) {
            if (is.null(b)) {
                MASS::ginv(A)
            }
            else {
                MASS::ginv(A) %*% b
            }
        })
    })
    ans
}

jacobian_fd <- function(f, x, eps = 1e-05) {
    p <- length(x)
    f0 <- f(x)
    J <- matrix(NA_real_, nrow = length(f0), ncol = p)
    for (j in seq_len(p)) {
        step <- rep(0, p)
        step[j] <- eps * max(1, abs(x[j]))
        J[, j] <- (f(x + step) - f(x - step))/(2 * step[j])
    }
    J
}

make_stratified_folds <- function(A, K = 5, seed = 20260427) {
    set.seed(seed)
    A <- as.integer(A)
    fold_id <- integer(length(A))
    for (a in sort(unique(A))) {
        idx <- which(A == a)
        fold_id[idx] <- sample(rep(seq_len(K), length.out = length(idx)))
    }
    fold_id
}

source_score <- function(X, alpha0 = -0.2, ps_truth = c("linear", "nonlinear"), ps_nonlinearity = 1) {
    ps_truth <- match.arg(ps_truth)
    X <- as.matrix(X)
    x1 <- X[, 1]
    x2 <- X[, 2]
    x3 <- X[, 3]
    x4 <- X[, 4]
    x5 <- X[, 5]
    lin_eta <- 0.75 * x1 + 0.75 * x2 + 0.65 * x3 + 0.65 * x4 + 0.55 * x5
    nonlin_eta <- 0.7 * sin(1.25 * x1) + 0.45 * (x2^2 - 1) - 0.55 * (as.numeric(x3 > 0) - 0.5) + 0.35 * x4 * x5 + 0.25 * (cos(x1 + x2) - exp(-1))
    if (ps_truth == "linear") {
        eta <- alpha0 + lin_eta
    }
    else {
        eta <- alpha0 + lin_eta + ps_nonlinearity * nonlin_eta
    }
    clip01(expit(eta), 0.02, 0.98)
}

sample_source_covariates <- function(n1, n0, p = 10, alpha0 = -0.2, ps_truth = "linear", ps_nonlinearity = 1) {
    X1 <- matrix(NA_real_, nrow = 0, ncol = p)
    X0 <- matrix(NA_real_, nrow = 0, ncol = p)
    while (nrow(X1) < n1 || nrow(X0) < n0) {
        m <- max(4000, 4 * (n1 + n0))
        Xcand <- matrix(rnorm(m * p), nrow = m, ncol = p)
        pi <- source_score(Xcand, alpha0 = alpha0, ps_truth = ps_truth, ps_nonlinearity = ps_nonlinearity)
        A <- rbinom(m, size = 1, prob = pi)
        add1 <- Xcand[A == 1, , drop = FALSE]
        add0 <- Xcand[A == 0, , drop = FALSE]
        if (nrow(X1) < n1) 
            X1 <- rbind(X1, add1)
        if (nrow(X0) < n0) 
            X0 <- rbind(X0, add0)
    }
    X1 <- X1[seq_len(n1), , drop = FALSE]
    X0 <- X0[seq_len(n0), , drop = FALSE]
    colnames(X1) <- paste0("X", seq_len(p))
    colnames(X0) <- paste0("X", seq_len(p))
    list(X1 = X1, X0 = X0)
}

domain_shift <- function(X, type = c("constant", "nonlinear")) {
    type <- match.arg(type)
    X <- as.matrix(X)
    if (type == "constant") {
        return(rep(1, nrow(X)))
    }
    1 + rowSums(sin(2 * X[, 1:5, drop = FALSE]))/sqrt(5)
}

prognostic_lp <- function(X, outcome_truth = c("linear_ph", "nonlinear_ph"), or_nonlinearity = 1, p_signal_outcome = 10) {
    outcome_truth <- match.arg(outcome_truth)
    X <- as.matrix(X)
    p <- ncol(X)
    if (p < 5) {
        stop("p must be at least 5 because the nonlinear prognostic component uses X1,...,X5.")
    }
    b_base <- c(log(1.75), log(1.75), log(1.6), log(1.6), log(1.5), log(1.25), log(1.25), log(1.25), log(1.25), log(1.25))
    b <- rep(0, p)
    J <- min(length(b_base), p_signal_outcome, p)
    b[seq_len(J)] <- b_base[seq_len(J)]
    lin_lp <- drop(X %*% b)
    x1 <- X[, 1]
    x2 <- X[, 2]
    x3 <- X[, 3]
    x4 <- X[, 4]
    x5 <- X[, 5]
    nonlin_lp <- 0.45 * sin(x2) + 0.35 * (x3^2 - 1) + 0.3 * (as.numeric(x4 > 0) - 0.5) + 0.25 * x1 * x5 + 0.2 * (cos(x2 + x5) - exp(-1))
    if (outcome_truth == "linear_ph") {
        lp <- lin_lp
    }
    else {
        lp <- lin_lp + or_nonlinearity * nonlin_lp
    }
    lp
}

weibull_time <- function(lp, lambda0 = 8e-05, eta_shape = 2) {
    U <- runif(length(lp))
    (-log(U)/(lambda0 * exp(lp)))^(1/eta_shape)
}

generate_observed_data <- function(n1 = 250, n0 = 500, p = 10, alpha0 = -0.2, ps_truth = "linear", outcome_truth = "linear_ph", ps_nonlinearity = 1, or_nonlinearity = 1, beta_cond = log(0.7), xi = 0, h_type = "constant", lambda0 = 8e-05, eta_shape = 2, censor_rate = 8e-04) {
    SX <- sample_source_covariates(n1 = n1, n0 = n0, p = p, alpha0 = alpha0, ps_truth = ps_truth, ps_nonlinearity = ps_nonlinearity)
    X1 <- SX$X1
    X0 <- SX$X0
    lp0_trial_ext <- prognostic_lp(X1, outcome_truth = outcome_truth, or_nonlinearity = or_nonlinearity)
    lp0_trial_true <- lp0_trial_ext + xi * domain_shift(X1, h_type)
    lp1_trial_true <- lp0_trial_true + beta_cond
    lp0_external <- prognostic_lp(X0, outcome_truth = outcome_truth, or_nonlinearity = or_nonlinearity)
    T1_trial <- weibull_time(lp1_trial_true, lambda0, eta_shape)
    T0_external <- weibull_time(lp0_external, lambda0, eta_shape)
    C1 <- rexp(n1, rate = censor_rate)
    C0 <- rexp(n0, rate = censor_rate)
    dat1 <- data.frame(A = 1L, time = pmin(T1_trial, C1), delta = as.integer(T1_trial <= C1), X1)
    dat0 <- data.frame(A = 0L, time = pmin(T0_external, C0), delta = as.integer(T0_external <= C0), X0)
    dat <- rbind(dat1, dat0)
    dat$id <- seq_len(nrow(dat))
    dat
}

make_att_ipw_data_glm <- function(data, xvars, trim = c(0.01, 0.99), normalize_controls = TRUE) {
    d <- as.data.frame(data)
    ps_fit <- glm(reformulate(xvars, response = "A"), data = d, family = binomial())
    ehat <- clip01(predict(ps_fit, type = "response"), trim[1], trim[2])
    qhat <- ehat/(1 - ehat)
    if (normalize_controls) {
        n1 <- sum(d$A == 1)
        qhat[d$A == 0] <- qhat[d$A == 0] * n1/sum(qhat[d$A == 0])
    }
    W <- with(d, A + (1 - A) * qhat)
    d$ehat_glm <- ehat
    d$qhat_glm <- qhat
    d$W_att <- W
    list(data = d, ps_fit = ps_fit, ehat = ehat, qhat = qhat, W = W)
}

att_weights_from_gamma <- function(gamma, data, xvars, trim = c(0.01, 0.99), normalize_controls = TRUE) {
    d <- as.data.frame(data)
    Z <- model.matrix(reformulate(xvars, response = "A"), data = d)
    ehat <- clip01(expit(as.vector(Z %*% gamma)), trim[1], trim[2])
    qhat <- ehat/(1 - ehat)
    if (normalize_controls) {
        n1 <- sum(d$A == 1)
        qhat[d$A == 0] <- qhat[d$A == 0] * n1/sum(qhat[d$A == 0])
    }
    with(d, A + (1 - A) * qhat)
}

cox_score_eta_fixed_weight <- function(theta, data, W) {
    d <- as.data.frame(data)
    A <- d$A
    time <- d$time
    delta <- d$delta
    n <- nrow(d)
    r <- exp(theta * A)
    ord_desc <- order(time, decreasing = TRUE)
    time_desc <- time[ord_desc]
    wr_desc <- W[ord_desc] * r[ord_desc]
    wrA_desc <- W[ord_desc] * r[ord_desc] * A[ord_desc]
    cum_S0 <- cumsum(wr_desc)
    cum_S1 <- cumsum(wrA_desc)
    rle_time <- rle(time_desc)
    times_desc <- rle_time$values
    group_end <- cumsum(rle_time$lengths)
    S0_desc <- cum_S0[group_end]
    S1_desc <- cum_S1[group_end]
    time_index <- match(time, times_desc)
    S0_at <- S0_desc[time_index]
    S1_at <- S1_desc[time_index]
    Abar_at <- S1_at/S0_at
    direct <- W * delta * (A - Abar_at)
    Ucox <- sum(direct)
    event_times_asc <- sort(unique(time[delta == 1]))
    D0_event <- as.numeric(tapply(W * delta, factor(time, levels = event_times_asc), sum))
    event_index_desc <- match(event_times_asc, times_desc)
    S0_event <- S0_desc[event_index_desc]
    S1_event <- S1_desc[event_index_desc]
    Abar_event <- S1_event/S0_event
    q1_event <- D0_event/S0_event
    qbar_event <- D0_event * Abar_event/S0_event
    cum_q1 <- cumsum(q1_event)
    cum_qbar <- cumsum(qbar_event)
    pos_event <- findInterval(time, event_times_asc)
    Q1_i <- numeric(n)
    Qbar_i <- numeric(n)
    has_event_before <- pos_event > 0
    Q1_i[has_event_before] <- cum_q1[pos_event[has_event_before]]
    Qbar_i[has_event_before] <- cum_qbar[pos_event[has_event_before]]
    risk_correction <- W * r * (A * Q1_i - Qbar_i)
    eta <- direct - risk_correction
    list(U = Ucox, eta = eta)
}

manual_linwei_se <- function(theta_hat, data, W) {
    ce <- cox_score_eta_fixed_weight(theta_hat, data, W)
    D_theta <- as.numeric(jacobian_fd(function(th) {
        c(cox_score_eta_fixed_weight(th[1], data, W)$U)
    }, x = c(theta_hat)))
    phi <- -as.numeric(ce$eta/D_theta)
    sqrt(sum(phi^2))
}

fit_unweighted_cox <- function(data) {
    d <- as.data.frame(data)
    fit <- coxph(Surv(time, delta) ~ A, data = d, robust = FALSE, ties = "breslow")
    theta <- unname(coef(fit)["A"])
    se <- sqrt(as.numeric(vcov(fit)["A", "A"]))
    list(theta = theta, se = se)
}

fit_att_ipw_naive_linwei_shu <- function(data, xvars, trim = c(0.01, 0.99), normalize_controls = TRUE) {
    obj <- make_att_ipw_data_glm(data = data, xvars = xvars, trim = trim, normalize_controls = normalize_controls)
    d <- obj$data
    W <- obj$W
    fit_naive <- coxph(Surv(time, delta) ~ A, data = d, weights = W, robust = FALSE, ties = "breslow")
    theta_hat <- unname(coef(fit_naive)["A"])
    se_naive <- sqrt(as.numeric(vcov(fit_naive)["A", "A"]))
    se_lw <- manual_linwei_se(theta_hat, d, W)
    ps_fit <- obj$ps_fit
    gamma_hat <- coef(ps_fit)
    Z <- model.matrix(reformulate(xvars, response = "A"), data = d)
    A <- d$A
    ehat <- obj$ehat
    ce <- cox_score_eta_fixed_weight(theta_hat, d, W)
    eta_theta <- ce$eta
    score_gamma <- sweep(Z, 1, A - ehat, "*")
    D_theta_theta <- as.numeric(jacobian_fd(function(th) {
        c(cox_score_eta_fixed_weight(th[1], d, W)$U)
    }, x = c(theta_hat)))
    D_theta_gamma <- jacobian_fd(function(gam) {
        W_gam <- att_weights_from_gamma(gamma = gam, data = d, xvars = xvars, trim = trim, normalize_controls = normalize_controls)
        c(cox_score_eta_fixed_weight(theta_hat, d, W_gam)$U)
    }, x = gamma_hat)
    D_gamma_gamma <- -crossprod(Z * (ehat * (1 - ehat)), Z)
    K <- length(gamma_hat)
    D <- rbind(c(D_theta_theta, as.numeric(D_theta_gamma)), cbind(rep(0, K), D_gamma_gamma))
    U_i <- cbind(eta_theta, score_gamma)
    IF <- -U_i %*% t(solve_safe(D))
    phi_theta <- IF[, 1]
    se_shu <- sqrt(sum(phi_theta^2))
    list(theta = theta_hat, se_naive = se_naive, se_lw = se_lw, se_shu = se_shu, W = W, ehat = ehat)
}

.to_matrix <- function(X) {
    if (is.null(dim(X))) 
        matrix(as.numeric(X), ncol = 1)
    else as.matrix(X)
}

.bregman_family <- function(name, alpha = 1/2) {
    name <- tolower(name)
    switch(name, quadratic = list(g = function(w) w, ginv = function(nu) nu, gprime = function(w) rep(1, length(w)), ok = function(nu) all(is.finite(nu))), kl = list(g = function(w) log(pmax(w, 1e-16)) + 1, ginv = function(nu) exp(nu - 1), gprime = function(w) 1/pmax(w, 1e-16), ok = function(nu) all(is.finite(nu))), el = list(g = function(w) -1/pmax(w, 1e-16), ginv = function(nu) -1/nu, gprime = function(w) 1/pmax(w, 1e-16)^2, ok = function(nu) all(is.finite(nu)) && all(nu < 0)), hellinger = list(g = function(w) 1 - 
        1/sqrt(pmax(w, 1e-16)), ginv = function(nu) 1/(1 - nu)^2, gprime = function(w) 0.5/pmax(w, 1e-16)^(3/2), ok = function(nu) all(is.finite(nu)) && all(nu < 1)), renyi = list(g = function(w) pmax(w, 1e-16)^alpha, ginv = function(nu) pmax(nu, 1e-16)^(1/alpha), gprime = function(w) alpha * pmax(w, 1e-16)^(alpha - 1), ok = function(nu) all(is.finite(nu)) && all(nu > 0)), stop("Unknown divergence."))
}

.weight_diag <- function(w) {
    n <- length(w)
    s1 <- sum(w)
    s2 <- sum(w^2)
    ess <- (s1^2)/s2
    rel_ess <- ess/n
    cvw <- sd(w)/mean(w)
    list(ESS = ess, Rel_ESS = rel_ess, W_CV = cvw, W_Min = min(w), W_Max = max(w), SumW = s1)
}

calibrate_weights_bregman_from_design_return <- function(Xs, Xp, d, divergence = c("kl", "quadratic", "el", "hellinger", "renyi"), alpha = 1/2, maxit = 100, tol = 1e-10, ridge = 1e-08) {
    divergence <- match.arg(divergence)
    fam <- .bregman_family(divergence, alpha)
    Xs <- .to_matrix(Xs)
    Xp <- .to_matrix(Xp)
    Tx <- colSums(Xp)
    if (length(d) == 1) 
        d <- rep(d, nrow(Xs))
    gd <- fam$g(d)
    lam <- rep(0, ncol(Xs))
    step_shrink <- 0.5
    converged <- FALSE
    for (it in seq_len(maxit)) {
        nu <- as.numeric(gd + Xs %*% lam)
        bt <- 0
        while (!fam$ok(nu) && bt < 30) {
            lam <- lam * step_shrink
            nu <- as.numeric(gd + Xs %*% lam)
            bt <- bt + 1
        }
        w <- fam$ginv(nu)
        grad <- as.numeric(crossprod(Xs, w) - Tx)
        grad_norm <- sqrt(sum(grad^2))
        if (grad_norm < tol) {
            converged <- TRUE
            break
        }
        dwdnu <- 1/fam$gprime(w)
        Hmat <- crossprod(Xs * dwdnu, Xs) + diag(ridge, ncol(Xs))
        step <- tryCatch(solve(Hmat, grad), error = function(e) MASS::ginv(Hmat) %*% grad)
        t <- 1
        repeat {
            lam_new <- lam - t * as.numeric(step)
            nu_new <- as.numeric(gd + Xs %*% lam_new)
            if (fam$ok(nu_new)) {
                lam <- lam_new
                break
            }
            t <- t * step_shrink
            if (t < 1e-08) {
                lam <- lam_new
                break
            }
        }
    }
    nu <- as.numeric(gd + Xs %*% lam)
    w <- fam$ginv(nu)
    grad <- as.numeric(crossprod(Xs, w) - Tx)
    grad_norm <- sqrt(sum(grad^2))
    dwdnu <- 1/fam$gprime(w)
    D_lambda_lambda <- crossprod(Xs * dwdnu, Xs)
    list(w = as.numeric(w), lambda = as.numeric(lam), grad = grad, grad_norm = grad_norm, D_lambda_lambda = D_lambda_lambda, converged = converged, iterations = it, diagnostics = .weight_diag(w))
}

make_balanced_tuning_indices <- function(A, tune_n = 100) {
    A <- as.integer(A)
    n <- length(A)
    if (n <= tune_n) {
        return(seq_len(n))
    }
    idx1 <- which(A == 1)
    idx0 <- which(A == 0)
    if (length(idx1) == 0 || length(idx0) == 0) {
        return(seq_len(min(tune_n, n)))
    }
    n1_target <- min(length(idx1), floor(tune_n/2))
    n0_target <- min(length(idx0), tune_n - n1_target)
    left <- tune_n - n1_target - n0_target
    if (left > 0 && length(idx1) > n1_target) {
        add <- min(left, length(idx1) - n1_target)
        n1_target <- n1_target + add
        left <- left - add
    }
    if (left > 0 && length(idx0) > n0_target) {
        add <- min(left, length(idx0) - n0_target)
        n0_target <- n0_target + add
    }
    sort(c(head(idx1, n1_target), head(idx0, n0_target)))
}

make_stratified_train_valid_split <- function(A, train_frac = 0.7, seed = 20260427) {
    set.seed(seed)
    A <- as.integer(A)
    n <- length(A)
    idx_train <- integer(0)
    idx_valid <- integer(0)
    for (a in sort(unique(A))) {
        idx_a <- which(A == a)
        if (length(idx_a) <= 2) {
            idx_train <- c(idx_train, idx_a)
        }
        else {
            n_train_a <- max(1, floor(train_frac * length(idx_a)))
            train_a <- sample(idx_a, size = n_train_a, replace = FALSE)
            valid_a <- setdiff(idx_a, train_a)
            idx_train <- c(idx_train, train_a)
            idx_valid <- c(idx_valid, valid_a)
        }
    }
    if (length(idx_valid) == 0) {
        idx_all <- seq_len(n)
        n_train <- max(1, floor(train_frac * n))
        idx_train <- sample(idx_all, size = n_train, replace = FALSE)
        idx_valid <- setdiff(idx_all, idx_train)
    }
    list(train = sort(idx_train), valid = sort(idx_valid))
}

binary_logloss <- function(y, p, eps = 1e-06) {
    p <- pmin(pmax(as.numeric(p), eps), 1 - eps)
    y <- as.integer(y)
    -mean(y * log(p) + (1 - y) * log(1 - p))
}

safe_metric_min <- function(x) {
    x <- as.numeric(x)
    if (all(!is.finite(x))) {
        return(NA_integer_)
    }
    which.min(ifelse(is.finite(x), x, Inf))
}

safe_metric_max <- function(x) {
    x <- as.numeric(x)
    if (all(!is.finite(x))) {
        return(NA_integer_)
    }
    which.max(ifelse(is.finite(x), x, -Inf))
}

tune_ps_ml <- function(data, xvars, ps_method = c("glm", "glm_lasso", "nnet", "rf", "bart"), trim = c(0.01, 0.99), tune_n = 100, tune_seed = 20260427, tune_verbose = FALSE) {
    ps_method <- match.arg(ps_method)
    if (ps_method %in% c("glm", "glm_lasso")) {
        return(list())
    }
    d <- as.data.frame(data)
    idx_tune <- make_balanced_tuning_indices(d$A, tune_n = tune_n)
    d_tune <- d[idx_tune, , drop = FALSE]
    if (length(unique(d_tune$A)) < 2 || nrow(d_tune) < 20) {
        if (isTRUE(tune_verbose)) {
            message("PS auto-tuning skipped: too few observations or only one class.")
        }
        return(list())
    }
    sp <- make_stratified_train_valid_split(A = d_tune$A, train_frac = 0.7, seed = tune_seed)
    d_train <- d_tune[sp$train, , drop = FALSE]
    d_valid <- d_tune[sp$valid, , drop = FALSE]
    if (length(unique(d_train$A)) < 2 || length(unique(d_valid$A)) < 2) {
        if (isTRUE(tune_verbose)) {
            message("PS auto-tuning skipped: train/validation split has one class.")
        }
        return(list())
    }
    p_dim <- length(xvars)
    if (ps_method == "nnet") {
        grid <- expand.grid(nnet_size = c(2, 4, 6), nnet_decay = c(1e-04, 0.001, 0.01), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    }
    if (ps_method == "rf") {
        grid <- expand.grid(rf_mtry = sort(unique(pmax(1, c(floor(sqrt(p_dim)), ceiling(p_dim/2), p_dim)))), rf_min_node_size = c(5, 10, 20), rf_num_trees = c(200), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    }
    if (ps_method == "bart") {
        grid <- expand.grid(bart_ntree = c(25, 50, 100), bart_ndpost = c(100), bart_nskip = c(50), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    }
    loss <- rep(Inf, nrow(grid))
    for (gg in seq_len(nrow(grid))) {
        pars <- as.list(grid[gg, , drop = FALSE])
        ans <- tryCatch({
            pred <- estimate_ps_ml(data = d_train, xvars = xvars, ps_method = ps_method, newdata = d_valid, trim = trim, nnet_size = if (!is.null(pars$nnet_size)) 
                pars$nnet_size
            else 4, nnet_decay = if (!is.null(pars$nnet_decay)) 
                pars$nnet_decay
            else 0.001, rf_num_trees = if (!is.null(pars$rf_num_trees)) 
                pars$rf_num_trees
            else 300, rf_mtry = if (!is.null(pars$rf_mtry)) 
                pars$rf_mtry
            else NULL, rf_min_node_size = if (!is.null(pars$rf_min_node_size)) 
                pars$rf_min_node_size
            else 5, bart_ntree = if (!is.null(pars$bart_ntree)) 
                pars$bart_ntree
            else 50, bart_ndpost = if (!is.null(pars$bart_ndpost)) 
                pars$bart_ndpost
            else 200, bart_nskip = if (!is.null(pars$bart_nskip)) 
                pars$bart_nskip
            else 100, auto_tune = FALSE)
            binary_logloss(d_valid$A, pred$ehat)
        }, error = function(e) Inf)
        loss[gg] <- ans
    }
    best_id <- safe_metric_min(loss)
    if (is.na(best_id)) {
        if (isTRUE(tune_verbose)) {
            message("PS auto-tuning failed for all candidates; using defaults.")
        }
        return(list())
    }
    best <- as.list(grid[best_id, , drop = FALSE])
    if (isTRUE(tune_verbose)) {
        message("PS auto-tuning selected for ", ps_method, ": ", paste(names(best), unlist(best), sep = "=", collapse = ", "), "; validation log-loss = ", signif(loss[best_id], 4))
    }
    best
}

estimate_ps_ml <- function(data, xvars, ps_method = c("glm", "glm_lasso", "nnet", "rf", "bart"), newdata = NULL, trim = c(0.01, 0.99), nnet_size = 4, nnet_decay = 0.001, nnet_maxit = 200, rf_num_trees = 300, rf_mtry = NULL, rf_min_node_size = 5, bart_ntree = 50, bart_ndpost = 200, bart_nskip = 100, ps_lasso_nfolds = 5, ps_lasso_lambda_rule = c("lambda.1se", "lambda.min"), ps_lasso_standardize = TRUE, auto_tune = FALSE, tune_n = 100, tune_seed = 20260427, tune_verbose = FALSE) {
    ps_method <- match.arg(ps_method)
    ps_lasso_lambda_rule <- match.arg(ps_lasso_lambda_rule)
    d_train <- as.data.frame(data)
    d_pred <- if (is.null(newdata)) 
        d_train
    else as.data.frame(newdata)
    tuned_params <- list()
    if (isTRUE(auto_tune) && ps_method %in% c("nnet", "rf", "bart")) {
        tuned_params <- tune_ps_ml(data = d_train, xvars = xvars, ps_method = ps_method, trim = trim, tune_n = tune_n, tune_seed = tune_seed, tune_verbose = tune_verbose)
        if (!is.null(tuned_params$nnet_size)) 
            nnet_size <- tuned_params$nnet_size
        if (!is.null(tuned_params$nnet_decay)) 
            nnet_decay <- tuned_params$nnet_decay
        if (!is.null(tuned_params$rf_num_trees)) 
            rf_num_trees <- tuned_params$rf_num_trees
        if (!is.null(tuned_params$rf_mtry)) 
            rf_mtry <- tuned_params$rf_mtry
        if (!is.null(tuned_params$rf_min_node_size)) 
            rf_min_node_size <- tuned_params$rf_min_node_size
        if (!is.null(tuned_params$bart_ntree)) 
            bart_ntree <- tuned_params$bart_ntree
        if (!is.null(tuned_params$bart_ndpost)) 
            bart_ndpost <- tuned_params$bart_ndpost
        if (!is.null(tuned_params$bart_nskip)) 
            bart_nskip <- tuned_params$bart_nskip
    }
    X_train <- as.matrix(d_train[, xvars, drop = FALSE])
    X_pred <- as.matrix(d_pred[, xvars, drop = FALSE])
    A_train <- as.integer(d_train$A)
    if (ps_method == "glm") {
        fit <- glm(reformulate(xvars, response = "A"), data = d_train, family = binomial())
        ehat <- predict(fit, newdata = d_pred, type = "response")
    }
    if (ps_method == "glm_lasso") {
        if (!requireNamespace("glmnet", quietly = TRUE)) {
            stop("Package {glmnet} is required for ps_method='glm_lasso'.")
        }
        if (length(unique(A_train)) < 2) {
            stop("The lasso-logistic PS fit needs both A=0 and A=1 in the training data.")
        }
        nfolds_use <- min(as.integer(ps_lasso_nfolds), length(A_train), sum(A_train == 1), sum(A_train == 0))
        if (!is.finite(nfolds_use) || nfolds_use < 3) {
            stop("Not enough treated/control observations for cv.glmnet lasso-logistic PS tuning.")
        }
        foldid <- make_stratified_folds(A_train, K = nfolds_use, seed = tune_seed)
        cvfit <- glmnet::cv.glmnet(x = X_train, y = A_train, family = "binomial", alpha = 1, foldid = foldid, type.measure = "deviance", standardize = ps_lasso_standardize, maxit = 1e+05, parallel = FALSE)
        lambda_selected <- if (ps_lasso_lambda_rule == "lambda.min") {
            cvfit$lambda.min
        }
        else {
            cvfit$lambda.1se
        }
        eta_hat <- as.numeric(stats::predict(cvfit, newx = X_pred, s = lambda_selected, type = "link"))
        ehat <- stats::plogis(eta_hat)
        beta_hat <- as.matrix(stats::coef(cvfit, s = lambda_selected))
        tuned_params <- c(tuned_params, list(ps_lasso_lambda_selected = lambda_selected, ps_lasso_nfolds = nfolds_use, ps_lasso_n_selected = sum(abs(beta_hat[, 1]) > 1e-12)))
    }
    if (ps_method == "nnet") {
        if (!requireNamespace("nnet", quietly = TRUE)) {
            stop("Package {nnet} is required for ps_method='nnet'.")
        }
        mu <- colMeans(X_train)
        sdv <- apply(X_train, 2, sd)
        sdv[sdv == 0] <- 1
        Xs_train <- scale(X_train, center = mu, scale = sdv)
        Xs_pred <- scale(X_pred, center = mu, scale = sdv)
        fit <- nnet::nnet(x = Xs_train, y = A_train, size = nnet_size, decay = nnet_decay, maxit = nnet_maxit, entropy = TRUE, trace = FALSE, MaxNWts = 10000)
        ehat <- as.numeric(predict(fit, Xs_pred, type = "raw"))
    }
    if (ps_method == "rf") {
        if (!requireNamespace("ranger", quietly = TRUE)) {
            stop("Package {ranger} is required for ps_method='rf'.")
        }
        df_rf <- d_train[, c("A", xvars), drop = FALSE]
        df_rf$A <- factor(df_rf$A, levels = c(0, 1))
        rf_mtry_use <- if (is.null(rf_mtry)) {
            max(1, floor(sqrt(length(xvars))))
        }
        else {
            max(1, min(length(xvars), as.integer(rf_mtry)))
        }
        fit <- ranger::ranger(A ~ ., data = df_rf, probability = TRUE, num.trees = rf_num_trees, mtry = rf_mtry_use, min.node.size = rf_min_node_size, num.threads = 1, seed = 20260427)
        pred <- predict(fit, data = d_pred[, xvars, drop = FALSE])$predictions
        ehat <- pred[, "1"]
    }
    if (ps_method == "bart") {
        if (!requireNamespace("dbarts", quietly = TRUE)) {
            stop("Package {dbarts} is required for ps_method='bart'.")
        }
        pA <- clip01(mean(A_train), 0.01, 0.99)
        binary_offset <- qnorm(pA)
        fit <- tryCatch(dbarts::bart(x.train = X_train, y.train = A_train, x.test = X_pred, ntree = bart_ntree, ndpost = bart_ndpost, nskip = bart_nskip, binaryOffset = binary_offset, verbose = FALSE, keeptrees = FALSE, nthread = 1), error = function(e) {
            dbarts::bart(x.train = X_train, y.train = A_train, x.test = X_pred, ntree = bart_ntree, ndpost = bart_ndpost, nskip = bart_nskip, verbose = FALSE, keeptrees = FALSE)
        })
        if (!is.null(fit$prob.test.mean)) {
            ehat <- as.numeric(fit$prob.test.mean)
        }
        else if (!is.null(fit$prob.test)) {
            ehat <- as.numeric(colMeans(fit$prob.test))
        }
        else if (!is.null(fit[["yhat.test"]])) {
            ehat <- as.numeric(colMeans(pnorm(fit[["yhat.test"]])))
        }
        else if (!is.null(fit$yhat.test.mean)) {
            ehat <- pnorm(as.numeric(fit$yhat.test.mean))
        }
        else {
            stop("Could not extract BART predicted probabilities.")
        }
    }
    ehat <- clip01(ehat, trim[1], trim[2])
    list(ehat = ehat, ps_method = ps_method, tuned_params = tuned_params)
}

fit_att_ipw_linwei_ml <- function(data, xvars, ps_method = c("glm", "glm_lasso", "nnet", "rf", "bart"), trim = c(0.01, 0.99), normalize_controls = TRUE, auto_tune = FALSE, tune_n = 100, tune_seed = 20260427, tune_verbose = FALSE, ps_lasso_nfolds = 5, ps_lasso_lambda_rule = c("lambda.1se", "lambda.min"), ps_lasso_standardize = TRUE) {
    ps_method <- match.arg(ps_method)
    ps_lasso_lambda_rule <- match.arg(ps_lasso_lambda_rule)
    d <- as.data.frame(data)
    ps_obj <- estimate_ps_ml(data = d, xvars = xvars, ps_method = ps_method, newdata = d, trim = trim, auto_tune = auto_tune, tune_n = tune_n, tune_seed = tune_seed, tune_verbose = tune_verbose, ps_lasso_nfolds = ps_lasso_nfolds, ps_lasso_lambda_rule = ps_lasso_lambda_rule, ps_lasso_standardize = ps_lasso_standardize)
    ehat <- ps_obj$ehat
    qhat <- ehat/(1 - ehat)
    if (normalize_controls) {
        n1 <- sum(d$A == 1)
        qhat[d$A == 0] <- n1 * qhat[d$A == 0]/sum(qhat[d$A == 0])
    }
    W <- with(d, A + (1 - A) * qhat)
    fit <- coxph(Surv(time, delta) ~ A, data = d, weights = W, robust = TRUE, ties = "breslow")
    theta_hat <- unname(coef(fit)["A"])
    se_lw <- manual_linwei_se(theta_hat = theta_hat, data = d, W = W)
    idx0 <- which(d$A == 0)
    wd <- .weight_diag(W[idx0])
    list(theta = theta_hat, se_lw = se_lw, W = W, ehat = ehat, qhat = qhat, ps_method = ps_method, ESS = wd$ESS, Rel_ESS = wd$Rel_ESS, W_CV = wd$W_CV, W_Min = wd$W_Min, W_Max = wd$W_Max)
}

cox_lp_basis_label <- function(basis_method = c("cox_lp", "lasso_cox_lp")) {
    basis_method <- match.arg(basis_method)
    if (basis_method == "cox_lp") {
        "Cox LP"
    }
    else {
        "lasso-Cox LP"
    }
}

standardize_lp_basis <- function(lp, eps = 1e-10) {
    lp <- as.numeric(lp)
    if (any(!is.finite(lp))) {
        stop("Cox linear-predictor basis contains non-finite values.")
    }
    s <- stats::sd(lp)
    if (!is.finite(s) || s < eps) {
        stop("Cox linear-predictor basis has near-zero empirical variance.")
    }
    z <- as.numeric((lp - mean(lp))/s)
    H <- cbind(Intercept = 1, Cox_LP = z)
    H
}

predict_cox_linear_predictor_train <- function(train_data, newdata, xvars, cox_iter_max = 50) {
    d_train <- as.data.frame(train_data)
    d_pred <- as.data.frame(newdata)
    d0 <- d_train[d_train$A == 0, , drop = FALSE]
    if (nrow(d0) < 20 || sum(d0$delta == 1) < 5) {
        stop("Too few external-control rows or events for the Cox LP basis.")
    }
    xvars_use <- xvars[vapply(d0[, xvars, drop = FALSE], function(z) {
        stats::sd(z, na.rm = TRUE) > 1e-10
    }, logical(1))]
    if (length(xvars_use) == 0) {
        stop("No nonconstant covariates are available for the Cox LP basis.")
    }
    form <- stats::as.formula(paste0("Surv(time, delta) ~ ", paste(xvars_use, collapse = " + ")))
    fit <- survival::coxph(form, data = d0, ties = "breslow", x = TRUE, singular.ok = TRUE, control = survival::coxph.control(iter.max = cox_iter_max))
    lp <- as.numeric(stats::predict(fit, newdata = d_pred, type = "lp", reference = "zero"))
    if (any(!is.finite(lp))) {
        stop("Ordinary Cox LP produced non-finite predictions. ", "Try mec_basis_lasso <- TRUE for the lasso-Cox LP basis.")
    }
    list(lp = lp, fit = fit, basis_method = "cox_lp", lambda_selected = NA_real_, n_selected = sum(is.finite(stats::coef(fit)) & abs(stats::coef(fit)) > 1e-12))
}

predict_lasso_cox_linear_predictor_train <- function(train_data, newdata, xvars, nfolds = 5, lambda_rule = c("lambda.1se", "lambda.min"), standardize = TRUE, tune_seed = 20260427, tune_verbose = FALSE) {
    lambda_rule <- match.arg(lambda_rule)
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package {glmnet} is required when mec_basis_lasso = TRUE.")
    }
    d_train <- as.data.frame(train_data)
    d_pred <- as.data.frame(newdata)
    d0 <- d_train[d_train$A == 0, , drop = FALSE]
    if (nrow(d0) < 20 || sum(d0$delta == 1) < 5) {
        stop("Too few external-control rows or events for the lasso-Cox LP basis.")
    }
    x0 <- as.matrix(d0[, xvars, drop = FALSE])
    x_new <- as.matrix(d_pred[, xvars, drop = FALSE])
    y0 <- survival::Surv(d0$time, d0$delta)
    nfolds_use <- min(as.integer(nfolds), nrow(d0), sum(d0$delta == 1))
    if (!is.finite(nfolds_use) || nfolds_use < 3) {
        stop("Not enough external-control events for cv.glmnet lasso-Cox tuning.")
    }
    set.seed(tune_seed)
    cvfit <- glmnet::cv.glmnet(x = x0, y = y0, family = "cox", alpha = 1, nfolds = nfolds_use, type.measure = "deviance", standardize = standardize, maxit = 1e+05, parallel = FALSE)
    lambda_selected <- if (lambda_rule == "lambda.min") {
        cvfit$lambda.min
    }
    else {
        cvfit$lambda.1se
    }
    lp <- as.numeric(stats::predict(cvfit, newx = x_new, s = lambda_selected, type = "link"))
    if (any(!is.finite(lp))) {
        stop("Lasso-Cox LP produced non-finite predictions.")
    }
    beta_hat <- as.matrix(stats::coef(cvfit, s = lambda_selected))
    n_selected <- sum(abs(beta_hat[, 1]) > 1e-12)
    if (isTRUE(tune_verbose)) {
        message("lasso-Cox basis selected ", lambda_rule, " = ", signif(lambda_selected, 4), "; nonzero coefficients = ", n_selected, "; nfolds = ", nfolds_use)
    }
    list(lp = lp, fit = cvfit, basis_method = "lasso_cox_lp", lambda_selected = lambda_selected, n_selected = n_selected)
}

predict_mec_cox_lp_basis_train <- function(train_data, newdata, xvars, basis_method = c("cox_lp", "lasso_cox_lp"), lasso_nfolds = 5, lasso_lambda_rule = c("lambda.1se", "lambda.min"), lasso_standardize = TRUE, tune_seed = 20260427, tune_verbose = FALSE) {
    basis_method <- match.arg(basis_method)
    lasso_lambda_rule <- match.arg(lasso_lambda_rule)
    if (basis_method == "cox_lp") {
        predict_cox_linear_predictor_train(train_data = train_data, newdata = newdata, xvars = xvars)
    }
    else {
        predict_lasso_cox_linear_predictor_train(train_data = train_data, newdata = newdata, xvars = xvars, nfolds = lasso_nfolds, lambda_rule = lasso_lambda_rule, standardize = lasso_standardize, tune_seed = tune_seed, tune_verbose = tune_verbose)
    }
}

estimate_mec_nuisance <- function(data, xvars, ps_method = c("glm", "glm_lasso", "nnet", "rf", "bart"), or_method = c("cox_lp", "lasso_cox_lp"), use_crossfit = TRUE, Kfold = 5, trim = c(0.01, 0.99), fold_seed = 20260427, auto_tune = FALSE, tune_n = 100, tune_seed = 20260427, tune_verbose = FALSE, ps_lasso_nfolds = 5, ps_lasso_lambda_rule = c("lambda.1se", "lambda.min"), ps_lasso_standardize = TRUE, lasso_nfolds = 5, lasso_lambda_rule = c("lambda.1se", "lambda.min"), lasso_standardize = TRUE) {
    ps_method <- match.arg(ps_method)
    or_method <- match.arg(or_method)
    ps_lasso_lambda_rule <- match.arg(ps_lasso_lambda_rule)
    lasso_lambda_rule <- match.arg(lasso_lambda_rule)
    d <- as.data.frame(data)
    n <- nrow(d)
    if (!use_crossfit) {
        ps_obj <- estimate_ps_ml(data = d, xvars = xvars, ps_method = ps_method, newdata = d, trim = trim, auto_tune = auto_tune, tune_n = tune_n, tune_seed = tune_seed, tune_verbose = tune_verbose, ps_lasso_nfolds = ps_lasso_nfolds, ps_lasso_lambda_rule = ps_lasso_lambda_rule, ps_lasso_standardize = ps_lasso_standardize)
        pred <- predict_mec_cox_lp_basis_train(train_data = d, newdata = d, xvars = xvars, basis_method = or_method, lasso_nfolds = lasso_nfolds, lasso_lambda_rule = lasso_lambda_rule, lasso_standardize = lasso_standardize, tune_seed = tune_seed, tune_verbose = tune_verbose)
        H <- standardize_lp_basis(pred$lp)
        return(list(ehat = ps_obj$ehat, H = H, basis_method = or_method, basis_label = cox_lp_basis_label(or_method), basis_dim = ncol(H), lambda_selected = pred$lambda_selected, n_selected = pred$n_selected, use_crossfit = FALSE, fold_id = NULL))
    }
    fold_id <- make_stratified_folds(d$A, K = Kfold, seed = fold_seed)
    ehat_cf <- rep(NA_real_, n)
    lp_cf <- rep(NA_real_, n)
    lambda_selected <- rep(NA_real_, Kfold)
    n_selected <- rep(NA_real_, Kfold)
    for (k in seq_len(Kfold)) {
        idx_valid <- which(fold_id == k)
        idx_train <- which(fold_id != k)
        d_train <- d[idx_train, , drop = FALSE]
        d_valid <- d[idx_valid, , drop = FALSE]
        ps_obj_k <- estimate_ps_ml(data = d_train, xvars = xvars, ps_method = ps_method, newdata = d_valid, trim = trim, auto_tune = auto_tune, tune_n = tune_n, tune_seed = tune_seed + k, tune_verbose = tune_verbose, ps_lasso_nfolds = ps_lasso_nfolds, ps_lasso_lambda_rule = ps_lasso_lambda_rule, ps_lasso_standardize = ps_lasso_standardize)
        ehat_cf[idx_valid] <- ps_obj_k$ehat
        pred_k <- predict_mec_cox_lp_basis_train(train_data = d_train, newdata = d_valid, xvars = xvars, basis_method = or_method, lasso_nfolds = lasso_nfolds, lasso_lambda_rule = lasso_lambda_rule, lasso_standardize = lasso_standardize, tune_seed = tune_seed + k, tune_verbose = tune_verbose)
        lp_cf[idx_valid] <- pred_k$lp
        lambda_selected[k] <- pred_k$lambda_selected
        n_selected[k] <- pred_k$n_selected
    }
    if (any(!is.finite(ehat_cf))) {
        stop("Cross-fitted PS contains non-finite values.")
    }
    if (any(!is.finite(lp_cf))) {
        stop("Cross-fitted Cox linear-predictor basis contains non-finite values.")
    }
    H <- standardize_lp_basis(lp_cf)
    list(ehat = ehat_cf, H = H, basis_method = or_method, basis_label = cox_lp_basis_label(or_method), basis_dim = ncol(H), lambda_selected = lambda_selected, n_selected = n_selected, use_crossfit = TRUE, fold_id = fold_id)
}

fit_mec_cox_from_nuisance <- function(data, nuis, divergence = c("kl", "quadratic", "el", "hellinger", "renyi"), normalize_baseline = TRUE, use_crossfit = TRUE, Kfold = 5) {
    divergence <- match.arg(divergence)
    d <- as.data.frame(data)
    n <- nrow(d)
    ehat <- nuis$ehat
    H <- nuis$H
    qhat <- ehat/(1 - ehat)
    idx1 <- which(d$A == 1)
    idx0 <- which(d$A == 0)
    n1 <- length(idx1)
    d0 <- qhat[idx0]
    if (normalize_baseline) {
        d0 <- n1 * d0/sum(d0)
    }
    H1 <- H[idx1, , drop = FALSE]
    H0 <- H[idx0, , drop = FALSE]
    cal <- calibrate_weights_bregman_from_design_return(Xs = H0, Xp = H1, d = d0, divergence = divergence, alpha = 1/2)
    w0_mec <- cal$w
    W_mec <- numeric(n)
    W_mec[idx1] <- 1
    W_mec[idx0] <- w0_mec
    if (any(!is.finite(W_mec)) || any(W_mec <= 0)) {
        stop("MEC calibration produced nonpositive or nonfinite weights for generator = '", divergence, "'. This can happen for the quadratic generator under weak overlap or an aggressive calibration basis.")
    }
    fit <- survival::coxph(survival::Surv(time, delta) ~ A, data = d, weights = W_mec, robust = TRUE, ties = "breslow")
    theta_hat <- unname(stats::coef(fit)["A"])
    ce <- cox_score_eta_fixed_weight(theta = theta_hat, data = d, W = W_mec)
    eta_theta <- ce$eta
    D_theta_theta <- as.numeric(jacobian_fd(function(th) {
        c(cox_score_eta_fixed_weight(th[1], d, W_mec)$U)
    }, x = c(theta_hat)))
    fam <- .bregman_family(divergence, alpha = 1/2)
    gd0 <- fam$g(d0)
    weights_from_lambda <- function(lambda) {
        nu0 <- as.numeric(gd0 + H0 %*% lambda)
        w0 <- fam$ginv(nu0)
        W <- numeric(n)
        W[idx1] <- 1
        W[idx0] <- w0
        W
    }
    D_theta_lambda <- if (divergence == "kl") {
        matrix(as.numeric(crossprod(eta_theta[idx0], H0)), nrow = 1L)
    }
    else {
        jacobian_fd(function(lambda) {
            W_lambda <- weights_from_lambda(lambda)
            c(cox_score_eta_fixed_weight(theta_hat, d, W_lambda)$U)
        }, x = cal$lambda)
    }
    rho <- matrix(0, nrow = n, ncol = ncol(H))
    rho[idx0, ] <- sweep(H0, 1, w0_mec, "*")
    rho[idx1, ] <- -H1
    D_lambda_lambda <- cal$D_lambda_lambda
    K <- ncol(H)
    D <- rbind(c(D_theta_theta, as.numeric(D_theta_lambda)), cbind(rep(0, K), D_lambda_lambda))
    U_i <- cbind(eta_theta, rho)
    IF <- -U_i %*% t(solve_safe(D))
    phi_theta <- IF[, 1]
    se_mec <- sqrt(sum(phi_theta^2))
    if (!is.finite(se_mec)) {
        se_mec <- sqrt(as.numeric(stats::vcov(fit)["A", "A"]))
    }
    wd <- .weight_diag(w0_mec)
    list(theta = theta_hat, se = se_mec, W = W_mec, ehat = ehat, qhat = qhat, H = H, lambda = cal$lambda, calibration_grad_norm = cal$grad_norm, calibration_converged = cal$converged, basis_method = nuis$basis_method, basis_label = nuis$basis_label, basis_dim = nuis$basis_dim, lasso_lambda_selected = nuis$lambda_selected, lasso_n_selected = nuis$n_selected, use_crossfit = use_crossfit, Kfold = Kfold, generator = divergence, ESS = wd$ESS, Rel_ESS = wd$Rel_ESS, W_CV = wd$W_CV, W_Min = wd$W_Min, W_Max = wd$W_Max)
}

fit_mec_cox <- function(data, xvars, ps_method = c("glm", "glm_lasso", "nnet", "rf", "bart"), or_method = c("cox_lp", "lasso_cox_lp"), use_crossfit = TRUE, Kfold = 5, divergence = c("kl", "quadratic", "el", "hellinger", "renyi"), trim = c(0.01, 0.99), normalize_baseline = TRUE, auto_tune = FALSE, tune_n = 100, tune_seed = 20260427, tune_verbose = FALSE, ps_lasso_nfolds = 5, ps_lasso_lambda_rule = c("lambda.1se", "lambda.min"), ps_lasso_standardize = TRUE, lasso_nfolds = 5, lasso_lambda_rule = c("lambda.1se", 
    "lambda.min"), lasso_standardize = TRUE) {
    ps_method <- match.arg(ps_method)
    or_method <- match.arg(or_method)
    ps_lasso_lambda_rule <- match.arg(ps_lasso_lambda_rule)
    divergence <- match.arg(divergence)
    lasso_lambda_rule <- match.arg(lasso_lambda_rule)
    d <- as.data.frame(data)
    nuis <- estimate_mec_nuisance(data = d, xvars = xvars, ps_method = ps_method, or_method = or_method, use_crossfit = use_crossfit, Kfold = Kfold, trim = trim, auto_tune = auto_tune, tune_n = tune_n, tune_seed = tune_seed, tune_verbose = tune_verbose, ps_lasso_nfolds = ps_lasso_nfolds, ps_lasso_lambda_rule = ps_lasso_lambda_rule, ps_lasso_standardize = ps_lasso_standardize, lasso_nfolds = lasso_nfolds, lasso_lambda_rule = lasso_lambda_rule, lasso_standardize = lasso_standardize)
    fit_mec_cox_from_nuisance(data = d, nuis = nuis, divergence = divergence, normalize_baseline = normalize_baseline, use_crossfit = use_crossfit, Kfold = Kfold)
}

compute_true_att_target <- function(n1_super = 30000, n0_super = 60000, p = 10, ps_truth = "linear", outcome_truth = "linear_ph", ps_nonlinearity = 1, or_nonlinearity = 1, beta_cond = log(0.7), xi = 0, h_type = "constant", lambda0 = 8e-05, eta_shape = 2, censor_rate = 8e-04, seed = 12345) {
    set.seed(seed)
    dat <- generate_observed_data(n1 = n1_super, n0 = n0_super, p = p, ps_truth = ps_truth, outcome_truth = outcome_truth, ps_nonlinearity = ps_nonlinearity, or_nonlinearity = or_nonlinearity, beta_cond = beta_cond, xi = xi, h_type = h_type, lambda0 = lambda0, eta_shape = eta_shape, censor_rate = censor_rate)
    xvars <- paste0("X", seq_len(p))
    X <- as.matrix(dat[, xvars, drop = FALSE])
    e_true <- source_score(X, ps_truth = ps_truth, ps_nonlinearity = ps_nonlinearity)
    q_true <- e_true/(1 - e_true)
    idx0 <- which(dat$A == 0)
    idx1 <- which(dat$A == 1)
    q_true[idx0] <- length(idx1) * q_true[idx0]/sum(q_true[idx0])
    W_true <- with(dat, A + (1 - A) * q_true)
    fit <- coxph(Surv(time, delta) ~ A, data = dat, weights = W_true, ties = "breslow")
    unname(coef(fit)["A"])
}

empty_row <- function(rep_id, n1, n0, method, theta_true, error_msg = NA_character_) {
    data.frame(rep = rep_id, n1 = n1, n0 = n0, Method = method, theta_true = theta_true, Estimate = NA_real_, SE = NA_real_, CI_L = NA_real_, CI_U = NA_real_, ESS = NA_real_, Rel_ESS = NA_real_, W_CV = NA_real_, W_Min = NA_real_, W_Max = NA_real_, Cal_Grad = NA_real_, Cal_Converged = NA, Error = error_msg, stringsAsFactors = FALSE)
}

fit_one_rep <- function(rep_id, n1, n0, theta_true, p = 10, ps_truth = "linear", outcome_truth = "linear_ph", ps_nonlinearity = 1, or_nonlinearity = 1, beta_cond = log(0.7), xi = 0, h_type = "constant", lambda0 = 8e-05, eta_shape = 2, censor_rate = 8e-04, mec_ps_method = "glm", mec_or_method = "cox_lp", Kfold_mec = 5, mec_generators = c("quadratic", "kl", "el", "hellinger", "renyi"), trim = c(0.01, 0.99), mec_lasso_nfolds = 5, mec_lasso_lambda_rule = "lambda.1se", mec_lasso_standardize = TRUE, include_mec_no_cf = FALSE, 
    include_unweighted_cox = TRUE, include_linwei_ml = TRUE, Lin_Wei_ps_method = "glm", auto_tune_ml = FALSE, auto_tune_n = 100, auto_tune_seed = 20260427, auto_tune_verbose = FALSE, ps_lasso_nfolds = 5, ps_lasso_lambda_rule = c("lambda.1se", "lambda.min"), ps_lasso_standardize = TRUE) {
    ps_lasso_lambda_rule <- match.arg(ps_lasso_lambda_rule)
    mec_generators <- match.arg(mec_generators, choices = c("quadratic", "kl", "el", "hellinger", "renyi"), several.ok = TRUE)
    mec_generators <- unique(mec_generators)
    dat <- generate_observed_data(n1 = n1, n0 = n0, p = p, ps_truth = ps_truth, outcome_truth = outcome_truth, ps_nonlinearity = ps_nonlinearity, or_nonlinearity = or_nonlinearity, beta_cond = beta_cond, xi = xi, h_type = h_type, lambda0 = lambda0, eta_shape = eta_shape, censor_rate = censor_rate)
    xvars <- paste0("X", seq_len(p))
    rows <- list()
    if (isTRUE(include_unweighted_cox)) {
        ans_unweighted <- tryCatch(fit_unweighted_cox(dat), error = function(e) e)
        if (inherits(ans_unweighted, "error")) {
            rows[[length(rows) + 1]] <- empty_row(rep_id, n1, n0, "Unweighted Cox", theta_true, conditionMessage(ans_unweighted))
        }
        else {
            rows[[length(rows) + 1]] <- data.frame(rep = rep_id, n1 = n1, n0 = n0, Method = "Unweighted Cox", theta_true = theta_true, Estimate = ans_unweighted$theta, SE = ans_unweighted$se, CI_L = ans_unweighted$theta - 1.96 * ans_unweighted$se, CI_U = ans_unweighted$theta + 1.96 * ans_unweighted$se, ESS = NA_real_, Rel_ESS = NA_real_, W_CV = NA_real_, W_Min = NA_real_, W_Max = NA_real_, Cal_Grad = NA_real_, Cal_Converged = NA, Error = NA_character_, stringsAsFactors = FALSE)
        }
    }
    ans_ipw <- tryCatch(fit_att_ipw_naive_linwei_shu(data = dat, xvars = xvars, trim = trim, normalize_controls = TRUE), error = function(e) e)
    if (inherits(ans_ipw, "error")) {
        rows[[length(rows) + 1]] <- empty_row(rep_id, n1, n0, "ATT-IPW Cox: Naive model-based", theta_true, conditionMessage(ans_ipw))
        rows[[length(rows) + 1]] <- empty_row(rep_id, n1, n0, "ATT-IPW Cox: Lin-Wei", theta_true, conditionMessage(ans_ipw))
        rows[[length(rows) + 1]] <- empty_row(rep_id, n1, n0, "ATT-IPW Cox: Shu", theta_true, conditionMessage(ans_ipw))
    }
    else {
        rows[[length(rows) + 1]] <- data.frame(rep = rep_id, n1 = n1, n0 = n0, Method = "ATT-IPW Cox: Naive model-based", theta_true = theta_true, Estimate = ans_ipw$theta, SE = ans_ipw$se_naive, CI_L = ans_ipw$theta - 1.96 * ans_ipw$se_naive, CI_U = ans_ipw$theta + 1.96 * ans_ipw$se_naive, ESS = NA_real_, Rel_ESS = NA_real_, W_CV = NA_real_, W_Min = NA_real_, W_Max = NA_real_, Cal_Grad = NA_real_, Cal_Converged = NA, Error = NA_character_, stringsAsFactors = FALSE)
        rows[[length(rows) + 1]] <- data.frame(rep = rep_id, n1 = n1, n0 = n0, Method = "ATT-IPW Cox: Lin-Wei", theta_true = theta_true, Estimate = ans_ipw$theta, SE = ans_ipw$se_lw, CI_L = ans_ipw$theta - 1.96 * ans_ipw$se_lw, CI_U = ans_ipw$theta + 1.96 * ans_ipw$se_lw, ESS = NA_real_, Rel_ESS = NA_real_, W_CV = NA_real_, W_Min = NA_real_, W_Max = NA_real_, Cal_Grad = NA_real_, Cal_Converged = NA, Error = NA_character_, stringsAsFactors = FALSE)
        rows[[length(rows) + 1]] <- data.frame(rep = rep_id, n1 = n1, n0 = n0, Method = "ATT-IPW Cox: Shu", theta_true = theta_true, Estimate = ans_ipw$theta, SE = ans_ipw$se_shu, CI_L = ans_ipw$theta - 1.96 * ans_ipw$se_shu, CI_U = ans_ipw$theta + 1.96 * ans_ipw$se_shu, ESS = NA_real_, Rel_ESS = NA_real_, W_CV = NA_real_, W_Min = NA_real_, W_Max = NA_real_, Cal_Grad = NA_real_, Cal_Converged = NA, Error = NA_character_, stringsAsFactors = FALSE)
    }
    if (isTRUE(include_linwei_ml)) {
        linwei_ml_label <- paste0("ATT-IPW Cox: Lin-Wei (PS=", Lin_Wei_ps_method, ", no CF)")
        ans_lw_ml <- tryCatch(fit_att_ipw_linwei_ml(data = dat, xvars = xvars, ps_method = Lin_Wei_ps_method, trim = trim, normalize_controls = TRUE, auto_tune = auto_tune_ml, tune_n = auto_tune_n, tune_seed = auto_tune_seed + rep_id, tune_verbose = auto_tune_verbose, ps_lasso_nfolds = ps_lasso_nfolds, ps_lasso_lambda_rule = ps_lasso_lambda_rule, ps_lasso_standardize = ps_lasso_standardize), error = function(e) e)
        if (inherits(ans_lw_ml, "error")) {
            rows[[length(rows) + 1]] <- empty_row(rep_id, n1, n0, linwei_ml_label, theta_true, conditionMessage(ans_lw_ml))
        }
        else {
            rows[[length(rows) + 1]] <- data.frame(rep = rep_id, n1 = n1, n0 = n0, Method = linwei_ml_label, theta_true = theta_true, Estimate = ans_lw_ml$theta, SE = ans_lw_ml$se_lw, CI_L = ans_lw_ml$theta - 1.96 * ans_lw_ml$se_lw, CI_U = ans_lw_ml$theta + 1.96 * ans_lw_ml$se_lw, ESS = ans_lw_ml$ESS, Rel_ESS = ans_lw_ml$Rel_ESS, W_CV = ans_lw_ml$W_CV, W_Min = ans_lw_ml$W_Min, W_Max = ans_lw_ml$W_Max, Cal_Grad = NA_real_, Cal_Converged = NA, Error = NA_character_, stringsAsFactors = FALSE)
        }
    }
    mec_basis_label <- cox_lp_basis_label(mec_or_method)
    nuis_mec_cf <- tryCatch(estimate_mec_nuisance(data = dat, xvars = xvars, ps_method = mec_ps_method, or_method = mec_or_method, use_crossfit = TRUE, Kfold = Kfold_mec, trim = trim, auto_tune = auto_tune_ml, tune_n = auto_tune_n, tune_seed = auto_tune_seed + rep_id, tune_verbose = auto_tune_verbose, ps_lasso_nfolds = ps_lasso_nfolds, ps_lasso_lambda_rule = ps_lasso_lambda_rule, ps_lasso_standardize = ps_lasso_standardize, lasso_nfolds = mec_lasso_nfolds, lasso_lambda_rule = mec_lasso_lambda_rule, 
        lasso_standardize = mec_lasso_standardize), error = function(e) e)
    for (mec_generator in mec_generators) {
        mec_label_cf <- paste0("MEC-Cox: PS=", mec_ps_method, ", basis=", mec_basis_label, ", dim=2", ", K=", Kfold_mec, ", Generator=", mec_generator)
        if (inherits(nuis_mec_cf, "error")) {
            rows[[length(rows) + 1]] <- empty_row(rep_id, n1, n0, mec_label_cf, theta_true, conditionMessage(nuis_mec_cf))
            next
        }
        ans_mec_cf <- tryCatch(fit_mec_cox_from_nuisance(data = dat, nuis = nuis_mec_cf, divergence = mec_generator, normalize_baseline = TRUE, use_crossfit = TRUE, Kfold = Kfold_mec), error = function(e) e)
        if (inherits(ans_mec_cf, "error")) {
            rows[[length(rows) + 1]] <- empty_row(rep_id, n1, n0, mec_label_cf, theta_true, conditionMessage(ans_mec_cf))
        }
        else {
            mec_label_cf <- paste0("MEC-Cox: PS=", mec_ps_method, ", basis=", ans_mec_cf$basis_label, ", dim=", ncol(ans_mec_cf$H), ", K=", Kfold_mec, ", Generator=", mec_generator)
            rows[[length(rows) + 1]] <- data.frame(rep = rep_id, n1 = n1, n0 = n0, Method = mec_label_cf, theta_true = theta_true, Estimate = ans_mec_cf$theta, SE = ans_mec_cf$se, CI_L = ans_mec_cf$theta - 1.96 * ans_mec_cf$se, CI_U = ans_mec_cf$theta + 1.96 * ans_mec_cf$se, ESS = ans_mec_cf$ESS, Rel_ESS = ans_mec_cf$Rel_ESS, W_CV = ans_mec_cf$W_CV, W_Min = ans_mec_cf$W_Min, W_Max = ans_mec_cf$W_Max, Cal_Grad = ans_mec_cf$calibration_grad_norm, Cal_Converged = ans_mec_cf$calibration_converged, Error = NA_character_, 
                stringsAsFactors = FALSE)
        }
    }
    if (isTRUE(include_mec_no_cf)) {
        nuis_mec_nocf <- tryCatch(estimate_mec_nuisance(data = dat, xvars = xvars, ps_method = mec_ps_method, or_method = mec_or_method, use_crossfit = FALSE, Kfold = Kfold_mec, trim = trim, auto_tune = auto_tune_ml, tune_n = auto_tune_n, tune_seed = auto_tune_seed + rep_id, tune_verbose = auto_tune_verbose, ps_lasso_nfolds = ps_lasso_nfolds, ps_lasso_lambda_rule = ps_lasso_lambda_rule, ps_lasso_standardize = ps_lasso_standardize, lasso_nfolds = mec_lasso_nfolds, lasso_lambda_rule = mec_lasso_lambda_rule, 
            lasso_standardize = mec_lasso_standardize), error = function(e) e)
        for (mec_generator in mec_generators) {
            mec_label_nocf <- paste0("MEC-Cox (no CF): PS=", mec_ps_method, ", basis=", mec_basis_label, ", dim=2", ", Generator=", mec_generator)
            if (inherits(nuis_mec_nocf, "error")) {
                rows[[length(rows) + 1]] <- empty_row(rep_id, n1, n0, mec_label_nocf, theta_true, conditionMessage(nuis_mec_nocf))
                next
            }
            ans_mec_nocf <- tryCatch(fit_mec_cox_from_nuisance(data = dat, nuis = nuis_mec_nocf, divergence = mec_generator, normalize_baseline = TRUE, use_crossfit = FALSE, Kfold = Kfold_mec), error = function(e) e)
            if (inherits(ans_mec_nocf, "error")) {
                rows[[length(rows) + 1]] <- empty_row(rep_id, n1, n0, mec_label_nocf, theta_true, conditionMessage(ans_mec_nocf))
            }
            else {
                mec_label_nocf <- paste0("MEC-Cox (no CF): PS=", mec_ps_method, ", basis=", ans_mec_nocf$basis_label, ", dim=", ncol(ans_mec_nocf$H), ", Generator=", mec_generator)
                rows[[length(rows) + 1]] <- data.frame(rep = rep_id, n1 = n1, n0 = n0, Method = mec_label_nocf, theta_true = theta_true, Estimate = ans_mec_nocf$theta, SE = ans_mec_nocf$se, CI_L = ans_mec_nocf$theta - 1.96 * ans_mec_nocf$se, CI_U = ans_mec_nocf$theta + 1.96 * ans_mec_nocf$se, ESS = ans_mec_nocf$ESS, Rel_ESS = ans_mec_nocf$Rel_ESS, W_CV = ans_mec_nocf$W_CV, W_Min = ans_mec_nocf$W_Min, W_Max = ans_mec_nocf$W_Max, Cal_Grad = ans_mec_nocf$calibration_grad_norm, Cal_Converged = ans_mec_nocf$calibration_converged, 
                  Error = NA_character_, stringsAsFactors = FALSE)
            }
        }
    }
    do.call(rbind, rows)
}

required_worker_packages <- function(mec_ps_method, mec_or_method, Lin_Wei_ps_method = "glm") {
    pkgs <- c("survival", "MASS", "dplyr")
    ps_methods_needed <- unique(c(mec_ps_method, Lin_Wei_ps_method))
    if ("bart" %in% ps_methods_needed) 
        pkgs <- c(pkgs, "dbarts")
    if ("nnet" %in% ps_methods_needed) 
        pkgs <- c(pkgs, "nnet")
    if ("rf" %in% ps_methods_needed) 
        pkgs <- c(pkgs, "ranger")
    if ("glm_lasso" %in% ps_methods_needed) 
        pkgs <- c(pkgs, "glmnet")
    if (mec_or_method == "lasso_cox_lp") 
        pkgs <- c(pkgs, "glmnet")
    unique(pkgs)
}

run_simulation_parallel <- function(sample_size_grid, R, theta_true, p, ps_truth, outcome_truth, ps_nonlinearity = 1, or_nonlinearity = 1, beta_cond, xi, h_type, lambda0, eta_shape, censor_rate, mec_ps_method, mec_or_method, Kfold_mec, mec_generators, trim, mec_lasso_nfolds = 5, mec_lasso_lambda_rule = "lambda.1se", mec_lasso_standardize = TRUE, Lin_Wei_ps_method = "glm", include_mec_no_cf = FALSE, include_unweighted_cox = TRUE, include_linwei_ml = TRUE, auto_tune_ml = FALSE, auto_tune_n = 100, auto_tune_seed = 20260427, 
    auto_tune_verbose = FALSE, ps_lasso_nfolds = 5, ps_lasso_lambda_rule = c("lambda.1se", "lambda.min"), ps_lasso_standardize = TRUE, cores = max(1, min(8, parallel::detectCores() - 1)), rng_seed = 20260427) {
    ps_lasso_lambda_rule <- match.arg(ps_lasso_lambda_rule)
    mec_generators <- match.arg(mec_generators, choices = c("quadratic", "kl", "el", "hellinger", "renyi"), several.ok = TRUE)
    mec_generators <- unique(mec_generators)
    worker_pkgs <- required_worker_packages(mec_ps_method = mec_ps_method, mec_or_method = mec_or_method, Lin_Wei_ps_method = Lin_Wei_ps_method)
    missing_pkgs <- worker_pkgs[!vapply(worker_pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), ". Install them first.")
    }
    job_grid <- expand.grid(ss = seq_len(nrow(sample_size_grid)), rep_id = seq_len(R))
    job_grid$n1 <- sample_size_grid$n1[job_grid$ss]
    job_grid$n0 <- sample_size_grid$n0[job_grid$ss]
    n_jobs <- nrow(job_grid)
    cl <- parallel::makeCluster(cores, type = "SOCK")
    on.exit({
        try(parallel::stopCluster(cl), silent = TRUE)
    }, add = TRUE)
    doSNOW::registerDoSNOW(cl)
    doRNG::registerDoRNG(rng_seed)
    parallel::clusterExport(cl, varlist = c("worker_pkgs"), envir = environment())
    parallel::clusterEvalQ(cl, {
        for (pkg in worker_pkgs) {
            suppressPackageStartupMessages(library(pkg, character.only = TRUE))
        }
        NULL
    })
    export_names <- c("expit", "clip01", "round_numeric_df", "solve_safe", "jacobian_fd", "make_stratified_folds", "make_balanced_tuning_indices", "make_stratified_train_valid_split", "binary_logloss", "safe_metric_min", "safe_metric_max", "tune_ps_ml", "source_score", "sample_source_covariates", "domain_shift", "prognostic_lp", "weibull_time", "generate_observed_data", "make_att_ipw_data_glm", "att_weights_from_gamma", "cox_score_eta_fixed_weight", "manual_linwei_se", "fit_unweighted_cox", "fit_att_ipw_naive_linwei_shu", 
        "fit_att_ipw_linwei_ml", ".to_matrix", ".bregman_family", ".weight_diag", "calibrate_weights_bregman_from_design_return", "estimate_ps_ml", "cox_lp_basis_label", "standardize_lp_basis", "predict_cox_linear_predictor_train", "predict_lasso_cox_linear_predictor_train", "predict_mec_cox_lp_basis_train", "estimate_mec_nuisance", "fit_mec_cox_from_nuisance", "fit_mec_cox", "empty_row", "fit_one_rep")
    parallel::clusterExport(cl, varlist = export_names, envir = globalenv())
    pb <- txtProgressBar(min = 0, max = n_jobs, style = 3)
    opts <- list(progress = function(n) setTxtProgressBar(pb, n))
    schema <- data.frame(rep = integer(1), n1 = integer(1), n0 = integer(1), Method = character(1), theta_true = numeric(1), Estimate = numeric(1), SE = numeric(1), CI_L = numeric(1), CI_U = numeric(1), ESS = numeric(1), Rel_ESS = numeric(1), W_CV = numeric(1), W_Min = numeric(1), W_Max = numeric(1), Cal_Grad = numeric(1), Cal_Converged = logical(1), Error = character(1), stringsAsFactors = FALSE)
    res <- foreach(jj = seq_len(n_jobs), .combine = "rbind", .inorder = FALSE, .options.snow = opts, .errorhandling = "pass", .packages = worker_pkgs) %dopar% {
        n1_now <- job_grid$n1[jj]
        n0_now <- job_grid$n0[jj]
        rep_now <- job_grid$rep_id[jj]
        ok <- tryCatch(fit_one_rep(rep_id = rep_now, n1 = n1_now, n0 = n0_now, theta_true = theta_true, p = p, ps_truth = ps_truth, outcome_truth = outcome_truth, ps_nonlinearity = ps_nonlinearity, or_nonlinearity = or_nonlinearity, beta_cond = beta_cond, xi = xi, h_type = h_type, lambda0 = lambda0, eta_shape = eta_shape, censor_rate = censor_rate, mec_ps_method = mec_ps_method, mec_or_method = mec_or_method, Kfold_mec = Kfold_mec, mec_generators = mec_generators, trim = trim, mec_lasso_nfolds = mec_lasso_nfolds, 
            mec_lasso_lambda_rule = mec_lasso_lambda_rule, mec_lasso_standardize = mec_lasso_standardize, Lin_Wei_ps_method = Lin_Wei_ps_method, include_mec_no_cf = include_mec_no_cf, include_unweighted_cox = include_unweighted_cox, include_linwei_ml = include_linwei_ml, auto_tune_ml = auto_tune_ml, auto_tune_n = auto_tune_n, auto_tune_seed = auto_tune_seed + rep_now, auto_tune_verbose = auto_tune_verbose, ps_lasso_nfolds = ps_lasso_nfolds, ps_lasso_lambda_rule = ps_lasso_lambda_rule, ps_lasso_standardize = ps_lasso_standardize), 
            error = function(e) e)
        if (inherits(ok, "error")) {
            er <- schema
            er$rep <- rep_now
            er$n1 <- n1_now
            er$n0 <- n0_now
            er$Method <- "ERROR"
            er$theta_true <- theta_true
            er$Error <- paste0(class(ok)[1], ": ", conditionMessage(ok))
            er
        }
        else {
            ok[, names(schema)]
        }
    }
    close(pb)
    res
}

is_mec_method <- function(x) {
    grepl("^MEC-Cox", x)
}

validate_mec_generators <- function(mec_generators = c("kl")) {
    allowed <- c("quadratic", "kl", "el", "hellinger", "renyi")
    if (is.null(mec_generators) || length(mec_generators) == 0) {
        stop("mec_generators must contain at least one generator.")
    }
    mec_generators <- unique(tolower(as.character(mec_generators)))
    bad <- setdiff(mec_generators, allowed)
    if (length(bad) > 0) {
        stop("Unknown MEC generator(s): ", paste(bad, collapse = ", "), ". Allowed choices are: ", paste(allowed, collapse = ", "), ".")
    }
    mec_generators
}

format_mec_generators <- function(mec_generators) {
    paste(validate_mec_generators(mec_generators), collapse = ", ")
}

extract_mec_generator <- function(x) {
    x <- as.character(x)
    m <- regexec("Generator=([^,]+)", x, perl = TRUE)
    z <- regmatches(x, m)
    vapply(z, function(xx) {
        if (length(xx) >= 2) {
            xx[2]
        }
        else {
            NA_character_
        }
    }, character(1))
}

pretty_generator_label <- function(generator) {
    generator <- tolower(as.character(generator))
    dplyr::case_when(generator == "quadratic" ~ "Quad.", generator == "kl" ~ "KL", generator == "el" ~ "EL", generator == "hellinger" ~ "Hell.", generator == "renyi" ~ "Renyi", TRUE ~ generator)
}

fmt_num <- function(x, digits = 4) {
    ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

summarize_simulation <- function(res) {
    weight_cols <- c("Avg_ESS", "Avg_Rel_ESS", "Avg_W_CV", "Avg_W_Min", "Avg_W_Max", "Avg_Cal_Grad", "Cal_Conv_Rate")
    res %>% filter(is.finite(Estimate), is.finite(SE)) %>% group_by(n1, n0, Method) %>% summarise(theta_true = mean(theta_true), Mean_Est = mean(Estimate), Bias = mean(Estimate - theta_true), Emp_SD = sd(Estimate), Avg_SE = mean(SE), MSE = mean((Estimate - theta_true)^2), RMSE = sqrt(MSE), Mean_CI_Len = mean(CI_U - CI_L), Coverage_95 = mean(theta_true >= CI_L & theta_true <= CI_U), Avg_ESS = mean(ESS, na.rm = TRUE), Avg_Rel_ESS = mean(Rel_ESS, na.rm = TRUE), Avg_W_CV = mean(W_CV, na.rm = TRUE), Avg_W_Min = mean(W_Min, 
        na.rm = TRUE), Avg_W_Max = mean(W_Max, na.rm = TRUE), Avg_Cal_Grad = mean(Cal_Grad, na.rm = TRUE), Cal_Conv_Rate = mean(Cal_Converged, na.rm = TRUE), Trials = dplyr::n(), .groups = "drop") %>% mutate(across(all_of(weight_cols), ~ifelse(is.nan(.x), NA_real_, .x)), Is_MEC = is_mec_method(Method)) %>% mutate(across(all_of(weight_cols), ~ifelse(Is_MEC, .x, NA_real_))) %>% select(-Is_MEC)
}

check_mec_success <- function(res) {
    mec_rows <- res %>% filter(grepl("^MEC-Cox", Method))
    if (nrow(mec_rows) == 0) {
        stop("No MEC-Cox rows were created.")
    }
    n_success <- sum(is.finite(mec_rows$Estimate) & is.finite(mec_rows$SE))
    if (n_success == 0) {
        cat("\nMEC-Cox failed in all replicates. Most common errors:\n")
        print(mec_rows %>% count(Method, Error, sort = TRUE) %>% head(10))
        stop("MEC-Cox failed in all replicates. See error table above.")
    }
    invisible(TRUE)
}

select_best_methods <- function(summ, cov_band = c(0.925, 0.975)) {
    summ %>% group_by(n1, n0) %>% group_modify(~{
        dd <- .x %>% filter(Method != "Unweighted Cox", is.finite(Coverage_95), is.finite(RMSE)) %>% mutate(cov_gap = abs(Coverage_95 - 0.95), cov_ok = Coverage_95 >= cov_band[1] & Coverage_95 <= cov_band[2])
        if (nrow(dd) == 0) {
            return(tibble(Method = character(0)))
        }
        if (any(dd$cov_ok)) {
            best <- dd %>% filter(cov_ok) %>% arrange(RMSE, cov_gap) %>% slice(1)
        }
        else {
            best <- dd %>% arrange(cov_gap, RMSE) %>% slice(1)
        }
        tibble(Method = best$Method)
    }) %>% ungroup() %>% mutate(best_key = paste(n1, n0, Method, sep = "||")) %>% select(best_key)
}

make_method_levels <- function(res, Lin_Wei_ps_method, include_unweighted_cox = TRUE, include_mec_no_cf = FALSE) {
    method_levels <- c(if (isTRUE(include_unweighted_cox)) "Unweighted Cox" else character(0), "ATT-IPW Cox: Naive model-based", "ATT-IPW Cox: Lin-Wei", paste0("ATT-IPW Cox: Lin-Wei (PS=", Lin_Wei_ps_method, ", no CF)"), "ATT-IPW Cox: Shu", res %>% filter(grepl("^MEC-Cox:", Method)) %>% pull(Method) %>% unique(), if (isTRUE(include_mec_no_cf)) {
        res %>% filter(grepl("^MEC-Cox \\(no CF\\):", Method)) %>% pull(Method) %>% unique()
    } else {
        character(0)
    })
    method_levels <- unique(method_levels)
    method_levels <- method_levels[method_levels %in% res$Method]
    method_levels
}

make_method_label_map <- function(res, Lin_Wei_ps_method, include_unweighted_cox = TRUE, include_mec_no_cf = FALSE) {
    method_levels <- make_method_levels(res = res, Lin_Wei_ps_method = Lin_Wei_ps_method, include_unweighted_cox = include_unweighted_cox, include_mec_no_cf = include_mec_no_cf)
    generator <- extract_mec_generator(method_levels)
    gen_short <- pretty_generator_label(generator)
    method_short <- dplyr::case_when(method_levels == "Unweighted Cox" ~ "Unweighted", method_levels == "ATT-IPW Cox: Naive model-based" ~ "IPW-Naive", method_levels == "ATT-IPW Cox: Lin-Wei" ~ "IPW-LW", grepl("^ATT-IPW Cox: Lin-Wei \\(PS=", method_levels) ~ "ML-IPW-LW", method_levels == "ATT-IPW Cox: Shu" ~ "IPW-Shu", grepl("^MEC-Cox \\(no CF\\):", method_levels) ~ paste0("MEC-", gen_short, " (no CF)"), grepl("^MEC-Cox:", method_levels) ~ paste0("MEC-", gen_short), TRUE ~ method_levels)
    tibble(Method = method_levels, Method_Short = method_short, Method_Order = seq_along(method_levels))
}

make_summary_kable <- function(summ, ps_truth, outcome_truth, ps_nonlinearity, or_nonlinearity, mec_ps_method, mec_or_method, Kfold_mec, mec_generators, mec_lasso_nfolds = 5, mec_lasso_lambda_rule = "lambda.1se", mec_lasso_standardize = TRUE, cores, Lin_Wei_ps_method, include_linwei_ml, include_mec_no_cf, include_unweighted_cox = TRUE, n0_over_n1, auto_tune_ml = FALSE, auto_tune_n = 100, best_cov_band = c(0.925, 0.975), overcov_red_threshold = 0.975) {
    mec_generator_txt <- format_mec_generators(mec_generators)
    best_keys <- select_best_methods(summ = summ, cov_band = best_cov_band) %>% pull(best_key)
    method_map <- make_method_label_map(res = summ, Lin_Wei_ps_method = Lin_Wei_ps_method, include_unweighted_cox = include_unweighted_cox, include_mec_no_cf = include_mec_no_cf)
    tab <- summ %>% left_join(method_map, by = "Method") %>% mutate(N_total = n1 + n0, Method_Order = ifelse(is.na(Method_Order), 9999, Method_Order), Method_Short = ifelse(is.na(Method_Short), Method, Method_Short), best_key = paste(n1, n0, Method, sep = "||"), Is_Best = best_key %in% best_keys, Is_Overcovered = Coverage_95 > overcov_red_threshold, Is_MEC = is_mec_method(Method)) %>% arrange(N_total, n1, n0, Method_Order, Method)
    group_info <- tab %>% count(n1, n0, N_total, name = "n_rows") %>% arrange(N_total, n1, n0) %>% mutate(start_row = cumsum(lag(n_rows, default = 0)) + 1, end_row = cumsum(n_rows), group_label = paste0("n1 = ", n1, ", n0 = ", n0, ", N = ", N_total))
    tab_fmt <- tab %>% transmute(n1 = n1, n0 = n0, Method = Method_Short, theta_true = fmt_num(theta_true, 4), Mean_Est = fmt_num(Mean_Est, 4), Bias = fmt_num(Bias, 4), Emp_SD = fmt_num(Emp_SD, 4), Avg_SE = fmt_num(Avg_SE, 4), MSE = fmt_num(MSE, 4), RMSE = fmt_num(RMSE, 4), Mean_CI_Len = fmt_num(Mean_CI_Len, 4), Coverage_95 = ifelse(Is_Overcovered, cell_spec(scales::percent(Coverage_95, accuracy = 0.1), color = "red", bold = TRUE), scales::percent(Coverage_95, accuracy = 0.1)), Avg_ESS = fmt_num(Avg_ESS, 
        4), Avg_Rel_ESS = fmt_num(Avg_Rel_ESS, 4), Avg_W_CV = fmt_num(Avg_W_CV, 4), Avg_W_Min = fmt_num(Avg_W_Min, 4), Avg_W_Max = fmt_num(Avg_W_Max, 4), Avg_Cal_Grad = fmt_num(Avg_Cal_Grad, 6), Cal_Conv_Rate = fmt_num(Cal_Conv_Rate, 4), Trials = Trials)
    caption_txt <- paste0("Simulation summary. ", "PS truth = ", ps_truth, ifelse(ps_truth == "nonlinear", paste0(" (degree = ", ps_nonlinearity, ")"), ""), "; outcome truth = ", outcome_truth, ifelse(outcome_truth == "nonlinear_ph", paste0(" (degree = ", or_nonlinearity, ")"), ""), "; p = ", p, "; domain shift: xi = ", xi, ", h_type = ", h_type, "; MEC PS = ", mec_ps_method, "; MEC basis = ", cox_lp_basis_label(mec_or_method), ifelse(mec_or_method == "lasso_cox_lp", paste0(" (", mec_lasso_lambda_rule, 
        ", nfolds = ", mec_lasso_nfolds, ")"), ""), "; basis dimension = 2", "; Kfold MEC = ", Kfold_mec, "; MEC generators = ", mec_generator_txt, "; auto-tune ML = ", auto_tune_ml, ifelse(isTRUE(auto_tune_ml), paste0(" (tune_n = ", auto_tune_n, ")"), ""), "; n1:n0 = 1:", n0_over_n1, "; cores = ", cores, ".")
    note_txt <- paste0("MEC-Cox denotes cross-fitted MEC-Cox by default. ", "MEC-Cox (no CF) is included only when include_mec_no_cf = ", include_mec_no_cf, ". ", "The additional ML-IPW-LW method uses the user-selected PS learner, PS = ", Lin_Wei_ps_method, ", without cross-fitting and treats the estimated ATT weights as fixed. ", "Standard IPW-Naive, IPW-LW, and IPW-Shu use the same logistic-GLM ATT-IPW point estimator and differ only in variance estimation. ", "Unweighted Cox is included only when include_unweighted_cox = ", 
        include_unweighted_cox, ". ", "MEC-Cox uses the two-dimensional basis h_i = (1, Cox_LP_i)^T. If mec_basis_lasso = TRUE, Cox_LP_i is obtained from lasso-Cox with cv.glmnet and the selected lambda rule shown in the caption. For generator sensitivity, the same estimated PS and Cox-LP basis are reused across selected generators within each replicate. ", "Weight-stability diagnostics, including ESS, relative ESS, weight CV, weight range, calibration gradient, and calibration convergence rate, are displayed only for MEC-Cox rows; non-MEC rows are left blank. ", 
        "When auto_tune_ml = TRUE, nnet/rf/bart PS learners are tuned using a small stratified tuning subset of size approximately auto_tune_n; glm and glm_lasso PS use their own model specifications/CV. ", "Best method is highlighted in green within each sample-size group, excluding Unweighted Cox. ", "Among remaining methods with coverage in [", scales::percent(best_cov_band[1], accuracy = 0.1), ", ", scales::percent(best_cov_band[2], accuracy = 0.1), "], the method with the smallest RMSE is selected; if no method falls in this band, ", 
        "the method closest to 95% coverage is selected, with RMSE as the tie-breaker. ", "Coverage greater than ", scales::percent(overcov_red_threshold, accuracy = 0.1), " is shown in red to flag conservative overcoverage.")
    kb <- kbl(tab_fmt, booktabs = TRUE, escape = FALSE, caption = caption_txt) %>% kable_classic_2(full_width = FALSE, html_font = "Arial") %>% footnote(general = note_txt, general_title = "Notes: ", footnote_as_chunk = TRUE, escape = FALSE)
    for (g in seq_len(nrow(group_info))) {
        kb <- kb %>% group_rows(group_label = group_info$group_label[g], start_row = group_info$start_row[g], end_row = group_info$end_row[g], bold = TRUE, italic = FALSE)
    }
    best_rows <- which(tab$Is_Best)
    for (rr in best_rows) {
        kb <- kb %>% row_spec(row = rr, bold = TRUE, background = "#E8F5E9")
    }
    kb
}

assign_aesthetic <- function(label, color, shape, linetype) {
    if (label %in% names(method_colors)) {
        method_colors[label] <<- color
        method_shapes[label] <<- shape
        method_linetypes[label] <<- linetype
    }
}

abbr_tag <- function(x) {
    x <- as.character(x)
    dict <- c(linear = "lin", nonlinear = "nl", linear_ph = "lph", nonlinear_ph = "nlph", glm = "glm", glm_lasso = "glmlas", nnet = "nn", rf = "rf", bart = "bart", cox_lp = "coxlp", lasso_cox_lp = "lcoxlp", quadratic = "quad", kl = "kl", el = "el", hellinger = "hel", renyi = "renyi05", lambda.1se = "1se", lambda.min = "min")
    out <- unname(dict[x])
    ifelse(is.na(out), x, out)
}

num_tag <- function(x) {
    gsub("\\.", "p", format(x, trim = TRUE, scientific = FALSE))
}

gen_tag <- function(x) {
    paste(abbr_tag(x), collapse = "-")
}

