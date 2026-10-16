expit <- function(x) {
    1/(1 + exp(-x))
}

clip01 <- function(x, lower = 0.01, upper = 0.99) {
    pmin(pmax(as.numeric(x), lower), upper)
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

