# Run the public breast-cancer case study and display the two supplementary tables.
# In RStudio, edit these settings and use Source. From a terminal:
#   Rscript path/to/mecCox/inst/reproduce/breast_cancer.R --seed=20260427
# Data are loaded from survival; no patient-data download or result export is used.
# Numerical results can vary with seeds, learner settings, and software versions.
seed <- 20260427L
n_folds <- 10L
n_landmarks <- 20L
mlp_epochs <- 200L
rsf_auto_tune <- TRUE

source_files <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
}, character(1))
source_files <- rev(source_files[nzchar(source_files)])
script_option <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_files <- c(source_files, sub("^--file=", "", script_option))
helper_directories <- unique(c(
  dirname(script_files), ".", file.path("inst", "reproduce"),
  file.path("mecCox", "inst", "reproduce"),
  system.file("reproduce", package = "mecCox")
))
for (helper in c("simulation_helpers.R", "breast_cancer_helpers.R")) {
  candidates <- file.path(helper_directories[nzchar(helper_directories)], helper)
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    stop("Cannot find ", helper, ". Keep both helpers beside breast_cancer.R ",
         "or install the current mecCox package.", call. = FALSE)
  }
  source(candidates[1L], local = TRUE)
}

arguments <- if (interactive() || length(source_files)) {
  character()
} else commandArgs(trailingOnly = TRUE)
if (length(arguments)) {
  if (length(arguments) != 1L || !grepl("^--seed=[0-9]+$", arguments)) {
    stop("The supported command-line option is --seed=<nonnegative integer>.",
         call. = FALSE)
  }
  seed <- as.numeric(sub("^--seed=", "", arguments))
}

prepare_breast_cancer_data <- function() {
  public_data <- new.env(parent = baseenv())
  utils::data("cancer", package = "survival", envir = public_data)
  gb <- public_data$gbsg[public_data$gbsg$hormon == 1L, ]
  rot <- public_data$rotterdam[public_data$rotterdam$hormon == 0L, ]

  # Match the original case study's RFS convention: first recorded recurrence
  # or death, otherwise the later recorded follow-up time. The survival help
  # page discusses the limitations of deaths recorded after recurrence follow-up.
  first_event <- pmin(ifelse(rot$recur == 1L, rot$rtime, Inf),
                      ifelse(rot$death == 1L, rot$dtime, Inf))
  rot_event <- as.integer(is.finite(first_event))
  rot_time <- ifelse(rot_event == 1L, first_event,
                     pmax(rot$rtime, rot$dtime, na.rm = TRUE))
  time_raw <- c(gb$rfstime, rot_time) / 365.25
  event_raw <- c(gb$status, rot_event)
  size_labels <- c("<=20 mm", "20-50 mm", ">50 mm")
  gb_size <- ifelse(gb$size <= 20, size_labels[1L],
                     ifelse(gb$size <= 50, size_labels[2L], size_labels[3L]))
  rot_size <- factor(rot$size, levels = c("<=20", "20-50", ">50"),
                      labels = size_labels)
  data <- data.frame(
    time = pmin(time_raw, 5),
    delta = as.integer(event_raw == 1L & time_raw <= 5),
    A = c(rep(1L, nrow(gb)), rep(0L, nrow(rot))),
    age = as.numeric(c(gb$age, rot$age)),
    meno = factor(c(gb$meno, rot$meno), levels = c(0, 1),
                    labels = c("Premenopausal", "Postmenopausal")),
    size_cat = factor(c(gb_size, as.character(rot_size)), levels = size_labels),
    grade = factor(c(gb$grade, rot$grade), levels = c(1, 2, 3)),
    log_nodes = log1p(c(gb$nodes, rot$nodes)),
    log_pgr = log1p(c(gb$pgr, rot$pgr)),
    log_er = log1p(c(gb$er, rot$er))
  )
  if (anyNA(data) || any(data$time <= 0)) {
    stop("The public case-study data contain missing or nonpositive inputs.",
         call. = FALSE)
  }
  if (sum(data$A == 1L) != 246L || sum(data$A == 0L) != 2643L) {
    stop("Public dataset cohort sizes differ from the case-study specification.",
         call. = FALSE)
  }
  data
}

breast_cancer_estimate_table <- function(fits) {
  row <- function(method, theta, se) {
    interval <- exp(theta + c(-1, 1) * stats::qnorm(0.975) * se)
    data.frame(Method = method, log_HR = unname(theta), SE_log_HR = unname(se),
               HR = exp(unname(theta)), CI_lower = interval[1L],
               CI_upper = interval[2L], row.names = NULL)
  }
  do.call(rbind, list(
    row("Unweighted Cox", fits$unweighted$theta, fits$unweighted$se),
    row("Naive", fits$ipw$theta, fits$ipw$se["naive"]),
    row("Robust sandwich", fits$ipw$theta, fits$ipw$se["robust"]),
    row("Corrected sandwich", fits$ipw$theta, fits$ipw$se["corrected"]),
    row("MEC-Cox (GLM/Cox)", fits$glm_cox$theta, fits$glm_cox$se),
    row("MEC-Cox (DL/RSF)", fits$dl_rsf$theta, fits$dl_rsf$se)
  ))
}

breast_cancer_balance_table <- function(data, fits) {
  # Include every reported category, including reference categories, as in S3.
  design <- cbind(
    age = data$age, meno = as.integer(data$meno == "Postmenopausal"),
    size1 = as.integer(data$size_cat == "<=20 mm"),
    size2 = as.integer(data$size_cat == "20-50 mm"),
    size3 = as.integer(data$size_cat == ">50 mm"),
    grade1 = as.integer(data$grade == "1"),
    grade2 = as.integer(data$grade == "2"),
    grade3 = as.integer(data$grade == "3"),
    log_nodes = data$log_nodes, log_pgr = data$log_pgr, log_er = data$log_er
  )
  target <- data$A == 1L
  external <- !target
  continuous <- c(1L, 9L, 10L, 11L)
  target_mean <- colMeans(design[target, , drop = FALSE])
  target_sd <- apply(design[target, , drop = FALSE], 2L, stats::sd)
  weights <- list(Unweighted = rep(1, nrow(data)), `ATT-IPW` = fits$ipw$weights,
                   `MEC-Cox (GLM/Cox)` = fits$glm_cox$weights,
                   `MEC-Cox (DL/RSF)` = fits$dl_rsf$weights)
  values <- lapply(seq_along(weights), function(index) {
    w <- weights[[index]][external]
    if (length(w) != sum(external) || any(!is.finite(w)) || any(w < 0) ||
        sum(w) <= 0) stop("Invalid external-control weights.", call. = FALSE)
    difference <- target_mean - colSums(design[external, , drop = FALSE] * w) / sum(w)
    difference[continuous] <- difference[continuous] / target_sd[continuous]
    c(difference, mean(abs(difference)), max(abs(difference)),
      sum(abs(difference) > 0.10),
      if (index == 1L) NA_real_ else sum(w)^2 / sum(w^2),
      if (index == 1L) NA_real_ else stats::sd(w) / mean(w))
  })
  names(values) <- names(weights)
  data.frame(
    Covariate = c("age", "meno: postmenopausal", "size: <=20 mm", "size: 20-50 mm",
                  "size: >50 mm", "grade 1", "grade 2", "grade 3", "log(1 + nodes)",
                  "log(1 + PGR)", "log(1 + ER)", "Mean |SMD|", "Max |SMD|",
                  "No. |SMD| > 0.10", "External-control ESS", "External-control weight CV"),
    Type = c("C", rep("B", 7L), rep("C", 3L), rep("-", 5L)),
    values, check.names = FALSE, row.names = NULL
  )
}

build_breast_cancer_report <- function(estimates, balance, settings) {
  fixed <- function(x) {
    # Avoid displaying numerical negative zero; keep full precision in R objects.
    x[!is.na(x) & abs(x) < 0.0005] <- 0
    ifelse(is.na(x), "-", sprintf("%.3f", x))
  }
  estimates_display <- data.frame(
    Method = estimates$Method, `Log HR` = fixed(estimates$log_HR),
    `SE (log HR)` = fixed(estimates$SE_log_HR), HR = fixed(estimates$HR),
    `95% CI of HR` = paste0("(", fixed(estimates$CI_lower), ", ",
                            fixed(estimates$CI_upper), ")"), check.names = FALSE
  )
  balance_display <- balance
  for (column in 3:ncol(balance_display)) {
    balance_display[[column]] <- fixed(balance[[column]])
    balance_display[14L, column] <- as.character(as.integer(balance[14L, column]))
  }
  style <- function(data, title) {
    table <- kableExtra::kbl(data, format = "html", row.names = FALSE,
                             caption = title, escape = TRUE)
    kableExtra::kable_styling(table, full_width = FALSE, position = "left",
                              bootstrap_options = c("striped", "hover", "condensed"))
  }
  tables <- list(
    estimates = style(estimates_display, "Table S2. Marginal hazard-ratio estimates"),
    balance = style(balance_display, "Table S3. Covariate balance and weight diagnostics")
  )
  report <- htmltools::browsable(htmltools::tagList(
    htmltools::tags$head(htmltools::tags$title("Public breast-cancer case study"),
      htmltools::tags$style(htmltools::HTML(paste(
        "body {font-family: Arial, sans-serif; color: #222; padding: 20px;}",
        "section {margin: 24px 0; overflow-x: auto;}",
        "table {border-collapse: collapse;} th, td {padding: 7px 10px;",
        "border-bottom: 1px solid #ddd; white-space: nowrap;}",
        "caption {text-align: left; font-weight: bold; padding: 10px 0;}",
        "tbody tr:nth-child(even) {background: #f5f7fa;}"
      )))),
    htmltools::tags$h1("Public breast-cancer case study"),
    htmltools::tags$p("GBSG hormonal therapy (n = 246) and Rotterdam without hormonal therapy ",
                      "(n = 2643); recurrence-free survival through 5 years. Data: survival package."),
    htmltools::tags$p(sprintf("Seed: %d; cross-fitting folds: %d; survival landmarks: %d.",
                              settings$seed, settings$n_folds, settings$n_landmarks)),
    htmltools::tags$section(htmltools::HTML(as.character(tables$estimates))),
    htmltools::tags$p("Unweighted Cox is descriptive. The three ATT-IPW rows share one ",
      "coefficient and differ in variance estimation. MEC-Cox uses KL calibration."),
    htmltools::tags$section(htmltools::HTML(as.character(tables$balance))),
    htmltools::tags$p("C: continuous (difference divided by target-cohort SD). ",
      "B: binary (difference in proportions). Differences are target minus weighted ",
      "external control. ESS and weight CV concern external controls only."),
    htmltools::tags$p("Grade 1 occurs in the treated cohort but not among controls; its ",
      "imbalance cannot be removed by weighting. The case-study helper retains the ",
      "original outcome-model extrapolation convention. See the accompanying README."),
    htmltools::tags$p("These tables are generated from the fitted models. Seeds, learner ",
      "settings, and software versions can change the numerical results.")
  ))
  list(tables = tables, report = report)
}

run_breast_cancer <- function(settings) {
  packages <- c("mecCox", "survival", "MASS", "brulee", "torch", "ranger",
                "kableExtra", "htmltools", "rstudioapi")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the following packages first: ",
                            paste(missing, collapse = ", "), call. = FALSE)
  if (!torch::torch_is_installed()) {
    stop("Install the Torch runtime with torch::install_torch() before running this study.",
         call. = FALSE)
  }
  if (length(settings$seed) != 1L || !is.finite(settings$seed) || settings$seed < 0 ||
      settings$seed > .Machine$integer.max - settings$n_folds ||
      settings$seed != as.integer(settings$seed)) stop("seed must be a nonnegative integer.")
  data <- prepare_breast_cancer_data()
  covariates <- c("age", "meno", "size_cat", "grade", "log_nodes", "log_pgr", "log_er")
  message("Fitting the unweighted and ATT-IPW Cox models ...")
  fits <- list(
    unweighted = mecCox::fit_unweighted_cox(data, "time", "delta", "A"),
    ipw = mecCox::fit_att_ipw_cox(data, "time", "delta", "A", covariates)
  )
  message("Fitting MEC-Cox with GLM/Cox ...")
  fits$glm_cox <- do.call(fit_breast_mec, c(list(data = data, covariates = covariates,
    ps_learner = "glm", survival_learner = "cox"), settings))
  message("Fitting MEC-Cox with DL/RSF ...")
  fits$dl_rsf <- do.call(fit_breast_mec, c(list(data = data, covariates = covariates,
    ps_learner = "mlp", survival_learner = "rsf"), settings))
  estimates <- breast_cancer_estimate_table(fits)
  balance <- breast_cancer_balance_table(data, fits)
  print(estimates, row.names = FALSE, digits = 4)
  print(balance, row.names = FALSE, digits = 4)
  list(data = data, fits = fits, estimates = estimates, balance = balance,
       metadata = list(settings = settings, session = utils::sessionInfo(),
                       cohort = c(GBSG = 246L, Rotterdam = 2643L)))
}

breast_cancer_settings <- list(seed = seed, n_folds = n_folds,
  n_landmarks = n_landmarks, mlp_epochs = mlp_epochs, rsf_auto_tune = rsf_auto_tune)
breast_cancer_results <- run_breast_cancer(breast_cancer_settings)
breast_cancer_estimates <- breast_cancer_results$estimates
breast_cancer_balance <- breast_cancer_results$balance
breast_cancer_display <- build_breast_cancer_report(
  breast_cancer_estimates, breast_cancer_balance, breast_cancer_settings)
breast_cancer_tables <- breast_cancer_display$tables
breast_cancer_report <- breast_cancer_display$report
show_simulation_report(breast_cancer_report)
