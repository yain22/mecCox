# These scripts are installed as reproduction examples, outside the package API.
helper_path <- testthat::test_path("..", "..", "inst", "reproduce",
                                  "simulation_helpers.R")
if (!file.exists(helper_path)) {
  helper_path <- system.file("reproduce", "simulation_helpers.R",
                             package = "mecCox")
}
helper_path <- normalizePath(helper_path, mustWork = TRUE)
simulation_helpers <- new.env(parent = baseenv())
source(helper_path, local = simulation_helpers)

testthat::test_that("simulation arguments have useful defaults and reject mistakes", {
  parse_arguments <- simulation_helpers$parse_simulation_arguments
  defaults <- parse_arguments("example-output", character())
  testthat::expect_identical(
    defaults,
    list(quick_run = FALSE, output_directory = "example-output", cores = 20L)
  )
  testthat::expect_identical(
    parse_arguments(NULL, character()),
    list(quick_run = FALSE, output_directory = NULL, cores = 20L)
  )

  chosen <- parse_arguments(
    "example-output", c("--quick", "--cores=2", "--output=a folder")
  )
  testthat::expect_identical(
    chosen,
    list(quick_run = TRUE, output_directory = "a folder", cores = 2L)
  )

  invalid_arguments <- list(
    "--cores=0", "--cores=-1", "--cores=1.5", "--cores=NaN",
    "--cores=Inf", "--cores=", "--cores=2147483648", "--output=",
    "--unknown", c("--cores=1", "--cores=2"),
    c("--output=a", "--output=b"), c("--quick", "--quick")
  )
  for (arguments in invalid_arguments) {
    testthat::expect_error(parse_arguments("example-output", arguments),
                           "cores|Use")
  }
})

testthat::test_that("both scripts can be sourced from another working directory", {
  reproduction_directory <- dirname(helper_path)
  scratch <- tempfile("simulation source ")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  previous_directory <- setwd(tempdir())
  on.exit(setwd(previous_directory), add = TRUE)

  for (scenario in c("scenario1", "scenario2")) {
    script <- file.path(reproduction_directory, paste0(scenario, ".R"))
    expressions <- as.list(parse(script))
    settings_index <- which(vapply(expressions, function(expression) {
      is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        identical(expression[[2L]], as.name("options"))
    }, logical(1)))
    # Exercise the real startup code without starting a full simulation in tests.
    startup <- expressions[seq_len(settings_index)]
    script_copy <- file.path(scratch, paste0(scenario, ".R"))
    writeLines(unlist(lapply(startup, deparse)), script_copy)
    file.copy(helper_path, file.path(scratch, "simulation_helpers.R"),
              overwrite = TRUE)

    execution <- new.env(parent = baseenv())
    execution$commandArgs <- function(trailingOnly = FALSE) {
      if (trailingOnly) "--unrelated-wrapper-option" else c("R", "--file=wrapper.R")
    }
    source(script_copy, local = execution)
    testthat::expect_identical(
      execution$options,
      list(quick_run = FALSE, output_directory = NULL,
           cores = 20L)
    )
    testthat::expect_identical(
      normalizePath(execution$helper_candidates[1L]),
      normalizePath(file.path(scratch, "simulation_helpers.R"))
    )

    # Settings edited in the file also work with RStudio's Source/chdir mode.
    startup[[1L]] <- quote(quick_run <- TRUE)
    startup[[2L]] <- quote(cores <- 2L)
    startup[[3L]] <- quote(output_directory <- "quick output")
    writeLines(unlist(lapply(startup, deparse)), script_copy)
    source(script_copy, local = execution, chdir = TRUE)
    testthat::expect_identical(
      execution$options,
      list(quick_run = TRUE, output_directory = "quick output", cores = 2L)
    )
  }
})

testthat::test_that("worker limits respect the request and number of replications", {
  available <- parallel::detectCores(logical = TRUE)
  expected <- 20L
  if (length(available) == 1L && is.finite(available) && available > 0L) {
    expected <- as.integer(min(expected, available))
  }
  testthat::expect_identical(
    simulation_helpers$choose_worker_count(20L, 1000L), expected
  )
  testthat::expect_identical(
    simulation_helpers$choose_worker_count(20L, 1L), 1L
  )
  testthat::expect_identical(
    simulation_helpers$choose_worker_count(1L, 1000L), 1L
  )
  testthat::expect_lte(simulation_helpers$choose_worker_count(3L, 2L), 2L)
  testthat::expect_null(
    simulation_helpers$start_simulation_cluster(1L, character())
  )
})

testthat::test_that("simulation plots preserve the active graphics device", {
  grDevices::pdf(NULL, width = 14, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  device <- grDevices::dev.cur()
  old <- graphics::par(c("mar", "oma", "mfrow", "las"))
  first_summary <- data.frame(
    ratio = "1:2", n1 = 200L,
    method = c("Naive", "Robust sandwich", "Corrected sandwich", "MEC-Cox"),
    coverage = c(0.8, 0.9, 0.95, 0.95),
    bias = c(0.1, 0.1, 0.1, 0.05), rmse = c(0.2, 0.2, 0.2, 0.1)
  )
  simulation_helpers$plot_scenario1_results(first_summary)
  testthat::expect_identical(grDevices::dev.cur(), device)
  testthat::expect_equal(graphics::par(names(old)), old)

  settings <- data.frame(setting = c("None", "Mild", "Severe"))
  second_summary <- expand.grid(
    setting = settings$setting,
    method = c("Naive", "Robust sandwich", "Corrected sandwich",
               "MEC-Cox (BART/Cox)", "MEC-Cox (BART/RSF)"),
    stringsAsFactors = FALSE
  )
  second_summary$n1 <- 200L
  second_summary$coverage <- 0.95
  second_summary$bias <- 0
  second_summary$rmse <- 0.1
  simulation_helpers$plot_scenario2_results(second_summary, settings)
  testthat::expect_identical(grDevices::dev.cur(), device)
  testthat::expect_equal(graphics::par(names(old)), old)

  # RStudio panes can be much smaller than the optional publication PDF.
  grDevices::dev.off()
  grDevices::pdf(NULL, width = 4, height = 3)
  device <- grDevices::dev.cur()
  old <- graphics::par(c("mar", "oma", "mfrow", "las", "cex", "mgp"))
  full_first_summary <- do.call(rbind, lapply(c("1:2", "1:3", "1:4"), function(ratio) {
    rows <- first_summary
    rows$ratio <- ratio
    rows
  }))
  simulation_helpers$plot_scenario1_results(full_first_summary)
  testthat::expect_identical(grDevices::dev.cur(), device)
  testthat::expect_equal(graphics::par(names(old)), old)
  simulation_helpers$plot_scenario2_results(second_summary, settings)
  testthat::expect_identical(grDevices::dev.cur(), device)
  testthat::expect_equal(graphics::par(names(old)), old)
})

testthat::test_that("simulation drivers retain results without creating files", {
  scratch <- tempfile("simulation display ")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  previous_directory <- setwd(scratch)
  on.exit(setwd(previous_directory), add = TRUE)

  for (scenario in c("scenario1", "scenario2")) {
    expressions <- as.list(parse(file.path(dirname(helper_path),
                                           paste0(scenario, ".R"))))
    execution <- new.env(parent = simulation_helpers)
    for (expression in expressions) {
      if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
          is.call(expression[[3L]]) &&
          identical(expression[[3L]][[1L]], as.name("function"))) {
        eval(expression, envir = execution)
      }
    }
    execution$compute_reference_target <- function(...) 0
    execution$options <- list(quick_run = TRUE, cores = 1L,
                               output_directory = NULL)
    execution$design <- list(replications = 1L, treated_sizes = 2L,
                              control_multipliers = 2L, control_multiplier = 4L)
    execution$settings <- data.frame(setting_id = 1L, setting = "None",
                                     kappa_pi = 0, kappa_m = 0)
    if (scenario == "scenario1") {
      execution$run_replication <- function(multiplier, treated_count, replicate,
                                             target, design) {
        make_result(multiplier, treated_count, replicate, "MEC-Cox", target,
                     estimate = 0.1, standard_error = 0.2)
      }
    } else {
      execution$run_replication <- function(setting, treated_count, replicate,
                                             target, design) {
        make_result(setting, treated_count, replicate, "MEC-Cox (BART/Cox)",
                     target, design, estimate = 0.1, standard_error = 0.2)
      }
    }
    environment(execution$run_replication) <- execution

    # Run the script's final assignments, with cheap deterministic fits, so
    # the test checks the real driver and the objects exposed to the user.
    result_names <- paste0(scenario, c("_results", "_summary", "_replications",
                                       "_targets"))
    for (expression in expressions) {
      if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
          as.character(expression[[2L]]) %in% result_names) {
        invisible(capture.output(suppressMessages(eval(expression, execution))))
      }
    }
    answer <- execution[[paste0(scenario, "_results")]]
    testthat::expect_s3_class(answer$summary, "data.frame")
    testthat::expect_s3_class(answer$replications, "data.frame")
    testthat::expect_type(answer$metadata, "list")
    testthat::expect_identical(execution[[paste0(scenario, "_summary")]],
                               answer$summary)
    testthat::expect_identical(execution[[paste0(scenario, "_replications")]],
                               answer$replications)
    testthat::expect_length(list.files(scratch, all.files = TRUE, no.. = TRUE), 0L)
  }
})

testthat::test_that("serial and socket workers preserve per-replication RNG results", {
  testthat::skip_on_cran()

  check_execution <- function() {
    previous_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      do.call(RNGkind, as.list(previous_kind))
      if (had_seed) {
        assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")

    thread_variables <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                          "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                          "BLIS_NUM_THREADS", "RCPP_PARALLEL_NUM_THREADS")
    original_threads <- Sys.getenv(thread_variables, unset = NA_character_)
    on.exit({
      unset <- is.na(original_threads)
      Sys.unsetenv(names(original_threads)[unset])
      if (any(!unset)) {
        do.call(Sys.setenv, as.list(original_threads[!unset]))
      }
    }, add = TRUE)
    Sys.setenv(OMP_NUM_THREADS = "3")
    Sys.unsetenv("OPENBLAS_NUM_THREADS")
    previous_threads <- Sys.getenv(thread_variables, unset = NA_character_)
    previous_connections <- rownames(showConnections(all = TRUE))

    worker_functions <- new.env(parent = baseenv())
    worker_functions$replication_seed <- function(seed, replicate) {
      seed + replicate
    }
    worker_functions$draw_replication <- function(seed, replicate) {
      set.seed(replication_seed(seed, replicate))
      list(replicate = replicate, kind = RNGkind(),
           normal = stats::rnorm(5L), binomial = stats::rbinom(5L, 1L, 0.4),
           sampled = sample.int(100L, 5L))
    }
    environment(worker_functions$draw_replication) <- worker_functions

    cluster <- simulation_helpers$start_simulation_cluster(
      2L, c("replication_seed", "draw_replication"), envir = worker_functions
    )
    on.exit({
      if (!is.null(cluster)) parallel::stopCluster(cluster)
    }, add = TRUE)

    testthat::expect_identical(
      Sys.getenv(thread_variables, unset = NA_character_), previous_threads
    )
    worker_state <- parallel::clusterCall(cluster, function(thread_variables) {
      list(kind = RNGkind(), threads = Sys.getenv(thread_variables),
           package_loaded = isNamespaceLoaded("mecCox"),
           helper_exported = exists("replication_seed", envir = .GlobalEnv,
                                    inherits = FALSE))
    }, thread_variables)
    for (state in worker_state) {
      testthat::expect_identical(state$kind, RNGkind())
      testthat::expect_true(all(state$threads == "1"))
      testthat::expect_true(state$package_loaded)
      testthat::expect_true(state$helper_exported)
    }

    serial <- simulation_helpers$run_simulation_replications(
      7L, worker_functions$draw_replication, arguments = list(seed = 8451L)
    )
    concurrent <- simulation_helpers$run_simulation_replications(
      7L, worker_functions$draw_replication, arguments = list(seed = 8451L),
      cluster = cluster
    )
    testthat::expect_identical(concurrent, serial)
    testthat::expect_identical(
      vapply(concurrent, function(result) result$replicate, integer(1)), 1:7
    )
    testthat::expect_false(identical(serial[[1L]]$normal, serial[[2L]]$normal))

    parallel::stopCluster(cluster)
    cluster <- NULL
    testthat::expect_setequal(rownames(showConnections(all = TRUE)),
                              previous_connections)
  }
  check_execution()
})

testthat::test_that("failed worker initialization restores environment and sockets", {
  testthat::skip_on_cran()

  thread_variables <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                        "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                        "BLIS_NUM_THREADS", "RCPP_PARALLEL_NUM_THREADS")
  previous_threads <- Sys.getenv(thread_variables, unset = NA_character_)
  previous_connections <- rownames(showConnections(all = TRUE))
  testthat::expect_error(
    simulation_helpers$start_simulation_cluster(
      2L, "missing_simulation_helper", envir = new.env(parent = emptyenv())
    ),
    "missing_simulation_helper"
  )
  testthat::expect_identical(
    Sys.getenv(thread_variables, unset = NA_character_), previous_threads
  )
  testthat::expect_setequal(rownames(showConnections(all = TRUE)),
                            previous_connections)
})
