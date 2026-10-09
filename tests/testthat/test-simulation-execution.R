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
  defaults <- parse_arguments()
  testthat::expect_identical(
    defaults,
    list(quick_run = FALSE, cores = 20L, replications = 1000L)
  )
  wrapper <- new.env(parent = simulation_helpers)
  wrapper$commandArgs <- function(...) "--unrelated-wrapper-option"
  isolated_parser <- parse_arguments
  environment(isolated_parser) <- wrapper
  testthat::expect_identical(isolated_parser(), defaults)

  chosen <- parse_arguments(
    c("--quick", "--cores=2", "--replications=7")
  )
  testthat::expect_identical(
    chosen,
    list(quick_run = TRUE, cores = 2L, replications = 7L)
  )
  testthat::expect_identical(
    parse_arguments(quick_run = TRUE, cores = 3L, replications = 11L),
    list(quick_run = TRUE, cores = 3L, replications = 11L)
  )
  testthat::expect_identical(
    parse_arguments(c("--cores=2", "--replications=5"),
                    cores = 3L, replications = 11L),
    list(quick_run = FALSE, cores = 2L, replications = 5L)
  )

  invalid_arguments <- list(
    "--cores=0", "--cores=-1", "--cores=1.5", "--cores=NaN",
    "--cores=Inf", "--cores=", "--cores=2147483648",
    "--replications=0", "--replications=-1", "--replications=1.5",
    "--replications=NaN", "--replications=Inf", "--replications=",
    "--replications=2147483648", "--output=", "--output=example-output",
    "--unknown", c("--cores=1", "--cores=2"),
    c("--replications=2", "--replications=3"), c("--quick", "--quick")
  )
  for (arguments in invalid_arguments) {
    testthat::expect_error(parse_arguments(arguments), "cores|replications|Use")
  }
})

testthat::test_that("quick checks retain explicitly selected replication counts", {
  parse_arguments <- simulation_helpers$parse_simulation_arguments
  testthat::expect_identical(parse_arguments("--quick")$replications, 2L)
  testthat::expect_identical(
    parse_arguments(quick_run = TRUE)$replications, 2L
  )
  testthat::expect_identical(
    parse_arguments("--quick", replications = 9L)$replications, 9L
  )
  testthat::expect_identical(
    parse_arguments(c("--quick", "--replications=1000"))$replications, 1000L
  )
  testthat::expect_identical(
    parse_arguments("--replications=4", quick_run = TRUE,
                    replications = 9L)$replications, 4L
  )
  testthat::expect_error(parse_arguments(cores = 0L), "cores")
  testthat::expect_error(parse_arguments(replications = 0L), "replications")
  testthat::expect_error(parse_arguments(replications = 1.5), "replications")
  testthat::expect_error(parse_arguments(replications = NA_integer_),
                         "replications")
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
      list(quick_run = FALSE, cores = 20L, replications = 1000L)
    )
    testthat::expect_identical(
      normalizePath(execution$helper_candidates[1L]),
      normalizePath(file.path(scratch, "simulation_helpers.R"))
    )

    # Settings edited in the file also work with RStudio's Source/chdir mode.
    startup[[1L]] <- quote(quick_run <- TRUE)
    startup[[2L]] <- quote(cores <- 2L)
    startup[[3L]] <- quote(replications <- 7L)
    writeLines(unlist(lapply(startup, deparse)), script_copy)
    source(script_copy, local = execution, chdir = TRUE)
    testthat::expect_identical(
      execution$options,
      list(quick_run = TRUE, cores = 2L, replications = 7L)
    )

    # Reduced quick-run designs still use the count edited by the user.
    design_index <- which(vapply(expressions, function(expression) {
      is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        identical(expression[[2L]], as.name("design"))
    }, logical(1)))
    first_function <- which(vapply(expressions, function(expression) {
      is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        is.call(expression[[3L]]) &&
        identical(expression[[3L]][[1L]], as.name("function"))
    }, logical(1)))[1L]
    design_expressions <- expressions[seq.int(design_index, first_function - 1L)]
    for (expression in design_expressions) {
      suppressMessages(eval(expression, envir = execution))
    }
    testthat::expect_identical(execution$design$replications, 7L)
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

  # RStudio panes can be much smaller than an explicitly opened device.
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

testthat::test_that("headless plotting does not open an automatic PDF device", {
  scratch <- tempfile("simulation headless ")
  dir.create(scratch)
  script <- tempfile("simulation headless check ", fileext = ".R")
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "arguments <- commandArgs(trailingOnly = TRUE)",
    "setwd(arguments[2L])",
    "helpers <- new.env(parent = baseenv())",
    "source(arguments[1L], local = helpers)",
    "stopifnot(grDevices::dev.cur() == 1L)",
    "first <- data.frame(ratio = '1:2', n1 = 200L, method = 'MEC-Cox',",
    "                    coverage = 0.95, bias = 0, rmse = 0.1)",
    "second <- data.frame(setting = 'None', n1 = 200L,",
    "                     method = 'MEC-Cox (BART/Cox)',",
    "                     coverage = 0.95, bias = 0, rmse = 0.1)",
    "first_drawn <- helpers$plot_scenario1_results(first)",
    "second_drawn <- helpers$plot_scenario2_results(",
    "  second, data.frame(setting = 'None'))",
    "stopifnot(identical(first_drawn, FALSE), identical(second_drawn, FALSE),",
    "          grDevices::dev.cur() == 1L,",
    "          length(list.files('.', all.files = TRUE, no.. = TRUE)) == 0L)",
    "cat('Headless plots left no files.\\n')"
  ), script)

  executable <- file.path(R.home("bin"),
                          if (.Platform$OS.type == "windows") {
                            "Rscript.exe"
                          } else "Rscript")
  output <- system2(executable,
                    c("--vanilla", shQuote(script), shQuote(helper_path),
                      shQuote(scratch)), stdout = TRUE, stderr = TRUE)
  testthat::expect_null(attr(output, "status"))
  testthat::expect_true(any(grepl("Headless plots left no files", output,
                                 fixed = TRUE)))
  testthat::expect_length(list.files(scratch, all.files = TRUE, no.. = TRUE), 0L)
})

testthat::test_that("simulation summaries are formatted HTML tables", {
  testthat::skip_if_not_installed("kableExtra")
  summary <- data.frame(
    ratio = "1:2", n1 = 200L, n0 = 400L, method = "MEC-Cox",
    replications = 1000L, successful = 1000L, failed = 0L,
    coverage = 0.956789, bias = 0.0123456, rmse = 0.1234567
  )
  table <- simulation_helpers$format_simulation_table(
    summary, "Scenario 1: simulation summary"
  )
  testthat::expect_s3_class(table, "kableExtra")
  html <- paste(as.character(table), collapse = "\n")
  testthat::expect_match(html, "<table")
  testthat::expect_match(html, "Scenario 1: simulation summary", fixed = TRUE)
  testthat::expect_match(html, "MEC-Cox", fixed = TRUE)
  testthat::expect_match(html, "coverage", ignore.case = TRUE)
  testthat::expect_false(grepl("0.956789", html, fixed = TRUE))
})

testthat::test_that("simulation reports invoke the Viewer using temporary HTML", {
  testthat::skip_if_not_installed("kableExtra")
  testthat::skip_if_not_installed("htmltools")
  scratch <- tempfile("simulation viewer ")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  previous_directory <- setwd(scratch)
  on.exit(setwd(previous_directory), add = TRUE)

  table <- simulation_helpers$format_simulation_table(
    data.frame(Method = "MEC-Cox", Coverage = 0.95), "Simulation summary"
  )
  report <- htmltools::tagList(htmltools::h1("Scenario 1"),
                               htmltools::HTML(as.character(table)))
  viewed_path <- NULL
  viewed_html <- NULL
  viewer <- function(path) {
    viewed_path <<- normalizePath(path, winslash = "/", mustWork = TRUE)
    viewed_html <<- paste(readLines(path, warn = FALSE), collapse = "\n")
  }
  shown <- simulation_helpers$show_simulation_report(report, viewer = viewer)
  testthat::expect_true(shown)
  testthat::expect_type(viewed_path, "character")
  testthat::expect_true(startsWith(viewed_path,
                                  normalizePath(tempdir(), winslash = "/")))
  testthat::expect_match(viewed_html, "Scenario 1", fixed = TRUE)
  testthat::expect_match(viewed_html, "<table")
  testthat::expect_length(list.files(scratch, all.files = TRUE, no.. = TRUE), 0L)
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
    execution$options <- list(quick_run = TRUE, cores = 1L, replications = 1L)
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
                                       "_targets", "_display", "_tables",
                                       "_report"))
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
    testthat::expect_identical(answer$metadata$design$replications, 1L)
    tables <- execution[[paste0(scenario, "_tables")]]
    expected_tables <- c("summary", "targets", "configuration")
    if (scenario == "scenario2") expected_tables <- c(expected_tables, "settings")
    testthat::expect_setequal(names(tables), expected_tables)
    for (table in tables) testthat::expect_s3_class(table, "kableExtra")
    report <- execution[[paste0(scenario, "_report")]]
    report_html <- htmltools::renderTags(report)$html
    testthat::expect_match(report_html, "simulation summary", fixed = TRUE)
    testthat::expect_match(report_html, "reference target", fixed = TRUE)
    testthat::expect_match(report_html, "run configuration", fixed = TRUE)
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

testthat::test_that("completion callbacks preserve serial and socket RNG results", {
  testthat::skip_on_cran()
  draw <- function(replicate, seed) {
    # Different runtimes exercise receiving results out of task order.
    if (replicate == 1L) Sys.sleep(0.15)
    set.seed(seed + replicate)
    list(replicate = replicate, draws = stats::rnorm(3L))
  }
  environment(draw) <- baseenv()
  serial_received <- list()
  serial <- simulation_helpers$run_simulation_replications(
    5L, draw, list(seed = 341L),
    on_result = function(replicate, value) {
      serial_received[[as.character(replicate)]] <<- value
    }
  )
  cluster <- simulation_helpers$start_simulation_cluster(2L, character())
  on.exit(simulation_helpers$stop_simulation_cluster(cluster), add = TRUE)
  socket_received <- list()
  concurrent <- simulation_helpers$run_simulation_replications(
    5L, draw, list(seed = 341L), cluster = cluster,
    on_result = function(replicate, value) {
      socket_received[[as.character(replicate)]] <<- value
    }
  )
  testthat::expect_identical(concurrent, serial)
  testthat::expect_identical(unname(serial_received), serial)
  testthat::expect_setequal(names(socket_received), as.character(1:5))
  testthat::expect_identical(unname(socket_received[as.character(1:5)]), serial)
})

testthat::test_that("socket interruption preserves already received completions", {
  testthat::skip_on_cran()
  previous_connections <- rownames(showConnections(all = TRUE))
  cluster <- simulation_helpers$start_simulation_cluster(2L, character())
  on.exit(simulation_helpers$stop_simulation_cluster(cluster), add = TRUE)
  received <- list()
  draw <- function(replicate) replicate
  environment(draw) <- baseenv()
  interrupted <- tryCatch({
    simulation_helpers$run_simulation_replications(
      6L, draw, list(), cluster = cluster,
      on_result = function(replicate, value) {
        received[[as.character(replicate)]] <<- value
        if (length(received) == 2L) {
          stop(structure(list(message = "test interruption", call = NULL),
                         class = c("interrupt", "condition")))
        }
      }
    )
    FALSE
  }, interrupt = function(condition) TRUE)
  testthat::expect_true(interrupted)
  testthat::expect_length(received, 2L)
  testthat::expect_identical(as.integer(names(received)), unlist(received,
                                                               use.names = FALSE))
  simulation_helpers$stop_simulation_cluster(cluster)
  cluster <- NULL
  testthat::expect_setequal(rownames(showConnections(all = TRUE)),
                            previous_connections)
})

scenario2_execution_fixture <- function(interrupt_at = NA_integer_) {
  execution <- new.env(parent = simulation_helpers)
  expressions <- as.list(parse(file.path(dirname(helper_path), "scenario2.R")))
  for (expression in expressions) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        is.call(expression[[3L]]) &&
        identical(expression[[3L]][[1L]], as.name("function"))) {
      eval(expression, execution)
    }
  }
  execution$compute_reference_target <- function(...) 0
  execution$options <- list(quick_run = TRUE, cores = 1L, replications = 4L)
  execution$design <- list(replications = 4L, treated_sizes = c(2L, 3L),
                            control_multiplier = 4L)
  execution$settings <- data.frame(setting_id = 1L, setting = "None",
                                   kappa_pi = 0, kappa_m = 0)
  execution$interrupt_at <- interrupt_at
  execution$interrupt_n1 <- 2L
  execution$run_replication <- function(setting, treated_count, replicate,
                                         target, design) {
    if (!is.na(interrupt_at) && replicate == interrupt_at &&
        treated_count == interrupt_n1) {
      stop(structure(list(message = "test interruption", call = NULL),
                     class = c("interrupt", "condition")))
    }
    make_result(setting, treated_count, replicate, "MEC-Cox (BART/Cox)",
                 target, design, estimate = 0.1, standard_error = 0.2)
  }
  environment(execution$run_replication) <- execution
  execution$expressions <- expressions
  execution
}

testthat::test_that("Scenario 2 returns partial results on interruption without files", {
  scratch <- tempfile("simulation interrupt ")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  previous_directory <- setwd(scratch)
  on.exit(setwd(previous_directory), add = TRUE)
  execution <- scenario2_execution_fixture(interrupt_at = 3L)
  messages <- character()
  output <- capture.output(answer <- withCallingHandlers(
    execution$run_scenario2(execution$design, execution$settings,
                            execution$options),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    }
  ))
  testthat::expect_identical(answer$status, "interrupted")
  testthat::expect_identical(answer$replications$replicate, 1:2)
  testthat::expect_equal(answer$summary$replications, 2L)
  testthat::expect_equal(answer$summary$successful, 2L)
  testthat::expect_identical(answer$metadata$execution$completed_replications, 2L)
  testthat::expect_equal(answer$metadata$execution$planned_replications, 8L)
  testthat::expect_equal(answer$targets$target, 0)
  testthat::expect_true(any(grepl("2 cells; 4 runs per cell = 8 datasets",
                                 messages, fixed = TRUE)))
  testthat::expect_true(any(grepl("Completed 1/4", messages, fixed = TRUE)))
  testthat::expect_true(any(grepl("interrupted: 2/8", messages, fixed = TRUE)))
  testthat::expect_length(list.files(scratch, all.files = TRUE, no.. = TRUE), 0L)
})

testthat::test_that("Scenario 2 exposes empty results without reporting after an interrupt", {
  execution <- scenario2_execution_fixture(interrupt_at = 1L)
  execution$build_simulation_report <- function(...) stop("must not build report")
  execution$show_simulation_report <- function(...) stop("must not show report")
  execution$plot_scenario2_results <- function(...) stop("must not plot results")
  result_index <- which(vapply(execution$expressions, function(expression) {
    is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      identical(expression[[2L]], as.name("scenario2_results"))
  }, logical(1)))
  for (expression in execution$expressions[result_index:length(execution$expressions)]) {
    invisible(capture.output(suppressMessages(eval(expression, execution))))
  }
  testthat::expect_identical(execution$scenario2_results$status, "interrupted")
  testthat::expect_equal(nrow(execution$scenario2_replications), 0L)
  testthat::expect_equal(nrow(execution$scenario2_summary), 0L)
  testthat::expect_named(execution$scenario2_replications,
                         c("setting", "kappa_pi", "kappa_m", "n1", "n0",
                           "replicate", "method", "target", "estimate",
                           "standard_error", "error"))
  testthat::expect_null(execution$scenario2_report)
  testthat::expect_identical(execution$scenario2_tables, list())
})

testthat::test_that("partial reports identify incomplete simulation results", {
  execution <- scenario2_execution_fixture(interrupt_at = 2L)
  invisible(capture.output(answer <- suppressMessages(execution$run_scenario2(
    execution$design, execution$settings, execution$options
  ))))
  report <- simulation_helpers$build_simulation_report(answer, "Scenario 2")
  html <- htmltools::renderTags(report$report)$html
  testthat::expect_match(html, "partial results", fixed = TRUE)
  testthat::expect_match(html, "interrupted; 1 of 8 planned datasets", fixed = TRUE)
})

testthat::test_that("interruption retains completed cells and the current partial cell", {
  execution <- scenario2_execution_fixture(interrupt_at = 3L)
  execution$interrupt_n1 <- 3L
  invisible(capture.output(answer <- suppressMessages(execution$run_scenario2(
    execution$design, execution$settings, execution$options
  ))))
  testthat::expect_identical(answer$status, "interrupted")
  testthat::expect_identical(answer$replications$n1, c(rep(2L, 4L), rep(3L, 2L)))
  testthat::expect_identical(answer$replications$replicate, c(1:4, 1:2))
  testthat::expect_equal(answer$summary$replications, c(4L, 2L))
  testthat::expect_identical(answer$metadata$execution$completed_replications, 6L)
})

testthat::test_that("partial Scenario 2 plots omit unstarted settings", {
  execution <- scenario2_execution_fixture(interrupt_at = 2L)
  invisible(capture.output(answer <- suppressMessages(execution$run_scenario2(
    execution$design, execution$settings, execution$options
  ))))
  settings <- rbind(execution$settings,
                     data.frame(setting_id = 2L, setting = "Not started",
                                 kappa_pi = 1, kappa_m = 2))
  grDevices::pdf(NULL, width = 14, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  testthat::expect_no_error(
    simulation_helpers$plot_scenario2_results(answer$summary, settings)
  )
})

testthat::test_that("unexpected Scenario 2 errors also retain prior results", {
  execution <- scenario2_execution_fixture()
  execution$run_replication <- function(setting, treated_count, replicate,
                                         target, design) {
    if (replicate == 3L) stop("test replication failure")
    make_result(setting, treated_count, replicate, "MEC-Cox (BART/Cox)",
                 target, design, estimate = 0.1, standard_error = 0.2)
  }
  environment(execution$run_replication) <- execution
  invisible(capture.output(answer <- suppressMessages(execution$run_scenario2(
    execution$design, execution$settings, execution$options
  ))))
  testthat::expect_identical(answer$status, "failed")
  testthat::expect_identical(answer$replications$replicate, 1:2)
  testthat::expect_match(answer$metadata$execution$condition,
                         "test replication failure", fixed = TRUE)
})

testthat::test_that("cluster cleanup tolerates a broken owned socket", {
  testthat::skip_on_cran()
  previous_connections <- rownames(showConnections(all = TRUE))
  cluster <- simulation_helpers$start_simulation_cluster(2L, character())
  on.exit(simulation_helpers$stop_simulation_cluster(cluster), add = TRUE)
  close(cluster[[1L]]$con)
  testthat::expect_no_error(simulation_helpers$stop_simulation_cluster(cluster))
  cluster <- NULL
  testthat::expect_setequal(rownames(showConnections(all = TRUE)),
                            previous_connections)
})
