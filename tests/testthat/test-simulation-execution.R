# These scripts are installed as reproduction examples, outside the package API.
helper_path <- testthat::test_path("..", "..", "inst", "reproduce",
                                  "simulation_helpers.R")
if (!file.exists(helper_path)) {
  helper_path <- system.file("reproduce", "simulation_helpers.R",
                             package = "mecCox")
}
simulation_helpers <- new.env(parent = baseenv())
source(helper_path, local = simulation_helpers)

testthat::test_that("simulation arguments have useful defaults and reject mistakes", {
  parse_arguments <- simulation_helpers$parse_simulation_arguments
  defaults <- parse_arguments("example-output", character())
  testthat::expect_identical(
    defaults,
    list(quick_run = FALSE, output_directory = "example-output", cores = 20L)
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
