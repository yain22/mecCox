# Original learner environment for the supplementary breast-cancer tables.
.breast_study_versions <- c(brulee = "0.6.0", torch = "0.16.1",
                            ranger = "0.17.0", survival = "3.8-3")

.breast_study_library <- function() {
  file.path(tools::R_user_dir("mecCox", "cache"), "breast-cancer-20261009",
            paste0("R-", R.version$major, ".",
                   strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L]))
}

.breast_use_study_environment <- function() {
  versions <- .breast_study_versions
  matches <- function(package, library = NULL) {
    tryCatch(utils::packageVersion(package, lib.loc = library) ==
               package_version(versions[[package]]), error = function(e) FALSE)
  }
  loaded <- intersect(names(versions), loadedNamespaces())
  incompatible <- loaded[!vapply(loaded, function(package) {
    getNamespaceVersion(package) == package_version(versions[[package]])
  }, logical(1))]
  if (length(incompatible)) {
    stop("Restart R, then source breast_cancer.R before loading other packages. ",
         "Different study-package versions are already loaded: ",
         paste(incompatible, collapse = ", "), ".", call. = FALSE)
  }

  library <- .breast_study_library()
  if ("torch" %in% loaded &&
      !identical(normalizePath(getNamespaceInfo("torch", "path"), mustWork = FALSE),
                 normalizePath(file.path(library, "torch"), mustWork = FALSE)))
    stop("Restart R, then source breast_cancer.R before loading Torch.", call. = FALSE)
  old_libraries <- .libPaths()
  old_environment <- Sys.getenv(c("TORCH_INSTALL", "TORCH_LOAD", "TORCH_HOME",
    "TORCH_URL", "LANTERN_URL", "LANTERN_BASE_URL", "TORCH_COMMIT_SHA",
    "TORCH_VERIFY_LOAD", "CUDA"), unset = NA_character_)
  old_timeout <- getOption("timeout")
  restore <- function() {
    .libPaths(old_libraries)
    for (name in names(old_environment)) {
      if (is.na(old_environment[[name]])) Sys.unsetenv(name)
      else do.call(Sys.setenv, setNames(list(old_environment[[name]]), name))
    }
    options(timeout = old_timeout)
  }
  ready <- FALSE
  on.exit(if (!ready) restore(), add = TRUE)
  dir.create(library, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(library, .libPaths()))
  missing <- names(versions)[!vapply(names(versions), function(package) {
    matches(package, if (package == "survival") NULL else library)
  }, logical(1))]
  options(timeout = max(600, old_timeout))
  Sys.unsetenv(c("TORCH_LOAD", "TORCH_URL", "LANTERN_URL", "LANTERN_BASE_URL",
                "TORCH_COMMIT_SHA", "TORCH_VERIFY_LOAD"))
  Sys.setenv(TORCH_INSTALL = "0", CUDA = "cpu",
             TORCH_HOME = file.path(library, "torch"))
  if (length(missing)) {
    message("Installing the original breast-cancer learner versions in a separate library ...")
    repository <- "https://packagemanager.posit.co/cran/2025-10-25"
    type <- .Platform$pkgType
    available <- utils::available.packages(repos = repository, type = type)
    supported <- missing %in% rownames(available)
    supported[supported] <- vapply(missing[supported], function(package) {
      package_version(available[package, "Version"]) ==
        package_version(versions[[package]])
    }, logical(1))
    if (!all(supported)) {
      stop("The original study packages are unavailable for this R/platform combination. ",
           "The verified environment is R 4.5.1 on Windows. Required versions: ",
           paste(paste(names(versions), versions), collapse = ", "), ".",
           call. = FALSE)
    }
    utils::install.packages(missing, lib = library, repos = repository,
                             type = type, dependencies = NA)
  }
  if (!all(vapply(names(versions), matches, logical(1))))
    stop("Installation of the original breast-cancer package versions failed.",
         call. = FALSE)
  torch_directory <- file.path(library, "torch")
  runtime_version_file <- file.path(torch_directory, "build-version")
  if (!file.exists(runtime_version_file)) {
    message("Installing the original CPU Torch runtime ...")
    # Install in a child R session so the current session never loads a
    # partially installed Torch DLL. Its library path is explicit.
    setup <- tempfile(fileext = ".R")
    on.exit(unlink(setup), add = TRUE)
    writeLines(c(
      paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
      "Sys.setenv(TORCH_INSTALL = '0', TORCH_LOAD = '0', TORCH_VERIFY_LOAD = 'FALSE', CUDA = 'cpu')",
      "options(timeout = 600)",
      "torch::install_torch()"), setup)
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c("--vanilla", shQuote(setup)))
    if (status != 0L)
      stop("The original CPU Torch runtime could not be installed.", call. = FALSE)
  }
  if (.Platform$OS.type == "windows" &&
      (!file.exists(runtime_version_file) ||
       !identical(trimws(readLines(runtime_version_file, warn = FALSE)), "2.7.1+cpu")))
    stop("The breast-cancer study requires its original Torch 2.7.1 CPU runtime.",
         call. = FALSE)
  torch::torch_tensor(0)
  ready <- TRUE
  invisible(list(library = library, versions = versions, restore = restore))
}
