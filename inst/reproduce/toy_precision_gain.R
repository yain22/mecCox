# Reproduce the toy example where prognostic balancing improves precision.
# Run from a checkout (installing mecCox itself is not required):
#   Rscript inst/reproduce/toy_precision_gain.R
# Or source this file in R / RStudio after editing the settings below.
# --quick uses 20 datasets for a code-path check, not a manuscript result.
# --replications=N overrides the count; --output-dir=PATH optionally saves
# CSV tables, the PDF/PNG figure, and sessionInfo. No files are saved by default.
# Results remain in toy_precision_gain_results.

quick_run <- FALSE
replications <- 10000L
seed <- 20261007L
output_dir <- NULL

# Resolve the adjacent helper for Rscript, source(), or an editor selection.
source_files <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
}, character(1))
source_files <- rev(source_files[nzchar(source_files)])
script_option <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_files <- c(source_files, sub("^--file=", "", script_option))
helper_candidates <- c(
  file.path(dirname(script_files), "toy_helpers.R"),
  "toy_helpers.R",
  file.path("inst", "reproduce", "toy_helpers.R"),
  file.path("mecCox", "inst", "reproduce", "toy_helpers.R"),
  system.file("reproduce", "toy_helpers.R", package = "mecCox")
)
helper_candidates <- unique(helper_candidates[
  nzchar(helper_candidates) & file.exists(helper_candidates)
])
if (!length(helper_candidates)) {
  stop("Cannot find toy_helpers.R. Keep it beside this script, run from the ",
       "repository root, or install the current mecCox package.", call. = FALSE)
}
source(helper_candidates[1L], local = TRUE)

# Do not consume an enclosing Rscript's arguments when this file is sourced.
arguments <- if (interactive() || length(source_files)) character() else commandArgs(trailingOnly = TRUE)
toy_options <- parse_toy_arguments(arguments, quick_run = quick_run,
                                   replications = replications, seed = seed,
                                   output_dir = output_dir)
if (toy_options$quick_run) {
  message("Quick check requested: ", toy_options$replications,
          " datasets; manuscript design and estimators retained.")
}
toy_precision_gain_results <- run_toy_example(
  "toy_precision_gain", replications = toy_options$replications,
  seed = toy_options$seed, output_dir = toy_options$output_dir
)

