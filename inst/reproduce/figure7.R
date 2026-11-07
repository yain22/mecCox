# Legacy entry point for the double-robustness study, now manuscript Figure 5.
# For the current manuscript Figure 7, use scenario2.R.
legacy_entry <- function() {
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
  }, character(1))
  args <- commandArgs(FALSE)
  script <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  candidates <- c(file.path(dirname(c(rev(files[nzchar(files)]), script)), "figure5.R"),
                  file.path("inst", "reproduce", "figure5.R"),
                  system.file("reproduce", "figure5.R", package = "mecCox"))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("Cannot locate figure5.R beside the legacy script.")
  normalizePath(candidates[1L])
}
message("Double robustness is now Figure 5. Use scenario2.R for the current Figure 7.")
if (!interactive() && sys.nframe() == 0L) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c(shQuote(legacy_entry()), shQuote(commandArgs(trailingOnly = TRUE))))
  if (status != 0L) quit(status = status)
} else {
  source(legacy_entry(), local = TRUE)
  figure7_estimates <- figure5_estimates
  figure7_summary <- figure5_summary
}
