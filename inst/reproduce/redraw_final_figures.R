# Explicit export only. Sourcing defines redraw_final_figures() and does not run it.
# Rscript redraw_final_figures.R --reference=/path/final_figure_reference.rds --output=/new/empty/directory
source_paths <- vapply(sys.frames(), function(frame)
  if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L], character(1))
source_paths <- rev(source_paths[nzchar(source_paths)])
command_paths <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
directories <- c(dirname(c(source_paths, command_paths)), ".", file.path("inst", "reproduce"),
                   system.file("reproduce", package = "mecCox"))
helper_paths <- file.path(directories, "final_figure_helpers.R")
helper_paths <- helper_paths[file.exists(helper_paths)]
if (!length(helper_paths)) stop("Keep the complete reproduce folder together.")
source(helper_paths[1L], local = TRUE)
if (!interactive() && !length(source_paths)) {
  arguments <- commandArgs(TRUE)
  reference <- grep("^--reference=.+", arguments, value = TRUE)
  output <- grep("^--output=.+", arguments, value = TRUE)
  if (length(arguments) != 2L || length(reference) != 1L || length(output) != 1L)
    stop("Use --reference=/path/final_figure_reference.rds --output=/new/empty/directory")
  redraw_final_figures(sub("^--reference=", "", reference), sub("^--output=", "", output))
}
