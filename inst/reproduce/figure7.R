# In RStudio, click Source to display Figure 7 in the Plots pane.
# From the installed package, run:
#   source(system.file("reproduce", "figure7.R", package = "mecCox"))

# Choose one option: keep its settings line uncommented.
# Manuscript figure: plot the archived estimates (1,000 replications per target).
figure7_mode <- "archived"

# Publication replication: regenerate all datasets and refit the models (takes time).
# figure7_mode <- "full"   # 1,000 replications per target, ten workers

# Quick learning: run the same eight panels with fewer replications.
# figure7_mode <- "quick"  # 200 replications per target, ten workers

figure7_study_dir <- function() {
  # Locate the study beside the sourced script or the Rscript entry point.
  source_files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
  }, character(1))
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  scripts <- c(rev(source_files[nzchar(source_files)]), sub("^--file=", "", file_arg))
  candidates <- c(file.path(dirname(scripts), "figure7_study"),
                  "figure7_study", file.path("inst", "reproduce", "figure7_study"),
                  system.file("reproduce", "figure7_study", package = "mecCox"))
  found <- candidates[nzchar(candidates) &
                        file.exists(file.path(candidates, "specification.rds"))]
  if (!length(found)) {
    stop("Install the current package or keep figure7_study beside figure7.R.")
  }
  normalizePath(found[1L])
}

simulate_figure7 <- function(study_dir, run_dir, replications) {
  inputs <- c("run_simulation.R", "dgp.R", "estimator.R", "original_project_functions.R",
              "project_code_verification.csv", "validate_results.R")
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  previous_dir <- setwd(run_dir)
  on.exit(setwd(previous_dir))
  # Keep checkpoints tied to the scripts that generated them.
  existing <- inputs[file.exists(inputs)]
  if (any(unname(tools::md5sum(existing)) !=
          unname(tools::md5sum(file.path(study_dir, existing))))) {
    stop("Existing simulation files differ from the study scripts. Use a fresh output directory.")
  }
  missing <- setdiff(inputs, existing)
  if (length(missing) && !all(file.copy(file.path(study_dir, missing), ".", overwrite = FALSE))) {
    stop("Could not copy the Figure 7 simulation scripts.")
  }
  rscript <- file.path(R.home("bin"), "Rscript")
  if (system2(rscript, c("--vanilla", "run_simulation.R", as.character(replications), "10")) != 0L ||
      system2(rscript, c("--vanilla", "validate_results.R")) != 0L) {
    stop("Figure 7 simulation or validation failed.")
  }
  normalizePath("results")
}

read_figure7_results <- function(results_dir, replications) {
  raw <- readRDS(file.path(results_dir, "replications.rds"))
  spec <- readRDS(file.path(results_dir, "specification.rds"))
  methods <- c("Oracle", "Naive", "IPW", "MEC1", "MEC5", "MEC10", "MEC20")
  raw <- raw[raw$method %in% methods, ]
  stopifnot(spec$n1 == 200L, spec$n0 == 400L, spec$B == replications,
            isTRUE(all.equal(spec$configuration$tau, log(8))),
            all(raw$n1 == spec$n1), all(raw$n0 == spec$n0),
            setequal(unique(raw$method), methods),
            setequal(unique(raw$true_HR), c(1, 0.70)),
            setequal(unique(raw$PS), c("Correct", "Misspecified")),
            setequal(unique(raw$OR), c("Correct", "Misspecified")),
            setequal(unique(raw$replication), seq_len(spec$B)),
            !anyDuplicated(raw[, c("true_HR", "replication", "PS", "OR", "method")]),
            nrow(raw) == 8L * length(methods) * spec$B,
            all(with(raw, table(true_HR, PS, OR, method)) == spec$B))
  raw$method <- factor(raw$method, levels = methods)
  raw[raw$converged & is.finite(raw$HR) & raw$HR > 0, ]
}

draw_figure7 <- function(raw) {
  methods <- levels(raw$method)
  labels <- c("Oracle", "Naive", "IPW", "MEC-Cox\n(1)", "MEC-Cox\n(5)",
              "MEC-Cox\n(10)", "MEC-Cox\n(20)")
  colors <- c("#5C7836", "#626262", "#202020", "#0072B2", "#8B4789",
              "#D55E00", "#007F73")
  fills <- grDevices::adjustcolor(colors, alpha.f = 0.22)
  scenario <- function(ps, or) ifelse(ps == "Correct",
    ifelse(or == "Correct", 1L, 2L), ifelse(or == "Correct", 3L, 4L))
  titles <- c("PS correct, OR correct", "PS correct, OR incorrect",
              "PS incorrect, OR correct", "PS incorrect, OR incorrect")

  par(mfrow = c(4, 2), cex = 1, mar = c(3.65, 4.55, 2.2, 0.7),
      oma = c(0, 0, 0, 0), mgp = c(2.75, 0.7, 0), family = "sans", las = 1)
  panel <- 0L
  for (truth in c(1, 0.70)) {
    yrange <- range(c(truth, raw$HR[raw$true_HR == truth]))
    padding <- max(0.02, diff(yrange) * 0.05)
    yrange <- c(max(0, yrange[1] - padding), yrange[2] + padding)
    for (s in 1:4) {
      panel <- panel + 1L
      current <- raw[raw$true_HR == truth & scenario(raw$PS, raw$OR) == s, ]
      values <- lapply(methods, function(method)
        current$HR[current$method == method])
      stopifnot(all(lengths(values) > 0L))
      boxplot(values, names = labels, xaxt = "n", ylim = yrange, outline = TRUE,
              col = fills, border = colors, xlab = "",
              ylab = "Estimated ATT hazard ratio", main = NULL,
              cex.axis = 0.83, cex.lab = 0.90, boxwex = 0.62,
              pars = list(outpch = 16, outcex = 0.23, medlwd = 1.5,
                          whisklty = 1, staplewex = 0.6))
      title(main = paste0("(", letters[panel], ") ",
                         if (truth == 1) "Null: " else "Alternative: ", titles[s]),
            adj = 0, line = 0.8, cex.main = 0.94, font.main = 2)
      axis(1, at = seq_along(labels), labels = FALSE)
      mtext(labels, side = 1, at = seq_along(labels), line = 1.0, cex = 0.70)
      abline(h = truth, col = "#555555", lty = 2, lwd = 1.1)
      box()
    }
  }
}

# Lay out the eight panels at manuscript size. The Plots pane and Zoom can
# then scale the complete figure to fit, even in a short or narrow pane.
capture_figure7 <- function(raw) {
  capture <- ragg::agg_capture(width = 10.6, height = 12.6, units = "in",
                               res = 180, background = "white")
  on.exit(grDevices::dev.off())
  draw_figure7(raw)
  capture(native = TRUE)
}

save_figure7_pdf <- function(raw, file) {
  grDevices::pdf(file, width = 10.6, height = 12.6, family = "Helvetica",
                 useDingbats = FALSE, title = "MEC-Cox: null and alternative point estimates")
  on.exit(grDevices::dev.off())
  draw_figure7(raw)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) && !identical(args, "--simulate")) {
    stop("Usage: Rscript figure7.R [--simulate]")
  }
  if (length(args)) figure7_mode <- "full"
}
stopifnot(figure7_mode %in% c("archived", "full", "quick"))
study_dir <- figure7_study_dir()
output_dir <- "figure7-output"
replications <- if (figure7_mode == "quick") 200L else 1000L
results_dir <- study_dir
if (figure7_mode != "archived") {
  run_folder <- if (figure7_mode == "quick") "simulation-200" else "simulation"
  results_dir <- simulate_figure7(study_dir, file.path(output_dir, run_folder), replications)
}
figure7_estimates <- read_figure7_results(results_dir, replications)
figure7_summary <- stats::aggregate(list(median_HR = figure7_estimates$HR),
                                    figure7_estimates[c("true_HR", "PS", "OR", "method")],
                                    stats::median)
if (interactive() || grDevices::dev.cur() != 1L) {
  figure7_image <- capture_figure7(figure7_estimates)
  grid::grid.newpage()
  grid::grid.raster(figure7_image)
} else {
  dir.create(output_dir, showWarnings = FALSE)
  figure_file <- file.path(output_dir, if (figure7_mode == "quick") "figure7-quick.pdf" else "figure7.pdf")
  save_figure7_pdf(figure7_estimates, figure_file)
  message("Wrote ", figure_file)
}
