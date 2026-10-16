# Supplemental Figures S2-S4 redraw from the saved complete reference.
# Source supplement_S2_S4_helpers.R first, then explicitly call the redraw.
# Sourcing defines the function only; redraw performs no model or target fits.

redraw_supplement_S2_S4 <- function(reference_path, reference_dir, output_directory) {
  bundle <- load_supp234_bundle(reference_dir)
  reference <- read_supp234_reference(reference_path, bundle)
  plot_locale <- supp234_plot_locale()
  on.exit(supp234_restore_plot_locale(plot_locale), add = TRUE)
  built <- build_supp234_figures(reference, bundle)
  out <- supp234_new_directory(output_directory)
  for (fig in names(built$figures)) {
    for (format in c("png", "pdf")) ggplot2::ggsave(
      file.path(out, paste0(fig, "_final_fresh1000.", format)), built$figures[[fig]],
      width = 12, height = 7, units = "in", dpi = 300, bg = "white")
    utils::write.csv(built$plot_data[[fig]], file.path(out, paste0(fig, "_plot_data.csv")), row.names = FALSE)
  }
  metadata <- list(status = "RENDERED_PENDING_PIXEL_AND_VISUAL_QA", datasets = 30000L,
    method_rows = 240000L,
    plotted_method_cells = sum(vapply(built$plot_data, nrow, integer(1))),
    eligible_plotted_method_cells = sum(vapply(built$plot_data, function(d) sum(!is.na(d$Trials)), integer(1))),
    possible_plotted_method_cells = 150L, width_in = 12, height_in = 7,
    dpi = 300, plotting_LC_CTYPE = plot_locale, manifest_md5 = bundle$provenance$manifest_md5,
    reference_sha256 = supp234_sha(reference_path),
    render_versions = vapply(c("ggplot2", "patchwork", "gtable", "dplyr", "tidyr"), function(p)
      as.character(utils::packageVersion(p)), character(1)),
    session = capture.output(utils::sessionInfo()), capabilities = capabilities(),
    pdf_device_defaults = grDevices::pdf.options(),
    ragg_version = if (requireNamespace("ragg", quietly = TRUE)) as.character(utils::packageVersion("ragg")) else NA_character_,
    note = "Original 4/7/4 curves and six panels per figure; exact original finite Estimate/SE primary summary.")
  saveRDS(metadata, file.path(out, "render_metadata.rds"), version = 3)
  supp234_write_hashes(out)
  invisible(built)
}
