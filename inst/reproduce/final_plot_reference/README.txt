Original Figure 5/6 rendering sources

The two complete revised plotting sources are copied unchanged. Do not source
them directly: their historical data-loading/editor commands are preserved as
provenance only. final_figure_helpers.R evaluates exactly these eight original
assignments after checking SHA256.csv:
  method_labels, method_order, method_colors, method_shapes, method_linetypes,
  theme_pub, make_panel, build_Nx3.
No scientific results are included in this directory.

Final rendering environment (distinct from the fitting dependency snapshot):
  R 4.5.1; ggplot2 4.0.0; patchwork 1.3.2; gtable 0.3.6.
The fitting bundle records ggplot2 4.0.3; that audited fitting snapshot remains
unchanged. Use a separate R library/session for the pinned rendering versions.
The renderer refuses other versions of the three pinned plotting packages.
It also requires dplyr, labeling, digest and grid. Exact installed versions,
R/session/device environment and any ragg version are recorded on export.

Both main figures use the published 12 x 9 inch device, 300 dpi PNG and PDF,
white background. The original source's commented 12 x 14 example is stale;
the supplied publication PDFs establish 12 x 9. All original fonts, colours,
markers, lines, legends, panel titles and automatic finite data ranges remain
unchanged. No artificial axis truncation or numerical-value clipping is added.
The scripts leave font-family selection at the original device default (PDF
Helvetica/sans on the reference R/Windows installation). Byte or pixel identity
across different operating systems/devices/fonts is not assumed. Final release
verification compares the rendered artifacts in the pinned reference environment.

The public interactive drivers use an equivalent refactored plotting builder.
This explicit final redraw path uses the exact original plotting expressions.
Complete main results passed independent audit; public summaries agreed
exactly, and both PNGs and rendered PDF pages matched the accepted figures
in the recorded rendering environment.
