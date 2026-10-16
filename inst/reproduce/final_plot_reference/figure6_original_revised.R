setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(patchwork)

# ============================================================================
# Common settings
# ============================================================================

base_dir <- "."

method_labels <- c(
  "IPW-Naive"      = "Naive",
  "IPW-LW"         = "Robust sandwich",
  "IPW-Shu"        = "Corrected sandwich",
  "MEC-Cox (cox)"  = "MEC-Cox (KL-BART(PS)-Cox(OR))",
  "MEC-Cox (rsf)"  = "MEC-Cox (KL-BART(PS)-RSF(OR))"
)

method_order <- c(
  "Naive", "Robust sandwich", "Corrected sandwich",
  "MEC-Cox (KL-BART(PS)-Cox(OR))", "MEC-Cox (KL-BART(PS)-RSF(OR))"
)

method_colors <- c(
  "Naive"              = "gray25",
  "Robust sandwich"    = "#0072B2",
  "Corrected sandwich" = "#009E73",
  "MEC-Cox (KL-BART(PS)-Cox(OR))"      = "#D55E00",
  "MEC-Cox (KL-BART(PS)-RSF(OR))"      = "blue"
)

method_shapes <- c(
  "Naive"              = 4,
  "Robust sandwich"    = 1,
  "Corrected sandwich" = 2,
  "MEC-Cox (KL-BART(PS)-Cox(OR))"      = 16,
  "MEC-Cox (KL-BART(PS)-RSF(OR))"      = 17
)

method_linetypes <- c(
  "Naive"              = "dotted",
  "Robust sandwich"    = "dashed",
  "Corrected sandwich" = "longdash",
  "MEC-Cox (KL-BART(PS)-Cox(OR))"      = "solid",
  "MEC-Cox (KL-BART(PS)-RSF(OR))"      = "solid"
)

theme_pub <- theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.5, "cm"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 13, face = "bold", hjust = 0)
  )

# Helper: load multiple .RData files and combine plot_dat with a setting label
load_and_combine <- function(files, setting_var, setting_labels,
                             include_methods = method_order) {
  all_dat <- list()
  for (i in seq_along(files)) {
    e <- new.env()
    load(file.path(base_dir, files[i]), envir = e)
    dat <- e$plot_dat
    dat[[setting_var]] <- setting_labels[i]
    all_dat[[i]] <- dat
  }
  combined <- bind_rows(all_dat) |>
    mutate(
      Method_Label = factor(method_labels[as.character(Method_Short)],
                            levels = method_order)
    ) |>
    filter(Method_Label %in% include_methods) |>
    mutate(
      Method_Label = droplevels(Method_Label),
      .setting = factor(.data[[setting_var]], levels = setting_labels)
    )
  combined
}

# Helper: create one panel
make_panel <- function(data, metric, ylab, title, hline = NULL,
                       include_methods = NULL) {
  cols <- method_colors
  shps <- method_shapes
  lts  <- method_linetypes
  if (!is.null(include_methods)) {
    cols <- cols[names(cols) %in% include_methods]
    shps <- shps[names(shps) %in% include_methods]
    lts  <- lts[names(lts)  %in% include_methods]
  }

  p <- ggplot(data, aes(x = n1, y = .data[[metric]],
                         color = Method_Label, shape = Method_Label,
                         linetype = Method_Label)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.2) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shps) +
    scale_linetype_manual(values = lts) +
    labs(x = expression(n[1]), y = ylab, title = title) +
    theme_pub

  if (!is.null(hline)) {
    p <- p +
      geom_hline(yintercept = hline, linetype = "dashed",
                 color = "gray30", linewidth = 0.5) +
      scale_y_continuous(breaks = function(lims) {
        b <- labeling::extended(lims[1], lims[2], m = 5)
        sort(unique(c(b, hline)))
      })
  }
  p
}

# Helper: build Nx3 figure (N settings x 3 metrics)
build_Nx3 <- function(data, setting_var,
                      setting_titles = NULL,
                      panel_labels = NULL,
                      include_methods = NULL) {
  settings <- levels(data$.setting)
  n_settings <- length(settings)
  if (is.null(setting_titles)) setting_titles <- settings
  if (is.null(panel_labels)) {
    panel_labels <- paste0("(", letters[seq_len(n_settings * 3)], ")")
  }

  used_methods <- if (!is.null(include_methods)) include_methods
                  else levels(data$Method_Label)

  panels <- list()
  idx <- 1
  for (i in seq_along(settings)) {
    d <- data |> filter(.setting == settings[i])
    stitle <- setting_titles[i]

    panels[[idx]] <- make_panel(
      d, "Coverage_95", "Coverage",
      paste0(panel_labels[idx], " Coverage: ", stitle),
      hline = 0.95, include_methods = used_methods
    )
    idx <- idx + 1

    panels[[idx]] <- make_panel(
      d, "Bias", "Bias",
      paste0(panel_labels[idx], " Bias: ", stitle),
      hline = 0, include_methods = used_methods
    )
    idx <- idx + 1

    panels[[idx]] <- make_panel(
      d, "RMSE", "RMSE",
      paste0(panel_labels[idx], " RMSE: ", stitle),
      include_methods = used_methods
    )
    idx <- idx + 1
  }

  # Stack rows dynamically: each row is 3 panels for one setting
  rows <- list()
  for (i in seq_len(n_settings)) {
    row_start <- (i - 1) * 3 + 1
    rows[[i]] <- panels[[row_start]] + panels[[row_start + 1]] + panels[[row_start + 2]]
  }
  combined <- Reduce(`/`, rows) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal") &
    guides(
      color    = guide_legend(nrow = 2, byrow = TRUE),
      shape    = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 2, byrow = TRUE)
    )

  combined
}


# ============================================================================
# Scenario: Nonlinear PS & outcome, BART PS, Cox-RSF outcome, KL divergence
# ============================================================================

# --- Diagnostic: check available Method_Short values ---
# Run this block first to verify method names in your data
e_check <- new.env()
load(file.path(base_dir,
  "mcx_lmk_pt-nl_ot-nlph_pd-1_od-2_r-4_p-10_mps-bart_mor-cox-rsf_n_land-5_Kfold-10_g-kl_R-1000_0529.RData"),
  envir = e_check)
cat("Available Method_Short values:\n")
print(unique(e_check$plot_dat$Method_Short))
cat("\nAvailable columns:\n")
print(names(e_check$plot_dat))
rm(e_check)

# --- Three selected files ---
s1_files <- c(
  "mcx_lmk_pt-nl_ot-nlph_pd-0_od-0_r-4_p-10_mps-bart_mor-cox-rsf_n_land-5_Kfold-10_g-kl_R-1000_0529.RData",
  "mcx_lmk_pt-nl_ot-nlph_pd-1_od-2_r-4_p-10_mps-bart_mor-cox-rsf_n_land-5_Kfold-10_g-kl_R-1000_0529.RData",
  "mcx_lmk_pt-nl_ot-nlph_pd-2_od-5_r-4_p-10_mps-bart_mor-cox-rsf_n_land-5_Kfold-10_g-kl_R-1000_0529.RData"
)

s1_labels <- c(
  "Nonlinearity: none",
  "Nonlinearity: mild",
  "Nonlinearity: severe"
)

# Methods to include (adjust after checking diagnostic output above)
s1_methods <- c(
  "Naive", "Robust sandwich", "Corrected sandwich",
  "MEC-Cox (KL-BART(PS)-Cox(OR))", "MEC-Cox (KL-BART(PS)-RSF(OR))"
  #"MEC-Cox (KL-RSF)"
)

s1 <- load_and_combine(s1_files, "nonlinearity", s1_labels,
                        include_methods = s1_methods)

s1_titles <- s1_labels

fig_s1 <- build_Nx3(s1, "nonlinearity",
                     setting_titles = s1_titles,
                     include_methods = s1_methods)

fig_s1

# ggsave("fig_scenario_nonlinear_bart.pdf", fig_s1, width = 12, height = 14, dpi = 300)
# ggsave("fig_scenario_nonlinear_bart.png", fig_s1, width = 12, height = 14, dpi = 300)
