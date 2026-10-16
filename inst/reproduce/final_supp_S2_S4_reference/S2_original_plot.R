setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(patchwork)

# ============================================================================
# Common settings (same as Scenario 1)
# ============================================================================

base_dir <- "."

method_labels <- c(
  "IPW-Naive"  = "Naive",
  "IPW-LW"     = "Robust sandwich",
  "ML-IPW-LW"  = "Robust sandwich (flexible-PS)",
  "IPW-Shu"    = "Corrected sandwich",
  "MEC-KL"     = "MEC-Cox (KL)",
  "MEC-Hell."  = "MEC-Cox (Hellinger)",
  "MEC-EL"     = "MEC-Cox (EL)",
  "MEC-Renyi"  = "MEC-Cox (Rényi)"
)

method_order <- c(
  "Naive", "Robust sandwich", "Corrected sandwich",
  "Robust sandwich (flexible-PS)",
  "MEC-Cox (KL)", "MEC-Cox (Hellinger)", "MEC-Cox (EL)", "MEC-Cox (Rényi)"
)

method_colors <- c(
  "Naive"                          = "gray25",
  "Robust sandwich"                = "#0072B2",
  "Corrected sandwich"             = "#009E73",
  "Robust sandwich (flexible-PS)"  = "#56B4E9",
  "MEC-Cox (KL)"                   = "#D55E00",
  "MEC-Cox (Hellinger)"            = "#CC79A7",
  "MEC-Cox (EL)"                   = "#8B008B",
  "MEC-Cox (Rényi)"                = "#E69F00"
)

method_shapes <- c(
  "Naive"                          = 4,
  "Robust sandwich"                = 1,
  "Corrected sandwich"             = 2,
  "Robust sandwich (flexible-PS)"  = 0,
  "MEC-Cox (KL)"                   = 16,
  "MEC-Cox (Hellinger)"            = 17,
  "MEC-Cox (EL)"                   = 15,
  "MEC-Cox (Rényi)"                = 18
)

method_linetypes <- c(
  "Naive"                          = "dotted",
  "Robust sandwich"                = "dashed",
  "Corrected sandwich"             = "longdash",
  "Robust sandwich (flexible-PS)"  = "dotdash",
  "MEC-Cox (KL)"                   = "solid",
  "MEC-Cox (Hellinger)"            = "solid",
  "MEC-Cox (EL)"                   = "solid",
  "MEC-Cox (Rényi)"                = "solid"
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

# Helper: build 2x3 figure (2 settings x 3 metrics)
build_2x3 <- function(data, setting_var,
                      setting_titles = NULL,
                      panel_labels = paste0("(", letters[1:6], ")"),
                      include_methods = NULL) {
  settings <- levels(data$.setting)
  if (is.null(setting_titles)) setting_titles <- settings

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

  combined <- (panels[[1]] + panels[[2]] + panels[[3]]) /
              (panels[[4]] + panels[[5]] + panels[[6]]) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal")

  combined
}


# ============================================================================
# Scenario 2: Linear PS & OR, n1:n0 = 1:4, varying M
#             Fitting by GLM-lasso (PS) and Cox-lasso (OR)
# ============================================================================

# All available files (choose 2 by commenting/uncommenting)
# M = 50
# M = 100
# M = 150
# M = 200

s2_files <- c(
  # --- Pick 2 of the files below ---
  "mcx_lp_gen_pt-lin_ot-lph_pd-2_od-2_r-4_p-50_mps-glmlas_pslam-1se_bas-lcoxlp_blam-1se_g-kl-hel-el-renyi05_R-1000_0521.RData",
  "mcx_lp_gen_pt-lin_ot-lph_pd-2_od-2_r-4_p-200_mps-glmlas_pslam-1se_bas-lcoxlp_blam-1se_g-kl-hel-el-renyi05_R-1000_0522.RData"
  # "mcx_lp_gen_pt-lin_ot-lph_pd-2_od-2_r-4_p-100_mps-glmlas_pslam-1se_bas-lcoxlp_blam-1se_g-kl-hel-el-renyi05_R-1000_0521.RData"
  # "mcx_lp_gen_pt-lin_ot-lph_pd-2_od-2_r-4_p-150_mps-glmlas_pslam-1se_bas-lcoxlp_blam-1se_g-kl-hel-el-renyi05_R-1000_0524.RData"
)

s2_labels <- c(
  # --- Match labels to the 2 files chosen above ---
  "M = 50",
  "M = 200"
  # "M = 100"
  # "M = 150"
)

# Choose which methods to include (comment/uncomment as needed)
#s2_methods <- c(
#  "Naive", "Robust sandwich", "Corrected sandwich",
#  "Robust sandwich (flexible-PS)",
#  "MEC-Cox (KL)", "MEC-Cox (Hellinger)", "MEC-Cox (EL)", "MEC-Cox (Rényi)"
#)
# Example: include only KL among MEC variants
 s2_methods <- c("Naive", "Robust sandwich", "Corrected sandwich", "MEC-Cox (KL)")

s2 <- load_and_combine(s2_files, "M_label", s2_labels,
                        include_methods = s2_methods)

s2_titles <- s2_labels

fig_s2 <- build_2x3(s2, "M_label",
                     setting_titles = s2_titles,
                     include_methods = s2_methods)

fig_s2

# ggsave("fig_scenario2.pdf", fig_s2, width = 16, height = 12, dpi = 300)
# ggsave("fig_scenario2.png", fig_s2, width = 16, height = 12, dpi = 300)
