# =============================================================================
# Script:       07b_fraction_rf_plots.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Plots AUC boxplot from Kat Coyte's fraction-classification
#               random forest (07a). Binary RF predicting whether a MAG's
#               relative expression profile came from IgA-positive or
#               IgA-negative fraction.
# Dependencies: plotting.R
# Input:        accuracy_data.csv (output of 07a, https://doi.org/10.48420/33951214)
#               auc_summary.csv (output of 07a, https://doi.org/10.48420/33951214)
# Output:       A .svg of fraction RF AUC boxplot. Figure 2 panel C
# =============================================================================

# ---- Paths (update for your system) ----
fig_dir  <- "/path/to/figure_output_dir"
data_dir <- "/path/to/data"

library(ggplot2)
library(svglite)

source("statistics/functions/plotting.R")

## ---- Load data ----
auc_data    <- read.csv(file.path(data_dir, "accuracy_data.csv"))
auc_summary <- read.csv(file.path(data_dir, "auc_summary.csv"))

## ---- AUC boxplot ----
auc_box <- ggplot(auc_data, aes(x = Type, y = accuracy, fill = Type)) +
    geom_jitter(width = 0.1, height = 0, size = pt_small,
        alpha = 0.25, color = "grey75") +
    geom_boxplot(width = 0.5, linewidth = lw_sm, alpha = 0.8,
        colour = "black", outlier.shape = NA) +
    scale_fill_manual(values = c("DNA" = "#EE7733", "RNA" = "#117733")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
        color = "grey50", linewidth = lw_sm) +
    theme_for_figures(show_labels = TRUE) +
    labs(x = "Feature type", y = "AUC")
auc_box

## ---- Save ----
panel_w <- 41.5
panel_h <- 46.2

ggsave(file.path(fig_dir, "AUC_boxplot.svg"), auc_box,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")
