# =============================================================================
# Script:       13_lepagella_volcanos.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Bin-level volcano plots for the two Lepagella muris MAGs
#               (Bin.44 and Bin.86) comparing Positive vs Negative fraction
#               gene expression
# Dependencies: data_structuring.R, plotting.R
# Input:        Bin_DEGs_08-05-2026.csv (from 05_bin_deg_analysis.R)
# Output:       A .svg of Bin.44 and Bin.86 volcanos. Figure 4 Panel C
# =============================================================================

# ---- Paths (update for your system) ----
fig_dir  <- "/path/to/figure_output_dir"
data_dir <- "/path/to/data"

library(ggplot2)
library(ggrepel)
library(patchwork)
library(svglite)

source("statistics/functions/plotting.R")

## ---- Load DEGs ----
DEGs <- read.csv(file.path(data_dir, "Bin_DEGs_08-05-2026.csv"))

## ---- Plotting function ----
plot_bin_volcano <- function(bin_id, DEGs, highlight_products = NULL) {

    bin_data <- DEGs %>%
        dplyr::filter(Bin == bin_id,
            contrast == "Positive - Negative") %>%
        mutate(
            direction = case_when(
                significant & estimate > 0 ~ "Up in Bound",
                significant & estimate < 0 ~ "Down in Bound",
                TRUE ~ "NS"),
            label = ifelse(product %in% highlight_products, product, NA))

    ggplot(bin_data, aes(x = estimate, y = -log10(p_adjust),
            colour = direction)) +
        geom_point(size = pt_small, alpha = 0.6) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed",
            linewidth = lw_sm, colour = "grey50") +
        geom_vline(xintercept = 0, linetype = "dashed",
            linewidth = lw_sm, colour = "grey50") +
        ggrepel::geom_label_repel(aes(label = label),
            size = base_sz / .pt, family = "Arial",
            max.overlaps = 20,
            na.rm = TRUE) +
        scale_colour_manual(values = c(
            "Up in Bound"   = col_pos,
            "Down in Bound" = col_neg,
            "NS"            = "grey70")) +
        labs(x = "Δ CLR(RNA)",
            y = expression(-log[10](FDR))) +
        theme_for_figures(show_labels = TRUE) +
        theme(legend.position = "none",
            plot.title = element_text(face = "bold.italic", size = base_sz),
            plot.margin = ggplot2::margin(0, 0, 0, 0, "mm"))
}

## ---- Plot ----
p_volc_44 <- plot_bin_volcano("Bin.44", DEGs)
p_volc_86 <- plot_bin_volcano("Bin.86", DEGs)

volc_combined <- wrap_plots(list(p_volc_44, p_volc_86), ncol = 2)

## ---- Save ----
panel_w <- 70.1
panel_h <- 45

ggsave(file.path(fig_dir, "bin_volcanos.svg"), volc_combined,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")
