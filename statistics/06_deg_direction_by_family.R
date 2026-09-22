# =============================================================================
# Script:       06_deg_direction_by_family.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Summarises DEG direction (up/down in Positive vs Negative)
#               by bacterial family as a diverging bar plot with sqrt-scaled
#               axes
# Dependencies: data_structuring.R, plotting.R
# Input:        DEGs (from 05_bin_deg_analysis.R or cached CSV)
#               Taxa, Bin_map (from data_structuring.R)
# Output:       A .svg of family diverging bar plot. Figure 2 panel B
# =============================================================================

# ---- Paths (update for your system) ----
fig_dir  <- "/path/to/figure_output_dir"
data_dir <- "/path/to/data"

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(svglite)

source("statistics/functions/plotting.R")
source("statistics/functions/data_structuring.R")

## ---- Load DEGs ----
DEGs <- read.csv(file.path(data_dir, "Bin_DEGs_08-05-2026.csv"))

## ---- Data structuring ----
gene_direction <- DEGs %>%
    dplyr::filter(contrast == "Positive - Negative") %>%
    dplyr::distinct(Bin, Geneid, significant, estimate) %>%
    dplyr::mutate(direction = dplyr::case_when(
        significant & estimate > 0 ~ "Up",
        significant & estimate < 0 ~ "Down",
        TRUE ~ "Unchanged"))

gene_direction_family <- gene_direction %>%
    dplyr::left_join(
        Bin_map %>%
            dplyr::mutate(Bin = paste0("Bin.", as.integer(Bin))) %>%
            dplyr::select(Bin, Identifier),
        by = "Bin") %>%
    dplyr::left_join(
        Taxa %>% tibble::rownames_to_column("Identifier") %>%
            dplyr::select(Identifier, Family) %>%
            dplyr::mutate(Family = gsub("f__", "", Family)),
        by = "Identifier") %>%
    dplyr::filter(!is.na(Family))

family_bin_counts <- gene_direction_family %>%
    dplyr::distinct(Family, Bin) %>%
    dplyr::count(Family, name = "n_bins")

family_direction_summary <- gene_direction_family %>%
    dplyr::count(Family, direction) %>%
    dplyr::group_by(Family) %>%
    dplyr::mutate(
        total = sum(n),
        pct   = n / total * 100) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(family_bin_counts, by = "Family") %>%
    dplyr::mutate(
        pct_signed = dplyr::if_else(direction == "Down", -pct, pct),
        pct_sqrt   = sign(pct_signed) * sqrt(abs(pct_signed)))

family_plot_data <- family_direction_summary %>%
    dplyr::filter(direction != "Unchanged")

family_order <- family_plot_data %>%
    tidyr::pivot_wider(id_cols = Family, names_from = direction,
        values_from = pct, values_fill = 0) %>%
    dplyr::mutate(diff = Down - Up) %>%
    dplyr::arrange(desc(diff)) %>%
    dplyr::pull(Family)

family_order <- union(family_order, unique(family_plot_data$Family))

## ---- Plot ----
sqrt_breaks <- c(-50, -40, -30, -20, -10, 0, 10)
sqrt_break_pos <- sign(sqrt_breaks) * sqrt(abs(sqrt_breaks))

family_direction_plot <- ggplot(family_plot_data,
        aes(x = pct_sqrt,
            y = factor(Family, levels = family_order),
            fill = direction)) +
    geom_col() +
    geom_vline(xintercept = 0, colour = "black", linewidth = lw_axis) +
    scale_fill_manual(name = NULL,
        values = c("Up" = col_pos, "Down" = col_neg)) +
    scale_x_continuous(
        breaks = sqrt_break_pos,
        labels = paste0(abs(sqrt_breaks), "%"),
        limits = c(-6.75, 3.5),
        expand = expansion(mult = c(0.02, 0.02))) +
    coord_flip(clip = "off") +
    theme_for_figures(show_labels = TRUE) +
    theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = "Percent differentially\nexpressed genes", y = NULL)
family_direction_plot

## ---- Save ----
panel_w <- 83
panel_h <- 46.2

ggsave(file.path(fig_dir, "Family_diverging.svg"), family_direction_plot,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")
