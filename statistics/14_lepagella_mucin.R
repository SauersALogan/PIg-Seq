# =============================================================================
# Script:       14_lepagella_mucin.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Mucin-degrading enzyme stacked bar plot for Lepagella muris
#               MAGs showing proportion of genes significantly up, down,
#               or non-significant (with directional log2 FC > 1.5 threshold)
#               in Positive vs Negative fractions. 
# Dependencies: plotting.R
# Input:        mucin_master_table.csv (https://doi.org/10.48420/33951214)
#               Bin_DEGs_08-05-2026.csv (from 05_bin_deg_analysis.R)
# Output:       A .svg of mucin enzyme stacked bar. Figure 4 Panel B
# =============================================================================

# ---- Paths (update for your system) ----
fig_dir  <- "/path/to/figure_output_dir"
data_dir <- "/path/to/data"

library(ggplot2)
library(svglite)
library(readr)
library(dplyr)
library(tidyr)

source("statistics/functions/plotting.R")

## ---- Load data ----
mucin_master_table <- read_csv(file.path(data_dir, "mucin_master_table.csv"))
DEGs <- read.csv(file.path(data_dir, "Bin_DEGs_08-05-2026.csv"))

## ---- Clean category and assign substrate ----
mucin_master_table <- mucin_master_table %>%
    mutate(
        Category = trimws(Category),
        Substrate = case_when(
            Category %in% c("Sialate-OAE", "GH110", "GH88", "GH88?", "GH20",
                "GH33", "Sulfatase", "Peptide M64",
                "Exo-alpha-sialidase (CAZy unconfirmed)") ~ "mammalian",
            Category %in% c("GH36", "GH127", "GH146", "GH43", "GH16",
                "PL42") ~ "plant",
            Category %in% c("GH97", "GH27", "GH2", "GH92", "GH3",
                "GH35") ~ "broad",
            Category %in% c("SusC", "SusD", "SusE/F") ~ NA_character_,
            TRUE ~ NA_character_))

## ---- Merge with DEG stats ----
mucin_master_table <- mucin_master_table %>%
    left_join(
        DEGs %>%
            filter(contrast == "Positive - Negative") %>%
            select(Geneid, estimate, p_adjust, significant),
        by = c("Locus_Tag" = "Geneid"))

## ---- Mucin categories ----
enzyme_from_master <- mucin_master_table %>%
    filter(!Category %in% c("SusC", "SusD", "SusE/F")) %>%
    mutate(enzyme_category = Category) %>%
    select(Bin, enzyme_category, estimate, p_adjust, significant)

sus_from_master <- mucin_master_table %>%
    filter(Category %in% c("SusC", "SusD", "SusE/F")) %>%
    mutate(enzyme_category = Category) %>%
    select(Bin, enzyme_category, estimate, p_adjust, significant)

## ---- Combine and exclude plant-substrate categories ----
all_enzyme_degs <- bind_rows(enzyme_from_master, sus_from_master)

all_enzyme_degs_host_relevant <- all_enzyme_degs %>%
    left_join(
        mucin_master_table %>% select(Category, Substrate) %>% distinct(),
        by = c("enzyme_category" = "Category")) %>%
    filter(is.na(Substrate) | Substrate != "plant")

## ---- Enzyme stacked bars ----
fc_cutoff <- log2(1.5)

fill_colors <- c(
    "Significant down" = col_neg,
    "Down (NS)"        = "#A3CEE2",
    "Neutral"          = "grey85",
    "Up (NS)"          = "#EAABB3",
    "Significant up"   = col_pos)

category_order <- unique(all_enzyme_degs_host_relevant$enzyme_category) %>%
    setdiff(c("SusC", "SusD", "SusE/F"))
category_order <- c("SusC", "SusD", "SusE/F", sort(category_order))

enzyme_summary <- all_enzyme_degs_host_relevant %>%
    mutate(status = case_when(
        significant & estimate > 0              ~ "Sig_Up",
        significant & estimate < 0              ~ "Sig_Down",
        !significant & estimate > fc_cutoff     ~ "NS_Up",
        !significant & estimate < -fc_cutoff    ~ "NS_Down",
        TRUE                                    ~ "Neutral")) %>%
    group_by(Bin, enzyme_category) %>%
    summarise(
        n_genes = n(),
        n_sig_up   = sum(status == "Sig_Up"),
        n_ns_up    = sum(status == "NS_Up"),
        n_neutral  = sum(status == "Neutral"),
        n_ns_down  = sum(status == "NS_Down"),
        n_sig_down = sum(status == "Sig_Down"),
        .groups = "drop") %>%
    complete(Bin, enzyme_category,
        fill = list(n_genes = 0, n_sig_up = 0, n_ns_up = 0,
            n_neutral = 0, n_ns_down = 0, n_sig_down = 0))

enzyme_plot_data <- enzyme_summary %>%
    mutate(enzyme_category = factor(enzyme_category,
        levels = rev(category_order))) %>%
    pivot_longer(cols = c(n_sig_down, n_ns_down, n_neutral, n_ns_up, n_sig_up),
        names_to = "segment", values_to = "n") %>%
    group_by(Bin, enzyme_category) %>%
    mutate(pct = if (n_genes[1] > 0) n / n_genes[1] * 100 else 0) %>%
    ungroup() %>%
    mutate(segment = factor(segment,
        levels = c("n_sig_down", "n_ns_down", "n_neutral", "n_ns_up", "n_sig_up"),
        labels = c("Significant down", "Down (NS)", "Neutral",
            "Up (NS)", "Significant up")))

label_data <- enzyme_summary %>%
    mutate(enzyme_category = factor(enzyme_category,
        levels = rev(category_order))) %>%
    select(Bin, enzyme_category, n_genes)

enzyme_short_names <- c(
    "SusC"            = "SusC",
    "SusD"            = "SusD",
    "SusE/F"          = "SusE/F",
    "Sialate-OAE"     = "Sialate OAE",
    "GH110"           = "GH110",
    "GH88"            = "GH88",
    "GH88?"           = "GH88?",
    "GH20"            = "GH20",
    "GH33"            = "GH33",
    "Sulfatase"       = "Sulfatase",
    "Peptide M64"     = "M64",
    "Exo-alpha-sialidase (CAZy unconfirmed)" = "Sialidase",
    "GH97"            = "GH97",
    "GH27"            = "GH27",
    "GH2"             = "GH2",
    "GH92"            = "GH92",
    "GH3"             = "GH3",
    "GH35"            = "GH35")

## ---- Plot ----
enzyme_stacked_plot <- ggplot(enzyme_plot_data,
        aes(x = enzyme_category, y = pct, fill = segment)) +
    geom_col(width = 0.7) +
    geom_text(data = label_data,
        aes(x = enzyme_category, y = 103, label = paste0("n=", n_genes)),
        inherit.aes = FALSE, size = base_sz / .pt * 0.8, hjust = 0,
        colour = "black", family = "Arial") +
    facet_wrap(~ Bin, ncol = 2) +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = enzyme_short_names) +
    scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 25),
        labels = seq(0, 100, 25)) +
    scale_fill_manual(name = NULL, values = fill_colors) +
    labs(x = NULL, y = "Percent of annotated genes") +
    theme_for_figures(show_labels = TRUE) +
    theme(legend.position = "none",
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        plot.margin = ggplot2::margin(1, 1, 1, 1, "mm"))
enzyme_stacked_plot

## ---- Save ----
panel_w <- 70.1
panel_h <- 45

ggsave(file.path(fig_dir, "enzyme_stacked.svg"), enzyme_stacked_plot,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")
