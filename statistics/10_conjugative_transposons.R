# =============================================================================
# Script:       10_conjugative_transposons.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Conjugative transposon and prophage gene expression plot
#               for Bacteroidaceae and Muribaculaceae. Model estimates with
#               corresponding per-mouse, per-bin positive-minus-negative CLR
#               contrasts shown as background points.
# Dependencies: data_structuring.R, plotting.R
# Input:        Family_DEGs_17-09-2026.csv (from 08_family_deg_analysis.R)
#               family_clr_final (from 08_family_deg_analysis.R)
# Output:       A .svg of conjugation/prophage expression. Figure 3 panel A
# =============================================================================

# ---- Paths (update for your system) ----
fig_dir  <- "/path/to/figure_output_dir"
data_dir <- "/path/to/data"

library(ggplot2)
library(svglite)
library(stringr)

source("statistics/functions/plotting.R")
source("statistics/functions/data_structuring.R")

## ---- Load family DEG results ----
family_product_results <- read.csv(file.path(data_dir, "Family_DEGs_17-09-2026.csv"))

## ---- Tra genes ----
transposon_products <- family_product_results %>%
    dplyr::filter(
        Family %in% c("Bacteroidaceae", "Muribaculaceae"),
        grepl("^conjugative transposon protein tra[a-z]$", product,
            ignore.case = TRUE)) %>%
    dplyr::bind_rows(
        family_product_results %>%
            dplyr::filter(
                Family %in% c("Bacteroidaceae", "Muribaculaceae"),
                product %in% c("trag family conjugative transposon atpase",
                    "tral conjugative transposon family protein"))
    ) %>%
    dplyr::distinct(Family, product, .keep_all = TRUE) %>%
    dplyr::mutate(
        short_name = paste0("Tra", toupper(stringr::str_sub(
            stringr::str_extract(product, "\\btra[a-z]\\b"), -1))),
        Family = as.character(Family))

## ---- Conjugation accessory (deduplicate per Family x short_name) ----
conj_acc <- family_product_results %>%
    dplyr::filter(
        Family %in% c("Bacteroidaceae", "Muribaculaceae"),
        product %in% c("mobilization protein", "relaxase",
            "mobv family relaxase")) %>%
    dplyr::mutate(
        short_name = case_when(
            grepl("mobilization", product) ~ "Mob",
            grepl("relaxase", product) ~ "Relaxase"),
        Family = as.character(Family)) %>%
    dplyr::group_by(Family, short_name) %>%
    dplyr::slice_min(p_adjust, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

## ---- Prophage genes ----
prophage <- family_product_results %>%
    dplyr::filter(
        Family %in% c("Bacteroidaceae", "Muribaculaceae"),
        product %in% c("phage tail protein",
            "phage major capsid protein",
            "phage tail tape measure protein",
            "phage n-6-adenine-methyltransferase")) %>%
    dplyr::distinct(Family, product, .keep_all = TRUE) %>%
    dplyr::mutate(
        short_name = case_when(
            product == "phage tail protein" ~ "Tail",
            product == "phage major capsid protein" ~ "Capsid",
            product == "phage tail tape measure protein" ~ "TMP",
            product == "phage n-6-adenine-methyltransferase" ~ "MTase"),
        Family = as.character(Family))

## ---- Combine model results ----
all_products <- dplyr::bind_rows(transposon_products, conj_acc, prophage) %>%
    dplyr::mutate(Family = factor(Family,
        levels = c("Bacteroidaceae", "Muribaculaceae")))

gene_order <- c(sort(unique(transposon_products$short_name)),
    "Mob", "Relaxase",
    "Tail", "Capsid", "TMP", "MTase")

## ---- Raw observations: tra genes ----
# NOTE: family_clr_final must be reconstructed here or loaded
# from 08_family_deg_analysis.R environment

## ---- Reconstruct family_clr_final for raw observations ----
bin_family_map <- Taxa %>%
    as.data.frame() %>%
    rownames_to_column("Identifier") %>%
    left_join(Bin_map %>%
        dplyr::rename(Identifier = Identifier), by = "Identifier") %>%
    mutate(Bin = paste0("Bin.", Bin)) %>%
    dplyr::select(Bin, Family) %>%
    mutate(Family = gsub("f__", "", Family)) %>%
    distinct()

pseudobulk_with_family <- pseudobulk_filtered %>%
    left_join(bin_family_map, by = "Bin")

post_clr_family <- pseudobulk_with_family %>%
    dplyr::filter(!is.na(Family),
        Family %in% c("Bacteroidaceae", "Muribaculaceae"),
        Fraction %in% c("Negative", "Positive")) %>%
    mutate(Fraction = factor(Fraction, levels = c("Positive", "Negative")))

tra_raw <- post_clr_family %>%
    dplyr::filter(product %in% transposon_products$product) %>%
    dplyr::mutate(
        short_name = paste0("Tra", toupper(stringr::str_sub(
            stringr::str_extract(product, "\\btra[a-z]\\b"), -1))))

## ---- Raw observations: accessory + prophage ----
other_products <- c("mobilization protein", "phage tail protein",
    "phage major capsid protein", "phage tail tape measure protein",
    "phage n-6-adenine-methyltransferase")

other_raw <- post_clr_family %>%
    dplyr::filter(product %in% other_products |
        grepl("relaxase", product, ignore.case = TRUE)) %>%
    dplyr::mutate(
        short_name = case_when(
            grepl("mobilization", product) ~ "Mob",
            grepl("relaxase", product) ~ "Relaxase",
            product == "phage tail protein" ~ "Tail",
            product == "phage major capsid protein" ~ "Capsid",
            product == "phage tail tape measure protein" ~ "TMP",
            product == "phage n-6-adenine-methyltransferase" ~ "MTase"))

## ---- Raw observations matching plotted model products exactly ----
selected_products <- all_products %>%
    dplyr::mutate(Family = as.character(Family)) %>%
    dplyr::select(Family, product, short_name) %>%
    dplyr::distinct()

all_raw <- post_clr_family %>%
    dplyr::mutate(Family = as.character(Family)) %>%
    dplyr::inner_join(
        selected_products,
        by = c("Family", "product")
    )

## ---- Compute per-bin contrasts ----
all_contrasts <- all_raw %>%
    dplyr::group_by(Family, short_name, Mouse, Bin, Fraction) %>%
    dplyr::summarise(
        clr_rna = mean(clr_rna, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
        id_cols = c(Family, short_name, Mouse, Bin),
        names_from = Fraction,
        values_from = clr_rna
    ) %>%
    dplyr::filter(!is.na(Positive), !is.na(Negative)) %>%
    dplyr::mutate(
        raw_contrast = Positive - Negative,
        Family = factor(
            Family,
            levels = c("Bacteroidaceae", "Muribaculaceae")
        )
    )

## ---- X-axis break position ----
n_conj <- length(unique(transposon_products$short_name)) + 2
break_1 <- n_conj + 0.5

## ---- Plot ----
transposon_plot <- ggplot() +
    geom_vline(xintercept = break_1, colour = "grey70",
        linewidth = lw_sm, linetype = "dashed") +
    geom_jitter(data = all_contrasts,
        aes(x = short_name, y = raw_contrast, color = Family, group = Family),
        position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
        size = pt_small * 0.25, alpha = 0.3, na.rm = TRUE) +
    geom_point(data = all_products,
        aes(x = short_name, y = estimate, color = Family, group = Family),
        position = position_dodge(width = 0.6), size = pt_mean * 0.35,
        na.rm = TRUE) +
    geom_errorbar(data = all_products,
        aes(x = short_name, ymin = estimate - SE, ymax = estimate + SE,
            color = Family, group = Family),
        position = position_dodge(width = 0.6),
        width = 0.15, linewidth = lw_sm,
        na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed",
        color = "grey50", linewidth = lw_sm) +
    scale_color_manual(name = "Family",
        values = c("Bacteroidaceae"  = "#332288",
                   "Muribaculaceae"  = "#882255")) +
    scale_x_discrete(limits = gene_order) +
    scale_y_continuous(limits = c(NA, NA)) +
    theme_for_figures(show_labels = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = base_sz),
        legend.title = element_text(size = base_sz),
        legend.box.spacing = unit(1, "mm"),
        plot.margin = ggplot2::margin(1, 0, 0, 5, "mm")) +
    labs(x = NULL, y = "Estimate (CLR)", color = "Family")
transposon_plot

## ---- Save ----
panel_w <- 101.8
panel_h <- 19 * 2

ggsave(file.path(fig_dir, "transposon_plot.svg"), transposon_plot,
    device = svglite::svglite, width = panel_w,
    height = panel_h, units = "mm")
