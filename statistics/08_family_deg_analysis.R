# =============================================================================
# Script:       08_family_deg_analysis.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Family-level differential product expression analysis using
#               parallelized lmer models with Bin as random effect. Product
#               names deduplicated via sorted-token normalization before
#               modeling. Produces per-family volcano plots with highlighted
#               genes.
# Dependencies: data_structuring.R, helpers.R, plotting.R
# Input:        pseudobulk_filtered, Taxa, Bin_map (from data_structuring.R)
#               Or cached: Family_DEGs_17-09-2026.csv
# Output:       Family_DEGs_17-09-2026.csv (model results)
#               A .svg of Muribaculaceae volcano. Figure 2 panel D
#               A .svg of Oscillospiraceae volcano. Figure 2 panel E
#               A .svg of Bacteroidaceae volcano. Figure 2 panel F
#               A .svg of Erysipelotrichaceae volcano. Figure 2 panel G
# =============================================================================

# ---- Paths (update for your system) ----
fig_dir  <- "/path/to/figure_output_dir"
data_dir <- "/path/to/data"

library(lme4)
library(emmeans)
library(MuMIn)
library(future)
library(future.apply)
library(progressr)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(svglite)
library(purrr)

options(contrasts = c("contr.sum", "contr.poly"))
filter <- dplyr::filter

n_cores <- parallel::detectCores() - 1
options(future.globals.maxSize = 2000 * 1024^2)
plan(multisession, workers = n_cores)
handlers(global = TRUE)
handlers("progress")

source("statistics/functions/plotting.R")
source("statistics/functions/helpers.R")
source("statistics/functions/data_structuring.R")

## ---- Append family data ----
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

families_2plus <- pseudobulk_with_family %>%
    group_by(Family) %>%
    summarise(n_bins = n_distinct(Bin)) %>%
    dplyr::filter(n_bins >= 2) %>%
    pull(Family)

post_clr_family <- pseudobulk_with_family %>%
    dplyr::filter(!is.na(Family), Family %in% families_2plus,
        Fraction %in% c("Negative", "Positive")) %>%
    mutate(Fraction = factor(Fraction, levels = c("Positive", "Negative")))

## ---- Filter to products in 2+ bins per family ----
product_bin_counts <- post_clr_family %>%
    dplyr::filter(!is.na(clr_rna), !is.na(clr_dna)) %>%
    dplyr::group_by(Family, product) %>%
    dplyr::summarise(n_bins = dplyr::n_distinct(Bin), .groups = "drop")

valid_family_products <- product_bin_counts %>%
    dplyr::filter(n_bins >= 2)

family_clr_filtered <- post_clr_family %>%
    dplyr::semi_join(valid_family_products, by = c("Family", "product"))

## ---- Product name deduplication (sorted-token normalization) ----
normalised <- family_clr_filtered %>%
    dplyr::distinct(product) %>%
    dplyr::mutate(norm = sapply(product, function(p) {
        paste(sort(tolower(strsplit(trimws(p), "\\s+")[[1]])), collapse = " ")
    }))

canonical_map <- normalised %>%
    dplyr::group_by(norm) %>%
    dplyr::summarise(canonical = sort(product)[1],
        all_forms = list(product), .groups = "drop") %>%
    tidyr::unnest(all_forms) %>%
    dplyr::select(product = all_forms, canonical)

family_clr_deduped <- family_clr_filtered %>%
    dplyr::left_join(canonical_map, by = "product") %>%
    dplyr::mutate(product = dplyr::coalesce(canonical, product)) %>%
    dplyr::select(-canonical)

## ---- Filter to 4+ samples in 2+ fractions (post-dedup) ----
valid_products_deduped <- family_clr_deduped %>%
    group_by(Family, product, Fraction) %>%
    summarise(n_samples = n_distinct(Mouse), .groups = "drop") %>%
    dplyr::filter(n_samples >= 4) %>%
    group_by(Family, product) %>%
    summarise(n_fractions = n_distinct(Fraction), .groups = "drop") %>%
    dplyr::filter(n_fractions >= 2)

family_clr_final <- family_clr_deduped %>%
    semi_join(valid_products_deduped, by = c("Family", "product"))

final_bin_check <- family_clr_final %>%
    dplyr::filter(!is.na(clr_rna), !is.na(clr_dna)) %>%
    dplyr::group_by(Family, product) %>%
    dplyr::summarise(n_bins_complete = dplyr::n_distinct(Bin), .groups = "drop") %>%
    dplyr::filter(n_bins_complete < 2)

cat("\nRemaining family-product combos with <2 complete-case bins:",
    nrow(final_bin_check), "\n")

family_product_list <- split(family_clr_final, list(family_clr_final$Family,
    family_clr_final$product), drop = TRUE)
cat("Total family-product combinations:", length(family_product_list), "\n")

## ---- Run or load family models ----
family_degs_file <- file.path(data_dir, "Family_DEGs_17-09-2026.csv")
if (file.exists(family_degs_file)) {
    cat("Loading existing family DEGs file from:", family_degs_file, "\n")
    family_product_results <- read.csv(family_degs_file)
} else {
    cat("No existing family DEGs file, running models \n")
    family_product_results <- run_lmer_parallel(data_list = family_product_list,
        formula = clr_rna ~ clr_dna + Fraction + (1|Mouse) + (1|Bin),
        emmeans_term = ~ Fraction, fdr_group = c("contrast", "Family"),
        singular_action = "flag_keep",
        output_cols = list(Family = quote(unique(Family)[1]),
            product = quote(unique(product)[1]),
            n_bins = quote(dplyr::n_distinct(Bin)),
            n_mice = quote(dplyr::n_distinct(Mouse)),
            n_obs = quote(nrow(d))))
    write.csv(family_product_results, family_degs_file, row.names = FALSE)
}

## ---- Structure for plotting ----
family_product_df <- family_product_results %>%
    dplyr::filter(contrast == "Positive - Negative") %>%
    group_by(Family) %>%
    mutate(p_adjust = p.adjust(p.value, method = "BH"),
        significant = p_adjust < 0.05,
        direction = case_when(significant & estimate > 0 ~ "Up in Positive",
            significant & estimate < 0 ~ "Down in Positive",
            TRUE ~ "Not significant"),
        neg_log10_p = -log10(p_adjust)) %>%
    ungroup()

## ---- Muribaculaceae volcano ----
muri_highlight <- family_product_df %>%
    dplyr::filter(Family == "Muribaculaceae", significant == TRUE,
        (product == "thioredoxin" |
         grepl("\\bgld[elmn]\\b|type ix secretion system membrane",
            tolower(product)))) %>%
    dplyr::arrange(p_adjust) %>%
    dplyr::mutate(label = dplyr::case_when(
        product == "thioredoxin"        ~ "Thioredoxin",
        grepl("glde", tolower(product)) ~ "GldE",
        grepl("gldm", tolower(product)) ~ "GldM",
        grepl("gldl", tolower(product)) ~ "GldL",
        grepl("gldn", tolower(product)) ~ "GldN",
        grepl("membrane", tolower(product)) ~ "T9SS",
        TRUE ~ NA_character_))

muri_highlight <- dplyr::bind_rows(
    muri_highlight %>% dplyr::filter(product == "thioredoxin"),
    muri_highlight %>% dplyr::filter(product != "thioredoxin") %>%
        dplyr::slice_head(n = 3))

Muri_volcano <- plot_family_product_volcano(family_product_df, "Muribaculaceae") +
    geom_point(data = muri_highlight,
        aes(x = estimate, y = neg_log10_p),
        shape = 21, fill = "#117733", colour = "black",
        stroke = lw_sm, size = pt_mean) +
    theme_for_figures(show_labels = TRUE) +
    ggrepel::geom_text_repel(data = muri_highlight,
        aes(x = estimate, y = neg_log10_p, label = label),
        size = base_sz / .pt, family = "Arial",
        min.segment.length = 0, segment.size = lw_sm, segment.color = "grey50",
        box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf)

## ---- Oscillospiraceae volcano ----
oscil_highlight <- family_product_df %>%
    dplyr::filter(Family == "Oscillospiraceae", significant == TRUE,
        product %in% c(
            "tail sheath protein subtilisin-like domain-containing protein",
            "holin",
            "spoivb peptidase",
            "stage v sporulation protein ac")) %>%
    dplyr::mutate(label = dplyr::case_when(
        grepl("tail sheath protein subtilisin-like domain-containing protein",
            tolower(product)) ~ "Tail",
        grepl("holin", tolower(product)) ~ "Holin",
        grepl("spoivb peptidase", tolower(product)) ~ "SpoIVB",
        grepl("stage v sporulation protein ac", tolower(product)) ~ "SpoVAC",
        TRUE ~ NA_character_))

Oscil_volcano <- plot_family_product_volcano(family_product_df, "Oscillospiraceae") +
    geom_point(data = oscil_highlight,
        aes(x = estimate, y = neg_log10_p),
        shape = 21, fill = "#117733", colour = "black",
        stroke = lw_sm, size = pt_mean) +
    theme_for_figures(show_labels = TRUE) +
    ggrepel::geom_text_repel(data = oscil_highlight,
        aes(x = estimate, y = neg_log10_p, label = label),
        size = base_sz / .pt, family = "Arial",
        min.segment.length = 0, segment.size = lw_sm, segment.color = "grey50",
        box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf)

## ---- Bacteroidaceae volcano ----
bacto_pul_highlight <- family_product_df %>%
    dplyr::filter(Family == "Bacteroidaceae", significant == TRUE,
        product %in% c("glycoside hydrolase family 2 protein",
                       "glycoside hydrolase 97",
                       "glycoside hydrolase family 88 protein")) %>%
    dplyr::mutate(label = dplyr::case_when(
        grepl("family 2",  product) ~ "GH2",
        grepl("97",        product) ~ "GH97",
        grepl("family 88", product) ~ "GH88"))

Bacto_volcano <- plot_family_product_volcano(family_product_df, "Bacteroidaceae") +
    geom_point(data = bacto_pul_highlight,
        aes(x = estimate, y = neg_log10_p),
        shape = 21, fill = "#117733", colour = "black",
        stroke = lw_sm, size = pt_mean) +
    theme_for_figures(show_labels = TRUE) +
    ggrepel::geom_text_repel(data = bacto_pul_highlight,
        aes(x = estimate, y = neg_log10_p, label = label),
        size = base_sz / .pt, family = "Arial",
        min.segment.length = 0, segment.size = lw_sm, segment.color = "grey50",
        box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf)

## ---- Erysipelotrichaceae volcano ----
erysi_highlight <- family_product_df %>%
    dplyr::filter(Family == "Erysipelotrichaceae", significant == TRUE,
        product %in% c("sortase", "hemolysin", "adhesin")) %>%
    dplyr::mutate(label = dplyr::case_when(
        grepl("sortase",  product) ~ "Sortase",
        grepl("hemolysin", product) ~ "Hemolysin",
        grepl("adhesin", product) ~ "Adhesin"))

Erysi_volcano <- plot_family_product_volcano(family_product_df, "Erysipelotrichaceae") +
    geom_point(data = erysi_highlight,
        aes(x = estimate, y = neg_log10_p),
        shape = 21, fill = "#117733", colour = "black",
        stroke = lw_sm, size = pt_mean) +
    theme_for_figures(show_labels = TRUE) +
    ggrepel::geom_text_repel(data = erysi_highlight,
        aes(x = estimate, y = neg_log10_p, label = label),
        size = base_sz / .pt, family = "Arial",
        min.segment.length = 0, segment.size = lw_sm, segment.color = "grey50",
        box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf)

## ---- Save ----
panel_w <- 40
panel_h <- 43

ggsave(file.path(fig_dir, "Muri_volcano.svg"), Muri_volcano,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")

ggsave(file.path(fig_dir, "Oscil_volcano.svg"), Oscil_volcano,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")

ggsave(file.path(fig_dir, "Bacto_volcano.svg"), Bacto_volcano,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")

ggsave(file.path(fig_dir, "Erysi_volcano.svg"), Erysi_volcano,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")
