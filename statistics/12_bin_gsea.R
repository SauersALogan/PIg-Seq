# =============================================================================
# Script:       12_bin_gsea.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Bin-level gene set enrichment analysis with phylogenetic tree
#               and NES heatmap overlay. Module ordering by NES correlation
#               clustering, bin ordering by taxonomy. Includes bin metabolic
#               activity classification (active/average/inactive) via
#               residual regression.
# Dependencies: data_structuring.R, helpers.R, plotting.R, gene_sets.R
# Input:        Bin_DEGs_08-05-2026.csv (from 05_bin_deg_analysis.R)
#               bin.tree (https://doi.org/10.48420/33951214)
#               functional_modules_filtered (from gene_sets.R)
#               Data, rna_data (from data_structuring.R)
# Output:       A .svg of phylogenetic tree with NES heatmap. Figure 4 panel A
# =============================================================================

# ---- Paths (update for your system) ----
fig_dir  <- "/path/to/figure_output_dir"
data_dir <- "/path/to/data"

library(ape)
library(ggtree)
library(ggtreeExtra)
library(microbiome)
library(ggplot2)
library(cluster)
library(svglite)
library(patchwork)

source("statistics/functions/plotting.R")
source("statistics/functions/helpers.R")
source("statistics/functions/data_structuring.R")
source("statistics/functions/gene_sets.R")

## ---- Load DEGs ----
DEGs <- read.csv(file.path(data_dir, "Bin_DEGs_08-05-2026.csv"))
DEGs <- DEGs %>% dplyr::rename(clr_t = t.ratio)

## ---- Bin metabolic activity classification ----
Data_clr <- microbiome::transform(Data, "clr")
rna_clr <- microbiome::transform(rna_data, "clr")

dna_clr_mat <- as.data.frame(otu_table(Data_clr))
rna_clr_mat <- as.data.frame(otu_table(rna_clr))

dna_long_bins <- dna_clr_mat %>%
    rownames_to_column("Bin_id") %>%
    pivot_longer(-Bin_id, names_to = "Sample", values_to = "clr_dna")

rna_long_bins <- rna_clr_mat %>%
    rownames_to_column("Bin_id") %>%
    pivot_longer(-Bin_id, names_to = "Sample", values_to = "clr_rna")

bin_activity_data <- rna_long_bins %>%
    left_join(dna_long_bins, by = c("Bin_id", "Sample")) %>%
    left_join(sample_metadata, by = "Sample") %>%
    dplyr::filter(Fraction == "Native")

activity_model_regression <- lm(clr_rna ~ clr_dna,
    data = bin_activity_data)

bin_activity_data$residual <- residuals(activity_model_regression)

activity_regression <- bin_activity_data %>%
    group_by(Bin_id) %>%
    summarise(
        mean_residual = mean(residual),
        wilcox_p = wilcox.test(residual, mu = 0, exact = FALSE)$p.value,
        .groups = "drop") %>%
    mutate(
        p_adjust = p.adjust(wilcox_p, method = "fdr"),
        activity_status = case_when(
            p_adjust < 0.05 & mean_residual > 0 ~ "Active",
            p_adjust < 0.05 & mean_residual < 0 ~ "Inactive",
            TRUE ~ "Average"))

activity_regression_mapped <- activity_regression %>%
    left_join(Bin_map %>%
        mutate(Bin = paste0("Bin.", Bin)) %>%
        dplyr::select(Bin, Taxonomy, Identifier),
        by = c("Bin_id" = "Identifier"))

prop_deg_per_bin <- DEGs %>%
    dplyr::filter(contrast == "Positive - Negative") %>%
    group_by(Bin) %>%
    summarise(
        n_degs = sum(significant, na.rm = TRUE),
        n_total = n(),
        prop_deg = n_degs / n_total,
        .groups = "drop")

deg_by_activity <- prop_deg_per_bin %>%
    dplyr::left_join(
        activity_regression_mapped %>%
            dplyr::select(Bin, activity_status, Taxonomy),
        by = "Bin") %>%
    dplyr::mutate(pct_deg = prop_deg * 100)

## ---- Select bins for GSEA ----
gsea_bins <- deg_by_activity %>%
    dplyr::filter(pct_deg > 1) %>%
    pull(Bin)

## ---- Run GSEA ----
gsea_results <- run_fgsea_parallel(
    bin_list        = gsea_bins,
    DEGs            = DEGs,
    selected_contrast = "Positive - Negative",
    gene_set_data   = functional_modules_filtered,
    analysis_name   = "Bin-level GSEA",
    min_genes       = 100,
    min_set_size    = 5,
    max_set_size    = 500,
    n_cores         = parallel::detectCores() - 1)

## ---- Module correlation clustering ----
nes_matrix <- gsea_results %>%
    tidyr::pivot_wider(id_cols = pathway, names_from = Bin,
        values_from = NES) %>%
    tibble::column_to_rownames("pathway") %>%
    as.matrix()

nes_cor <- cor(t(nes_matrix), use = "pairwise.complete.obs")

## ---- Load and filter tree ----
tree <- read.tree(file.path(data_dir, "bin.tree"))
bin_patterns <- c("^S3_", "^S6_", "^S13_", "^S20_", "^S22_",
    "^S27_", "^S29_", "^S453_")
tree_bin_names <- tree$tip.label[grepl(paste(bin_patterns,
    collapse = "|"), tree$tip.label)]
tree_bins <- keep.tip(tree, tree_bin_names)
tree_bins_ladderized <- ladderize(tree_bins)

## ---- Supporting data ----
tip_bin_map <- Bin_map %>%
    dplyr::rename(tree_label = Identifier) %>%
    dplyr::select(tree_label, Bin) %>%
    dplyr::mutate(Bin = paste0("Bin.", Bin))

bin_species <- Bin_map %>%
    dplyr::rename(tree_label = Identifier) %>%
    dplyr::select(tree_label, Bin) %>%
    dplyr::mutate(Bin = paste0("Bin.", Bin)) %>%
    dplyr::left_join(
        Taxa %>% tibble::rownames_to_column("tree_label") %>%
            dplyr::select(tree_label, Species),
        by = "tree_label") %>%
    dplyr::mutate(
        Species_clean = gsub("_", " ", Species),
        Species_label = paste0(Species_clean, " (", Bin, ")"))

chlamydia_bin <- deg_by_activity %>%
    dplyr::filter(grepl("Chlamydia", Taxonomy)) %>%
    pull(Bin)

gsea_bins_plus_outgroup <- c(gsea_bins, chlamydia_bin)

tree_tips_to_keep <- tip_bin_map %>%
    dplyr::filter(Bin %in% gsea_bins_plus_outgroup) %>%
    pull(tree_label)

tree_bins_focused <- keep.tip(tree_bins, tree_tips_to_keep)
tree_bins_focused <- ladderize(tree_bins_focused)

## ---- Heatmap data ----
tip_family <- data.frame(tree_label = tree_bins_focused$tip.label) %>%
    dplyr::left_join(
        Taxa %>% tibble::rownames_to_column("tree_label") %>%
            dplyr::select(tree_label, Family) %>%
            dplyr::mutate(Family = gsub("f__", "", Family)),
        by = "tree_label")

custom_tested_count <- gsea_results %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::group_by(pathway) %>%
    dplyr::summarise(n_bins_tested = n_distinct(Bin), .groups = "drop")

modules_to_keep <- custom_tested_count %>%
    dplyr::filter(n_bins_tested > 2) %>%
    dplyr::pull(pathway)

nes_cor_filtered <- nes_cor[
    rownames(nes_cor) %in% modules_to_keep,
    colnames(nes_cor) %in% modules_to_keep]
nes_cor_filtered[is.na(nes_cor_filtered)] <- 0
module_dist <- as.dist(1 - nes_cor_filtered)
module_hc <- hclust(module_dist, method = "average")
module_order <- rownames(nes_cor_filtered)[module_hc$order]

nes_for_tree <- gsea_results %>%
    dplyr::filter(pathway %in% modules_to_keep) %>%
    dplyr::left_join(bin_species %>% dplyr::select(Bin, tree_label),
        by = "Bin") %>%
    dplyr::filter(!is.na(tree_label)) %>%
    dplyr::filter(tree_label %in% tree_tips_to_keep) %>%
    dplyr::select(tree_label, pathway, NES) %>%
    tidyr::pivot_wider(names_from = pathway, values_from = NES) %>%
    tibble::column_to_rownames("tree_label") %>%
    as.data.frame()

nes_for_tree <- nes_for_tree[, module_order[module_order %in%
    colnames(nes_for_tree)]]
nes_for_tree_clean <- nes_for_tree
colnames(nes_for_tree_clean) <- gsub("_", " ", colnames(nes_for_tree_clean))

## ---- Order color mapping ----
tip_order <- data.frame(tree_label = tree_bins_focused$tip.label) %>%
    dplyr::left_join(
        Taxa %>% tibble::rownames_to_column("tree_label") %>%
            dplyr::select(tree_label, Order) %>%
            dplyr::mutate(Order = gsub("o__", "", Order)),
        by = "tree_label") %>%
    dplyr::mutate(Order_group = ifelse(
        Order %in% c("Bacteroidales", "Lachnospirales", "Erysipelotrichales",
            "Lactobacillales", "Oscillospirales"),
        Order, "Other"))

order_colors <- c(
    "Bacteroidales"      = "#6644AA",
    "Lachnospirales"     = "#228833",
    "Erysipelotrichales" = "#EE6677",
    "Lactobacillales"    = "#EE7733",
    "Oscillospirales"    = "#66CCEE",
    "Other"              = "#BBBBBB")

## ---- Tree (grey branches, tip dots for order) ----
p_tree_focused <- ggtree(tree_bins_focused, color = "grey80",
        linewidth = lw_sm * 2) %<+%
    (bin_species %>% dplyr::select(tree_label, Species_label) %>%
        dplyr::left_join(tip_order, by = "tree_label")) +
    geom_tiplab(aes(label = gsub(" \\(Bin\\.\\d+\\)", "", Species_label)),
        size = base_sz / .pt, family = "Arial", align = TRUE, color = "black") +
    theme_tree() +
    theme(text = element_text(family = "Arial"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

## ---- Combined tree + heatmap ----
p_combined <- gheatmap(p_tree_focused, nes_for_tree_clean,
        offset = 3.6, width = 2,
        colnames = FALSE,
        color = "black") +
    scale_fill_gradient2(low = col_neg, mid = "white", high = col_pos,
        midpoint = 0, limits = c(-2, 2), oob = scales::squish,
        name = "NES", na.value = "#E8E8E8") +
    geom_tippoint(aes(color = Order_group), size = pt_mean * 0.5) +
    scale_color_manual(name = "Order", values = order_colors, drop = TRUE) +
    coord_cartesian(clip = "off") +
    theme(text = element_text(family = "Arial"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none")
p_combined

## ---- Save ----
panel_w <- 142
panel_h <- 84

ggsave(file.path(fig_dir, "tree_heatmap_combined.svg"), p_combined,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")
