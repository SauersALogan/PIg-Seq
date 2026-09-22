# =============================================================================
# Script:       09_family_gsea.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Family-level gene set enrichment analysis using fgsea on
#               product-level t-ratios from family DEG results. NES heatmap
#               ordered by module correlation clustering and taxonomic
#               phylogeny.
# Dependencies: data_structuring.R, helpers.R, plotting.R, gene_sets.R
# Input:        family_product_df (from 08_family_deg_analysis.R or cached CSV)
#               functional_modules_filtered (from gene_sets.R)
#               Annotations (from gene_sets.R)
#               Taxa (from data_structuring.R)
# Output:       A .svg of family GSEA NES heatmap. Figure 2 panel H
# =============================================================================

# ---- Paths (update for your system) ----
fig_dir  <- "/path/to/figure_output_dir"
data_dir <- "/path/to/data"

library(fgsea)
library(ggplot2)
library(svglite)
library(tibble)
library(tidyr)

source("statistics/functions/plotting.R")
source("statistics/functions/helpers.R")
source("statistics/functions/data_structuring.R")
source("statistics/functions/gene_sets.R")

## ---- Load family DEG results ----
family_product_results <- read.csv(file.path(data_dir, "Family_DEGs_17-09-2026.csv"))

family_product_df <- family_product_results %>%
    dplyr::filter(contrast == "Positive - Negative") %>%
    group_by(Family) %>%
    mutate(p.adjust = p.adjust(p.value, method = "BH"),
        significant = p.adjust < 0.05,
        direction = case_when(significant & estimate > 0 ~ "Up in Positive",
            significant & estimate < 0 ~ "Down in Positive",
            TRUE ~ "Not significant"),
        neg_log10_p = -log10(p.adjust)) %>%
    ungroup()

## ---- Build product-level modules ----
functional_modules_products <- lapply(functional_modules_filtered, function(locus_tags) {
    Annotations %>%
        dplyr::filter(Locus_Tag %in% locus_tags) %>%
        dplyr::mutate(Product_clean = trimws(tolower(URLdecode(Product)))) %>%
        dplyr::pull(Product_clean) %>%
        unique()
})

min_genes <- 5
max_genes <- 5000
functional_modules_products <- functional_modules_products[
    sapply(functional_modules_products, length) >= min_genes &
        sapply(functional_modules_products, length) <= max_genes
]

## ---- Run GSEA per family ----
families <- unique(family_product_df$Family)

gsea_family <- lapply(setNames(families, families), function(fam) {
    ranked <- family_product_df %>%
        dplyr::filter(Family == fam, contrast == "Positive - Negative") %>%
        arrange(desc(t.ratio)) %>%
        dplyr::select(product, t.ratio) %>%
        deframe()

    if (length(ranked) < 10) return(NULL)

    fgsea(pathways = functional_modules_products,
        stats = ranked, minSize = 3, maxSize = 500)
}) %>%
    dplyr::bind_rows(.id = "Family")

cat("Families with results:", n_distinct(gsea_family$Family), "\n")
cat("Significant (padj < 0.05):",
    sum(gsea_family$padj < 0.05, na.rm = TRUE), "\n")
cat("Significant (padj < 0.25):",
    sum(gsea_family$padj < 0.25, na.rm = TRUE), "\n")

## ---- Heatmap data ----
heatmap_data <- gsea_family %>%
    dplyr::select(Family, pathway, NES, padj) %>%
    mutate(significant = padj < 0.25)

all_combos <- expand.grid(
    Family  = unique(heatmap_data$Family),
    pathway = unique(heatmap_data$pathway), stringsAsFactors = FALSE)

heatmap_data <- all_combos %>%
    left_join(heatmap_data, by = c("Family", "pathway"))

nes_matrix <- heatmap_data %>%
    dplyr::select(Family, pathway, NES) %>%
    tidyr::pivot_wider(names_from = pathway, values_from = NES, values_fill = 0) %>%
    tibble::column_to_rownames("Family") %>%
    as.matrix()

pathway_clust <- hclust(dist(t(nes_matrix)))
pathway_order <- pathway_clust$labels[pathway_clust$order]

family_order <- Taxa %>%
    tibble::rownames_to_column("Identifier") %>%
    dplyr::select(Phylum, Class, Order, Family) %>%
    dplyr::mutate(Family = gsub("f__", "", Family),
        Phylum = gsub("p__", "", Phylum),
        Class  = gsub("c__", "", Class),
        Order  = gsub("o__", "", Order)) %>%
    dplyr::filter(Family %in% unique(heatmap_data$Family)) %>%
    dplyr::distinct(Family, .keep_all = TRUE) %>%
    dplyr::arrange(Phylum, Class, Order, Family) %>%
    dplyr::pull(Family)

family_order <- union(family_order, unique(heatmap_data$Family))

heatmap_data <- heatmap_data %>%
    mutate(pathway = gsub("_", " ", pathway),
        Family  = factor(Family, levels = family_order),
        pathway = factor(pathway, levels = gsub("_", " ", pathway_order)),
        sig_label = case_when(padj < 0.05 ~ "*",
            padj < 0.25 ~ "*", TRUE ~ ""))

## ---- Plot ----
family_gsea_heatmap <- ggplot(heatmap_data,
        aes(x = pathway, y = Family, fill = NES)) +
    geom_tile(colour = "white", linewidth = lw_sm) +
    geom_text(aes(label = sig_label), size = 7 / .pt, family = "Arial",
        vjust = 0.75) +
    scale_fill_gradient2(low = col_neg, mid = "white", high = col_pos,
        midpoint = 0, na.value = "grey85", name = "NES") +
    labs(x = NULL, y = NULL) +
    theme_for_figures(show_labels = TRUE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.text = element_text(size = base_sz),
        legend.title = element_text(size = base_sz),
        legend.key.width = unit(2, "mm"),
        plot.margin = ggplot2::margin(1, 0, 0, 0, "mm"))
family_gsea_heatmap

## ---- Save ----
panel_w <- 84.6
panel_h <- 90.0

ggsave(file.path(fig_dir, "GSEA_heatmap_legend.svg"), family_gsea_heatmap,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")
