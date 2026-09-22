# =============================================================================
# Script:       11_phage_analysis.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Phage differential gene expression analysis using lmer with
#               host bin-level CLR DNA abundance as covariate. FDR correction
#               per contrast across all phage genes. Produces per-provirus
#               genome maps showing estimate direction and functional regions.
# Dependencies: data_structuring.R, helpers.R, plotting.R
# Input:        DNA_phage_features.csv (https://doi.org/10.48420/33951214)
#               RNA_phage_features.csv (https://doi.org/10.48420/33951214)
#               phynteny_featurecounts_mapping.csv (https://doi.org/10.48420/33951214)
#               pharokka_cds_final_merged_output.tsv (https://doi.org/10.48420/33951214)
#               phage_combined_annotations_filtered.tsv (https://doi.org/10.48420/33951214)
#               Data (phyloseq, from data_structuring.R)
# Output:       A .svg of phage genome maps. Figure 4 panel B
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
library(ggnewscale)
library(gggenes)
library(patchwork)
library(svglite)
library(microbiome)

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

## ---- Load phage data ----
phage_gene_counts <- read.csv(file.path(data_dir, "DNA_phage_features.csv"),
    row.names = 1)
phage_gene_expression <- read.csv(file.path(data_dir, "RNA_phage_features.csv"),
    row.names = 1)
phynteny_map <- read.csv(file.path(data_dir, "phynteny_featurecounts_mapping.csv")) %>%
    dplyr::rename(viral_function = gene_function)
pharokka_cds <- read.delim(file.path(data_dir, "pharokka_cds_final_merged_output.tsv"),
    sep = "\t")
phage_annotations_combined <- read.table(
    file.path(data_dir, "phage_combined_annotations_filtered.tsv"),
    sep = "\t", header = TRUE)

## ---- Bin-level CLR abundance (covariate) ----
Data_clr <- microbiome::transform(Data, "clr")
bin_clr_abundance <- as.data.frame(otu_table(Data_clr))

## ---- Structure phage data ----
phage_expression_long <- phage_gene_expression %>%
    rownames_to_column("Geneid") %>%
    pivot_longer(cols = -c(Geneid, Chr, Length, product),
        names_to = "Sample", values_to = "tpm_rna")

phage_counts_long <- phage_gene_counts %>%
    rownames_to_column("Geneid") %>%
    pivot_longer(cols = -c(Geneid, Chr, Length, product),
        names_to = "Sample", values_to = "tpm_dna")

phage_merged <- phage_expression_long %>%
    left_join(phage_counts_long %>% dplyr::select(Geneid, Sample, tpm_dna),
        by = c("Geneid", "Sample"))

## ---- Join phynteny annotations + metadata ----
phage_meta_merged <- phage_merged %>%
    left_join(sample_metadata, by = "Sample") %>%
    left_join(phynteny_map, by = c("Geneid" = "featurecounts_id")) %>%
    dplyr::rename(bakta_product = product.x,
        phynteny_product = product.y) %>%
    mutate(Bin = str_extract(Chr, "Bin\\.\\d+"))

## ---- Variance filter ----
phage_initial_filter <- phage_meta_merged %>%
    dplyr::filter(!(tpm_dna == 0 & tpm_rna > 0)) %>%
    group_by(Geneid) %>%
    summarise(var_D = var(tpm_dna, na.rm = TRUE),
        var_R = var(tpm_rna, na.rm = TRUE)) %>%
    dplyr::filter(var_D > 0, var_R > 0)

phage_intermediate_filter <- phage_meta_merged %>%
    dplyr::filter(Geneid %in% phage_initial_filter$Geneid)

## ---- CLR transformation ----
phage_clr_data <- phage_intermediate_filter %>%
    dplyr::filter(tpm_rna > 0, tpm_dna > 0) %>%
    group_by(Mouse, Fraction) %>%
    mutate(clr_rna = log2(tpm_rna) - mean(log2(tpm_rna)),
        clr_dna = log2(tpm_dna) - mean(log2(tpm_dna))) %>%
    ungroup()

## ---- Post-CLR filtering ----
phage_clr_filtered <- phage_clr_data %>%
    group_by(Geneid, Fraction) %>%
    dplyr::filter(n() >= 2) %>%
    group_by(Geneid) %>%
    dplyr::filter(n_distinct(Fraction) >= 2) %>%
    ungroup()

## ---- Join host bin CLR DNA covariate ----
bin_clr_long <- bin_clr_abundance %>%
    rownames_to_column("Bin") %>%
    pivot_longer(cols = -Bin, names_to = "Sample", values_to = "bin_clr_dna")

phage_model_data <- phage_clr_filtered %>%
    left_join(bin_clr_long %>% dplyr::select(Sample, Bin, bin_clr_dna),
        by = c("Sample", "Bin"))

## ---- Run phage DEG models ----
phage_data_list <- phage_model_data %>%
    group_by(Geneid) %>%
    dplyr::filter(n_distinct(Fraction) >= 2) %>%
    dplyr::filter(n() >= 4) %>%
    ungroup() %>%
    split(., .$Geneid)

cat("Phage genes to model:", length(phage_data_list), "\n")

phage_DEGs <- run_lmer_parallel(
    data_list       = phage_data_list,
    formula         = clr_rna ~ bin_clr_dna + Fraction + (1|Mouse),
    emmeans_term    = ~ Fraction,
    fdr_group       = "contrast",
    singular_action = "flag_keep",
    output_cols     = list(
        Geneid = quote(Geneid[1]),
        phynteny_product = quote(phynteny_product[1]),
        viral_function   = quote(viral_function[1]),
        phrog            = quote(phrog[1]),
        provirus         = quote(provirus[1]),
        Bin              = quote(Bin[1])
    )
)

## ---- Genome maps ----
phage_deg_joined <- pharokka_cds %>%
    dplyr::select(gene, start, stop, strand, contig, annot) %>%
    dplyr::rename(locus_tag = gene, provirus = contig) %>%
    dplyr::left_join(
        phage_annotations_combined %>% dplyr::select(locus_tag, combined_category),
        by = "locus_tag") %>%
    dplyr::left_join(
        phynteny_map %>% dplyr::select(locus_tag, featurecounts_id),
        by = "locus_tag") %>%
    dplyr::left_join(
        phage_DEGs %>%
            dplyr::filter(contrast == "Positive - Negative") %>%
            dplyr::select(Geneid, estimate, p_adjust, significant),
        by = c("featurecounts_id" = "Geneid")) %>%
    dplyr::mutate(
        Bin = str_extract(provirus, "Bin\\.\\d+"),
        sig_label = case_when(p_adjust < 0.05 ~ "*"))

region_colours <- c(
    "Lysogeny"    = "#DDCC77",
    "Replication" = "#AA3377",
    "Structural"  = "#999933",
    "Lysis"       = "#44AA99",
    "Moron/AMG"   = "#CC6677")

## ---- Genome map plotting function ----
plot_phage_map <- function(provirus_id, show_legend = TRUE) {
    d <- phage_deg_joined %>%
        dplyr::filter(provirus == provirus_id) %>%
        mutate(gene_xmin = pmin(start, stop),
            gene_xmax = pmax(start, stop),
            mid = (start + stop) / 2,
            tested = !is.na(estimate),
            region = case_when(
                combined_category %in% c("integration and excision",
                    "transcription regulation") ~ "Lysogeny",
                combined_category == "DNA, RNA and nucleotide metabolism" ~ "Replication",
                combined_category %in% c("head and packaging", "tail",
                    "connector") ~ "Structural",
                combined_category == "lysis" ~ "Lysis",
                TRUE ~ NA_character_),
            xmin_rect = pmin(start, stop),
            xmax_rect = pmax(start, stop),
            y_pos = 1) %>%
        arrange(pmin(start, stop))

    p <- ggplot(d, aes(xmin = gene_xmin, xmax = gene_xmax, y = y_pos,
            forward = strand == "+")) +
        geom_gene_arrow(data = dplyr::filter(d, !tested),
            fill = "grey70", colour = NA,
            arrowhead_height  = unit(2, "mm"),
            arrowhead_width   = unit(1.5, "mm"),
            arrow_body_height = unit(1.5, "mm")) +
        geom_gene_arrow(data = dplyr::filter(d, tested),
            aes(fill = estimate), colour = NA,
            arrowhead_height  = unit(2, "mm"),
            arrowhead_width   = unit(1.5, "mm"),
            arrow_body_height = unit(1.5, "mm")) +
        geom_text(aes(x = mid, y = y_pos, label = sig_label),
            size = 7 / .pt, family = "Arial", vjust = -0.5,
            inherit.aes = FALSE) +
        scale_fill_gradient2(
            low = col_neg, mid = "grey90", high = col_pos,
            midpoint = 0, limits = c(-3, 3), oob = scales::squish,
            na.value = "grey70",
            name = "Estimate",
            guide = if (show_legend) guide_colorbar(
                direction = "vertical", title.position = "top",
                barwidth = unit(2, "mm"), barheight = unit(15, "mm"),
                order = 1) else "none") +
        new_scale_fill() +
        geom_rect(data = dplyr::filter(d, !is.na(region)),
            aes(xmin = xmin_rect, xmax = xmax_rect,
                ymin = 0.85, ymax = 0.90, fill = region),
            inherit.aes = FALSE) +
        scale_fill_manual(values = region_colours, name = "Region",
            guide = if (show_legend) guide_legend(
                order = 2, keywidth = unit(2, "mm"), keyheight = unit(2, "mm"),
                ncol = 1) else "none") +
        coord_cartesian(ylim = c(0.75, 1.15), clip = "off") +
        theme_genes() +
        theme_for_figures(show_labels = TRUE) +
        theme(panel.border = element_rect(colour = "black",
                fill = NA, linewidth = lw_axis),
            axis.text.y      = element_blank(),
            axis.title.y     = element_blank(),
            axis.line.y      = element_blank(),
            axis.ticks.y     = element_blank(),
            legend.position  = if (show_legend) "right" else "none",
            legend.box       = "vertical",
            legend.text      = element_text(size = base_sz, family = "Arial"),
            legend.title     = element_text(size = base_sz, family = "Arial"),
            legend.spacing.y = unit(3, "mm"),
            legend.margin    = ggplot2::margin(0, 0, 0, 0, "mm"),
            plot.margin      = ggplot2::margin(0, 0, 0, 0, "mm")) +
        labs(title = NULL, x = "Genomic position (bp)")

    p
}

## ---- Build genome maps ----
map_bin44   <- plot_phage_map("Bin.44_20|provirus_2_61998", show_legend = FALSE)
map_bin61_5 <- plot_phage_map("Bin.61_5|provirus_61_52543", show_legend = FALSE)
map_bin76   <- plot_phage_map("Bin.76_1|provirus_283925_325418", show_legend = TRUE)

phage_maps_selected <- wrap_plots(
    list(map_bin44, map_bin61_5, map_bin76),
    ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "none",
        legend.text      = element_text(size = base_sz, family = "Arial"),
        legend.title     = element_text(size = base_sz, family = "Arial"),
        legend.spacing.y = unit(3, "mm"),
        legend.margin    = ggplot2::margin(0, 0, 0, 0, "mm"))

## ---- Save ----
panel_w <- 95
panel_h <- 19 * 3

ggsave(file.path(fig_dir, "phage_no_legends.svg"), phage_maps_selected,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")
