# =============================================================================
# Script:       05_bin_deg_analysis.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Bin-level differential gene expression analysis using
#               parallelized lmer models on CLR-transformed gene counts,
#               controlling for DNA abundance. Volcano plot of Positive
#               vs Negative contrast.
# Dependencies: data_structuring.R, helpers.R, plotting.R
# Input:        clr_genes_filtered (from data_structuring.R)
#               Or cached: Bin_DEGs_08-05-2026.csv
# Output:       Bin_DEGs_08-05-2026.csv (model results)
#               A .svg of bin-level volcano plot. Figure 2 panel A
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

## ---- Data structuring ----
bin_list <- clr_genes_filtered %>%
    dplyr::select(Geneid, Bin, product, Mouse, Fraction, clr_rna, clr_dna) %>%
    mutate(bin_gene = paste0(Bin, "__", Geneid)) %>%
    group_by(Bin) %>%
    group_split() %>%
    setNames(sort(unique(clr_genes_filtered$Bin)))

cat("Total bin datasets to model:", length(bin_list), "\n")

cat("Genes per bin:\n")
purrr::map_dfr(bin_list, ~tibble(Bin = unique(.x$Bin),
    n_genes = n_distinct(.x$bin_gene))) %>%
    arrange(desc(n_genes)) %>%
    print(n = 100)

## ---- Run models ----
bin_degs_file <- file.path(data_dir, "Bin_DEGs_08-05-2026.csv")

if(file.exists(bin_degs_file)){
    cat("Loading existing DEGs file from:", bin_degs_file, "\n")
    DEGs <- read.csv(bin_degs_file)
} else {
    cat("No existing DEGs file, running models \n")
    DEGs <- run_lmer_parallel(data_list = bin_list, inner_split = "bin_gene",
        formula = clr_rna ~ clr_dna + Fraction + (1|Mouse),
        emmeans_term = ~ Fraction, fdr_group = c("contrast", "Bin"),
        singular_action = "flag_keep", output_cols = list(
            Geneid = quote(unique(Geneid)[1]), Bin = quote(unique(Bin)[1]),
            product = quote(unique(product)[1]), n_obs = quote(nrow(d))))

    write.csv(DEGs, bin_degs_file, row.names = FALSE)
}

## ---- Collect results ----
cat("DEGs per contrast:\n")
DEGs %>%
    dplyr::filter(significant) %>%
    count(contrast) %>%
    arrange(desc(n)) %>%
    print()

cat("\nDEGs per bin:\n")
DEGs %>%
    dplyr::filter(significant) %>%
    count(Bin) %>%
    arrange(desc(n))

cat("\nDEGs per bin per contrast:\n")
DEGs %>%
    filter(significant) %>%
    count(Bin, contrast) %>%
    arrange(Bin, contrast)

## ---- Volcano plot ----
pos_neg_volcano_degs <- DEGs %>%
    dplyr::filter(contrast == "Positive - Negative") %>%
    mutate(neg_log10_fdr = -log10(p_adjust), direction = case_when(
        significant & estimate > 0 ~ "Up", significant & estimate < 0 ~ "Down",
        TRUE ~ "Not significant"))

pos_neg_deg_volcano <- ggplot(pos_neg_volcano_degs, aes(x = estimate,
        y = neg_log10_fdr)) +
    geom_point_rast(aes(color = direction), alpha = 0.5, size = pt_small,
        raster.dpi = 1200) +
    geom_vline(xintercept = 0, linetype = "dashed",
        color = "grey50", linewidth = lw_sm) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
        color = "grey50", linewidth = lw_sm) +
    scale_color_manual(values = c("Up"   = col_pos,
        "Down" = col_neg, "Not significant" = "grey70")) +
    scale_x_continuous(limits = c(-20, 20)) +
    theme_for_figures(show_labels = TRUE) +
    labs(x = "Δ CLR(RNA)",
        y = "-log10(FDR)")
pos_neg_deg_volcano

## ---- Save ----
panel_w <- 41.5
panel_h <- 46.2

ggsave(file.path(fig_dir, "Bin_volcano_posneg.svg"), pos_neg_deg_volcano,
    device = svglite::svglite, width = panel_w, height = panel_h, units = "mm")
