# =============================================================================
# Script:       02_dna_abundance.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Metagenomic (DNA) community analysis — PCoA ordination
#               with PERMANOVA and Shannon alpha diversity with pairwise
#               contrasts
# Dependencies: data_structuring.R, plotting.R
# Input:        Data (phyloseq object from data_structuring.R)
# Output:       A .svg of DNA PCoA ordination. Figure 1 panel F
#               A .svg of DNA Shannon diversity. Figure 1 panel G
# =============================================================================

# ---- Paths (update for your system) ----
fig_dir <- "/path/to/figure_output_dir"

library(microbiome)
library(vegan)
library(lmerTest)
library(emmeans)
library(rstatix)
library(ggplot2)
library(svglite)

source("statistics/functions/plotting.R")
source("statistics/functions/data_structuring.R")

# ---- DNA ordination ----
Data_clr <- microbiome::transform(Data, "clr")
frac.dist <- phyloseq::distance(Data_clr, "euclidean", type = "samples")

ord.data <- as(sample_data(Data_clr), "data.frame")
ord.data$Fraction <- as.factor(ord.data$Fraction)

PERMAnova <- adonis2(frac.dist ~ Fraction, 
    data = ord.data, perm = 10000)
PERMAnova

ordination <- ordinate(Data_clr, method = "MDS", 
    distance = frac.dist)

ord.plot <- plot_ordination(Data_clr, ordination, 
    color = "Fraction", axes = c(1, 2))
ord.plot$layers <- ord.plot$layers[-1]

ord.plot <- ord.plot +
    geom_point(size = pt_small * 1.5) +
    theme_for_figures(show_labels = TRUE) +
    theme(plot.margin = ggplot2::margin(0, 0, 0, 5, "mm")) +
    scale_color_manual(name = "Fraction",
        limits = c("Native", "Negative", "Positive"),
        labels = c("Unsorted", "Negative", "Positive"),
        values = c("Native"   = col_native,
                   "Negative" = col_neg,
                   "Positive" = col_pos)) +
    labs(y = "PCoA 2 [15.2%]", x = "PCoA 1 [37.5%]") +
    scale_y_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 5),
        expand = expansion(mult = c(0.05, 0.05))) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)))
ord.plot

panel_w <- 84.33
panel_h <- 40.5

ggsave(file.path(fig_dir, "DNA_ordination.svg"), ord.plot,
    width = panel_w, height = panel_h, units = "mm",
    device = svglite::svglite, dpi = 600)

# ---- Shannon alpha diversity ----
abundance_shannons <- estimate_richness(Data, measures = "Shannon")
abundance_shannons$Fraction <- factor(sample_data(Data)$Fraction, 
    levels = c("Native", "Negative", "Positive"))
abundance_shannons$Mouse <- sample_data(Data)$Mouse

model_shannon <- lmerTest::lmer(Shannon ~ Fraction + (1|Mouse), 
    data = abundance_shannons)
anova(model_shannon)

abundance_pairwise_contrasts <- abundance_shannons %>%
    emmeans_test(Shannon ~ Fraction, model = model_shannon, 
        p.adjust.method = "fdr") %>%
    add_xy_position(x = "Fraction")

abundance.shannon.plot <- ggplot(abundance_shannons,
        aes(x = Fraction, y = Shannon, color = Fraction)) +
    geom_boxplot(linewidth = lw_sm, outlier.shape = NA) +
    geom_jitter(size = pt_small, width = 0.1, height = 0) +
    scale_color_manual(name = "Fraction",
        limits = c("Native", "Negative", "Positive"),
        labels = c("Unsorted", "Negative", "Positive"),
        values = c("Native"   = col_native,
                   "Negative" = col_neg,
                   "Positive" = col_pos)) +
    scale_x_discrete(limits = c("Native", "Negative", "Positive"),
        labels = c("Unsorted", "Negative", "Positive")) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1),
        expand = expansion(mult = c(0.02, 0.02))) +
    labs(y = "Shannon Diversity", x = "Fraction") +
    theme_for_figures(show_labels = TRUE)
abundance.shannon.plot

panel_w <- 41.37
panel_h <- 37

ggsave(file.path(fig_dir, "DNA_shannon.svg"), abundance.shannon.plot,
    width = panel_w, height = panel_h, units = "mm",
    device = svglite::svglite, dpi = 600)
