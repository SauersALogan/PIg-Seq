# =============================================================================
# Script:       03_transcriptional_activity.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Metatranscriptomic (RNA) community analysis — PCoA ordination
#               with PERMANOVA and Shannon alpha diversity with pairwise
#               contrasts
# Dependencies: data_structuring.R, plotting.R
# Input:        rna_data (phyloseq object from data_structuring.R)
# Output:       A .svg of RNA PCoA ordination. Figure 1 panel D
#               A .svg of RNA Shannon diversity. Figure 1 panel E
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

# ---- RNA ordination ----
rna_clr <- microbiome::transform(rna_data, "clr")
frac.dist.rna <- phyloseq::distance(rna_clr, "euclidean", type = "samples")

ord.data.rna <- as(sample_data(rna_clr), "data.frame")
ord.data.rna$Fraction <- as.factor(ord.data.rna$Fraction)

PERMAnova <- adonis2(frac.dist.rna ~ Fraction, 
    data = ord.data.rna, perm = 10000)
PERMAnova

ordination.rna <- ordinate(rna_clr, method = "MDS", 
    distance = frac.dist.rna)

rna.ord.plot <- plot_ordination(rna_clr, ordination.rna, 
    color = "Fraction", axes = c(1, 2))
rna.ord.plot$layers <- rna.ord.plot$layers[-1]

rna.ord.plot <- rna.ord.plot +
    geom_point(size = pt_small * 1.5) +
    theme_for_figures(show_labels = TRUE) +
    theme(plot.margin = ggplot2::margin(0, 0, 0, 5, "mm")) +
    scale_color_manual(name = "Fraction",
        limits = c("Native", "Negative", "Positive"),
        labels = c("Unsorted", "Negative", "Positive"),
        values = c("Native"   = col_native,
                   "Negative" = col_neg,
                   "Positive" = col_pos)) +
    labs(y = "PCoA 2 [17.5%]", x = "PCoA 1 [23.4%]") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))
rna.ord.plot

panel_w <- 84.33
panel_h <- 40.5

ggsave(file.path(fig_dir, "RNA_ordination.svg"), rna.ord.plot,
    width = panel_w, height = panel_h, units = "mm",
    device = svglite::svglite, dpi = 600)

# ---- Shannon alpha diversity ----
transcript_shannons <- estimate_richness(rna_data, measures = "Shannon")
transcript_shannons$Fraction <- factor(sample_data(rna_data)$Fraction, 
    levels = c("Native", "Negative", "Positive"))
transcript_shannons$Mouse <- sample_data(rna_data)$Mouse

transcript_model_shannon <- lmerTest::lmer(Shannon ~ Fraction + (1|Mouse), 
    data = transcript_shannons)
anova(transcript_model_shannon)

transcript_pairwise_contrasts <- transcript_shannons %>%
    emmeans_test(Shannon ~ Fraction, model = transcript_model_shannon, 
        p.adjust.method = "fdr") %>%
    add_xy_position(x = "Fraction")

transcript.shannon.plot <- ggplot(transcript_shannons,
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
transcript.shannon.plot

panel_w <- 41.37
panel_h <- 37

ggsave(file.path(fig_dir, "RNA_shannon.svg"), transcript.shannon.plot,
    width = panel_w, height = panel_h, units = "mm",
    device = svglite::svglite, dpi = 600)

