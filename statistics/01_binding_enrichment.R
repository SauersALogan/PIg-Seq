# =============================================================================
# Script:       01_binding_enrichment.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Tests whether IgA-bound bacteria are enriched in the
#               positive fraction using beta-regression
# Dependencies: data_structuring.R, plotting.R
# Input:        metadata (from data_structuring.R)
# Output:       A .svg of binding enrichment. Figure 1 panel C
# =============================================================================
library(glmmTMB)
library(emmeans)
library(car)
library(ggplot2)
library(scales)
library(svglite)

# ---- Paths (update for your system) ----
data_dir <- "/path/to/data"
fig_dir  <- "/path/to/figure_output_dir"

source("statistics/functions/plotting.R")
source("statistics/functions/data_structuring.R")

# =====================================================================
# Beta-regression model
# =====================================================================
binding_model <- glmmTMB(IgA_bound ~ Fraction + (1|Sex), 
    data = metadata, family = beta_family())

Anova(binding_model, test = c("Chisq"), type = c("III"), component = "cond")

binding_emmeans <- as.data.frame(emmeans(binding_model, 
    ~Fraction, type = "response", component = "cond"))

pairs_result <- pairs(emmeans(binding_model, ~Fraction), adjust = "fdr")
pairs_df <- as.data.frame(pairs_result)

# =====================================================================
# Plotting
# =====================================================================
Binding_plot <- ggplot(binding_emmeans, 
        aes(x = Fraction, y = response, color = Fraction)) +
    geom_point(data = metadata, aes(x = Fraction, y = IgA_bound, color = Fraction),
        alpha = 0.3, size = pt_small,
        position = position_jitter(width = 0.13, height = 0, seed = 1)) +
    geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
        linewidth = lw_main, size = pt_mean / 6, fatten = 2) +
    theme_for_figures(show_labels = TRUE) +
    scale_x_discrete(limits = c("Native", "Negative", "Positive"),
        labels = c("Unsorted", "Negative", "Positive"), name = "Fraction",
        expand = expansion(add = 0.65)) +
    scale_y_continuous(limits = c(0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0),
        expand = expansion(mult = c(0.02, 0.02)),
        oob = scales::oob_keep) +
    scale_color_manual(values = c("Native"   = col_native,
                                  "Negative" = col_neg,
                                  "Positive" = col_pos)) +
    labs(y = "Proportion IgA Bound")
Binding_plot

# =====================================================================
# Save
# =====================================================================
panel_w <- 41.37
panel_h <- 37

ggsave(file.path(fig_dir, "binding_enrichment.svg"), Binding_plot,
    width = panel_w, height = panel_h, units = "mm",
    device = svglite::svglite, dpi = 600)
