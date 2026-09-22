# =============================================================================
# Script:       15_iga_scores.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  IgA binding score estimation via probability ratio method
#               (Jackson et al. 2021, IgAScores), rescaled to 0-1 beta
#               distribution. Beta regression with emmeans testing each
#               taxon against the midpoint (0.5). Produces per-taxon
#               binding classification (Preferentially Bound/Unbound/No
#               Preference) and a dotplot of significant taxa.
# Dependencies: data_structuring.R, plotting.R
# Input:        pos_taxa, neg_taxa, pos_binding, neg_binding, Taxa, Bin_map,
#               sample_metadata (from data_structuring.R)
# Output:       bin_binding_class (binding classifications)
#               A .svg of IgA binding scores. Figure 5 Panel A
# =============================================================================

# ---- Paths (update for your system) ----
fig_dir  <- "/path/to/figure_output_dir"
data_dir <- "/path/to/data"

library(IgAScores)
library(glmmTMB)
library(emmeans)
library(ggplot2)
library(svglite)

source("statistics/functions/plotting.R")
source("statistics/functions/data_structuring.R")

## ---- IgA scores calculation ----
taxon_scores <- igascores(posabunds = pos_taxa, negabunds = neg_taxa,
    possizes = pos_binding, negsizes = neg_binding,
    method = "probratio",
    scaleratio = TRUE, nazeros = TRUE, pseudo = 1e-9)

beta_scores <- (taxon_scores + 1) / 2

beta_long <- beta_scores %>%
    rownames_to_column("Taxon") %>%
    pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "IgA_score")

## ---- Model data ----
binding_errors <- sample_metadata %>%
    dplyr::select(Mouse, Fraction, Binding_error) %>%
    pivot_wider(names_from = Fraction, values_from = Binding_error,
        names_prefix = "Error_")

iga_sample_data <- sample_metadata %>%
    dplyr::filter(Fraction == "Native") %>%
    dplyr::select(Mouse, Sex, Age, Cage, Batch) %>%
    left_join(binding_errors, by = "Mouse")

model_data <- beta_long %>%
    left_join(iga_sample_data, by = c("Sample" = "Mouse"))

## ---- Beta regression ----
beta_model <- glmmTMB(IgA_score ~ Taxon + (1|Sample),
    data = model_data, family = beta_family())

drop1(beta_model, test = "Chi")

em_taxon <- emmeans(beta_model, ~ Taxon, type = "response")
em_summary <- as.data.frame(summary(em_taxon))
test_results <- test(em_taxon, null = qlogis(0.5))

taxon_results <- as.data.frame(summary(test_results))
taxon_results$p.adj <- p.adjust(test_results$p.value, method = "fdr")

taxon_results <- taxon_results %>%
    left_join(em_summary %>%
        dplyr::select(Taxon, asymp.LCL, asymp.UCL), by = "Taxon")

## ---- Binding classification ----
bin_binding_class <- taxon_results %>%
    mutate(mean_score = (response * 2) - 1,
        binding_class_beta = case_when(
            p.adj < 0.05 & mean_score > 0 ~ "Preferentially Bound",
            p.adj < 0.05 & mean_score < 0 ~ "Preferentially Unbound",
            TRUE ~ "No Preference")) %>%
    dplyr::select(Taxon, binding_class_beta) %>%
    left_join(Bin_map %>%
        mutate(Bin = paste0("Bin.", Bin)) %>%
        dplyr::select(Identifier, Bin),
        by = c("Taxon" = "Identifier"))

## ---- Plotting ----
sig_taxa <- taxon_results %>%
    dplyr::filter(p.adj < 0.05)

sig_taxa$mean_score <- (sig_taxa$response * 2) - 1
sig_taxa$lower <- (sig_taxa$asymp.LCL * 2) - 1
sig_taxa$upper <- (sig_taxa$asymp.UCL * 2) - 1

sig_taxa$Species <- Taxa[sig_taxa$Taxon, "Species"]
sig_taxa$Sign <- ifelse(sig_taxa$mean_score >= 0, "Positive", "Negative")

IgA_scores_plot <- ggplot(sig_taxa,
        aes(y = reorder(Taxon, mean_score), x = mean_score)) +
    geom_point(aes(color = Sign), size = pt_mean) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.3,
        linewidth = lw_sm) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = lw_sm) +
    labs(y = "Taxa", x = "Binding probability") +
    scale_color_manual(values = c("Positive" = "#008B8B",
        "Negative" = "#D55E00")) +
    theme_for_figures(show_labels = TRUE) +
    scale_y_discrete(labels = setNames(gsub("_",
        " ", sig_taxa$Species), sig_taxa$Taxon)) +
    scale_x_continuous(limits = c(-1, 1))
IgA_scores_plot

## ---- Save ----
ggsave(file.path(fig_dir, "IgA_scores.svg"), IgA_scores_plot,
    device = svglite::svglite, width = 100, height = 100, units = "mm",
    dpi = 600)
