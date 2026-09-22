# =============================================================================
# Script:       16_binding_rf_classifier.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Three-class random forest classifier predicting per-Mouse x
#               Bin IgA binding category (Low < 0.25, Medium 0.25-0.75,
#               High > 0.75) from three independent feature sets: DNA
#               presence/absence, RNA Positive fraction expression, and RNA
#               Negative fraction expression. 100 bootstrap iterations with
#               bin-level holdout and inverse class weighting. Macro-averaged
#               AUC with pooled ROC curves.
# Dependencies: data_structuring.R, plotting.R
# Input:        DNA_gene_count_controlled_subsampled.csv (https://doi.org/10.48420/33951214)
#               RNA_gene_count_controlled_subsampled.csv (https://doi.org/10.48420/33951214)
#               pos_taxa, neg_taxa, pos_binding, neg_binding (from data_structuring.R)
# Output:       A .svg of AUC boxplot. Figure 5 Panel B
#               A .svg of pooled ROC curves. Figure 5 Panel C
#               A .svg of IgA score distribution. Figure 5 Panel A
# =============================================================================

# ---- Paths (update for your system) ----
fig_dir  <- "/path/to/figure_output_dir"
data_dir <- "/path/to/data"

## ---- Libraries ----
library(tidyverse)
library(ranger)
library(pROC)
library(furrr)
library(future)
library(progressr)
library(IgAScores)
library(glmmTMB)
library(emmeans)
library(pracma)
library(stringr)
library(MuMIn)
 
options(contrasts = c("contr.sum", "contr.poly"))
filter <- dplyr::filter
 
source("statistics/functions/plotting.R")
 
## ---- Data loading ----
dna_subsample_counts <- read.csv(file.path(data_dir, "DNA_gene_count_controlled_subsampled.csv"), sep = ",", row.names = 1)
rna_subsample_counts <- read.csv(file.path(data_dir, "RNA_gene_count_controlled_subsampled.csv"), sep = ",", row.names = 1)
 
metadata <- read.csv(file.path(data_dir, "metadata_30-09-2025.csv"), sep = ",", row.names = 1)
Taxa     <- read.csv(file.path(data_dir, "Taxonomy_30-09-2025.csv"), sep = ",", row.names = 1)
Bin_map  <- read.csv(file.path(data_dir, "Bin_number_to_taxonomy.csv"), sep = ",")
 
pos_taxa <- read.csv(file.path(data_dir, "Pos_rel_abund_30-09-2025.csv"), sep = ",", row.names = 1)
neg_taxa <- read.csv(file.path(data_dir, "Neg_rel_abund_30-09-2025.csv"), sep = ",", row.names = 1)
 
pos_binding <- c(0.485, 0.325, 0.2620, 0.231, 0.248, 0.285, 0.379, 0.351)
neg_binding <- c((1 - 0.485), (1 - 0.325), (1 - 0.2620), (1 - 0.231),
                 (1 - 0.248), (1 - 0.285), (1 - 0.379), (1 - 0.351))
 
## ---- Factor and type coercion ----
Bin_map$Bin  <- as.character(Bin_map$Bin)
metadata$Fraction <- factor(metadata$Fraction,
    levels = c("Positive", "Negative", "Native"))
metadata$Mouse <- as.character(metadata$Mouse)
metadata$Batch <- as.factor(metadata$Batch)
metadata$Sex   <- as.factor(metadata$Sex)
metadata$Cage  <- as.factor(metadata$Cage)
 
sample_metadata <- as.data.frame(metadata)
sample_metadata$Sample <- rownames(sample_metadata)
 
## ---- Gene-level merge, filter, and product-level aggregation ----
 
## ---- Long-format gene expression and counts ----
expression_long <- rna_subsample_counts %>%
    rownames_to_column("Geneid") %>%
    pivot_longer(cols = -c(Geneid, Chr, Length, product),
        names_to = "Sample", values_to = "tpm_rna")
 
counts_long <- dna_subsample_counts %>%
    rownames_to_column("Geneid") %>%
    pivot_longer(cols = -c(Geneid, Chr, Length, product),
        names_to = "Sample", values_to = "tpm_dna")
 
## ---- Merge RNA and DNA ----
merged <- expression_long %>%
    left_join(counts_long %>% dplyr::select(Geneid, Sample, tpm_dna),
        by = c("Geneid", "Sample"))
 
meta_merged <- merged %>%
    left_join(sample_metadata, by = "Sample")
 
cat("Total genes:", n_distinct(merged$Geneid), "\n")
 
## ---- Initial gene filter ----
# Remove genes where RNA present but no DNA, and no variance
initial_gene_filter <- meta_merged %>%
    filter(!(tpm_dna == 0 & tpm_rna > 0)) %>%
    group_by(Geneid) %>%
    summarise(var_D = var(tpm_dna, na.rm = TRUE),
        var_R = var(tpm_rna, na.rm = TRUE)) %>%
    filter(var_D > 0, var_R > 0)
 
intermediate_data_filtered <- meta_merged %>%
    filter(Geneid %in% initial_gene_filter$Geneid) %>%
    mutate(Bin = str_extract(Geneid, "Bin\\.\\d+"))
 
cat("Genes after first filter:",
    n_distinct(intermediate_data_filtered$Geneid), "\n")
 
## ---- Product-level aggregation (pseudobulk) ----
pseudobulk_data <- intermediate_data_filtered %>%
    mutate(Bin = str_extract(Geneid, "Bin\\.\\d+"),
        product = tolower(as.character(product)),
        product = gsub("%2c", ",",  product, ignore.case = TRUE),
        product = gsub("%2e", ".",  product, ignore.case = TRUE),
        product = gsub("%28", "(",  product, ignore.case = TRUE),
        product = gsub("%29", ")",  product, ignore.case = TRUE),
        product = gsub("%27", "'",  product, ignore.case = TRUE),
        product = str_squish(product),
        product = str_trim(product)) %>%
    group_by(Mouse, Fraction, Bin, product) %>%
    summarise(tpm_rna = sum(tpm_rna, na.rm = TRUE),
        tpm_dna = sum(tpm_dna, na.rm = TRUE),
        n_genes = n_distinct(Geneid), .groups = "drop")
 
cat("Products (including hypotheticals):",
    n_distinct(pseudobulk_data$product), "\n")
 
## ---- IgA score calculation and label assignment ----
 
## ---- IgA scores ----
taxon_scores <- igascores(posabunds = pos_taxa, negabunds = neg_taxa,
    possizes = pos_binding, negsizes = neg_binding,
    method = "probratio", scaleratio = TRUE, nazeros = TRUE, pseudo = 1e-9)
 
beta_scores <- (taxon_scores + 1) / 2
 
beta_long <- beta_scores %>%
    rownames_to_column("Taxon") %>%
    pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "IgA_score")
 
## ---- Beta-regression binding classification (bin-level) ----
binding_errors <- sample_metadata %>%
    dplyr::select(Mouse, Fraction, Binding_error) %>%
    pivot_wider(names_from = Fraction, values_from = Binding_error,
        names_prefix = "Error_")
 
sample_data <- sample_metadata %>%
    dplyr::filter(Fraction == "Native") %>%
    dplyr::select(Mouse, Sex, Mouse, Age, Cage, Batch) %>%
    left_join(binding_errors, by = "Mouse")
 
model_data <- beta_long %>%
    left_join(sample_data, by = c("Sample" = "Mouse"))
 
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
 
## ---- Per Mouse x Bin labels from raw IgA score (for RF) ----
sample_binding_join <- beta_long %>%
    left_join(Bin_map %>%
                  mutate(Bin = paste0("Bin.", Bin)) %>%
                  dplyr::select(Identifier, Bin),
              by = c("Taxon" = "Identifier")) %>%
    mutate(binding_class_sample = case_when(
        IgA_score <  0.25 ~ "Low",
        IgA_score >  0.75 ~ "High",
        TRUE              ~ "Medium"
    )) %>%
    dplyr::select(Mouse = Sample, Bin, binding_class_sample) %>%
    filter(!is.na(Bin)) %>%
    mutate(binding_class_sample = factor(binding_class_sample,
                                         levels = c("Low", "Medium", "High")))
 
cat("Label distribution:\n")
table(sample_binding_join$binding_class_sample)
 
sample_binding_join %>%
    group_by(Bin) %>%
    summarise(n_classes = n_distinct(binding_class_sample), .groups = "drop") %>%
    count(n_classes)
 
## ---- Feature matrix construction ----
 
## ---- Common filtered base ----
posneg_base <- pseudobulk_data %>%
    filter(Fraction %in% c("Positive", "Negative"), tpm_rna > 0, tpm_dna > 0) %>%
    filter(!grepl("hypothetical protein", product, ignore.case = TRUE)) %>%
    filter(!grepl("ribosomal rna|rrna|16s|23s|5s", product, ignore.case = TRUE))
 
## ---- Shared row structure and feature set ----
sample_index <- posneg_base %>%
    dplyr::select(Mouse, Bin, Fraction) %>%
    distinct() %>%
    group_by(Mouse, Bin) %>%
    filter(n_distinct(Fraction) == 2) %>%
    ungroup() %>%
    dplyr::select(Mouse, Bin) %>%
    distinct()
 
cat("Shared row structure:", nrow(sample_index), "Mouse x Bin combinations\n")
 
product_prevalence <- posneg_base %>%
    semi_join(sample_index, by = c("Mouse", "Bin")) %>%
    group_by(Fraction, product) %>%
    summarise(prevalence = n_distinct(paste0(Mouse, Bin)) / nrow(sample_index),
              .groups = "drop") %>%
    group_by(product) %>%
    filter(all(prevalence >= 0.10)) %>%
    ungroup()
 
shared_products <- unique(product_prevalence$product)
cat("Shared feature set:", length(shared_products), "products\n")
 
## ---- 1. DNA presence/absence matrix (genome-level) ----
dna_genome_level <- pseudobulk_data %>%
    filter(product %in% shared_products) %>%
    group_by(Bin, product) %>%
    summarise(present = as.integer(any(tpm_dna > 0)), .groups = "drop") %>%
    pivot_wider(names_from = product, values_from = present, values_fill = 0)
 
dna_presence_rf <- sample_index %>%
    left_join(dna_genome_level, by = "Bin") %>%
    mutate(across(all_of(intersect(shared_products, names(.))),
                  ~replace_na(.x, 0)))
 
cat("DNA matrix:", nrow(dna_presence_rf), "x",
    ncol(dna_presence_rf) - 2, "features\n")
 
## ---- 2. RNA Positive fraction matrix (relative expression) ----
rna_pos_rf <- posneg_base %>%
    filter(Fraction == "Positive") %>%
    semi_join(sample_index, by = c("Mouse", "Bin")) %>%
    filter(product %in% shared_products) %>%
    group_by(Mouse, Bin) %>%
    mutate(rel_expr = tpm_rna / sum(tpm_rna)) %>%
    ungroup() %>%
    dplyr::select(Mouse, Bin, product, rel_expr) %>%
    pivot_wider(names_from = product, values_from = rel_expr, values_fill = 0) %>%
    right_join(sample_index, by = c("Mouse", "Bin")) %>%
    mutate(across(all_of(intersect(shared_products, names(.))),
                  ~replace_na(.x, 0)))
 
cat("RNA Positive matrix:", nrow(rna_pos_rf), "x",
    ncol(rna_pos_rf) - 2, "features\n")
 
## ---- 3. RNA Negative fraction matrix (relative expression) ----
rna_neg_rf <- posneg_base %>%
    filter(Fraction == "Negative") %>%
    semi_join(sample_index, by = c("Mouse", "Bin")) %>%
    filter(product %in% shared_products) %>%
    group_by(Mouse, Bin) %>%
    mutate(rel_expr = tpm_rna / sum(tpm_rna)) %>%
    ungroup() %>%
    dplyr::select(Mouse, Bin, product, rel_expr) %>%
    pivot_wider(names_from = product, values_from = rel_expr, values_fill = 0) %>%
    right_join(sample_index, by = c("Mouse", "Bin")) %>%
    mutate(across(all_of(intersect(shared_products, names(.))),
                  ~replace_na(.x, 0)))
 
cat("RNA Negative matrix:", nrow(rna_neg_rf), "x",
    ncol(rna_neg_rf) - 2, "features\n")
 
## ---- Attach per-sample labels ----
relabel <- function(df) {
    df %>%
        left_join(sample_binding_join, by = c("Mouse", "Bin")) %>%
        filter(!is.na(binding_class_sample))
}
 
dna_labelled     <- relabel(dna_presence_rf)
rna_pos_labelled <- relabel(rna_pos_rf)
rna_neg_labelled <- relabel(rna_neg_rf)
 
cat("DNA labelled:", nrow(dna_labelled), "rows\n")
cat("RNA Pos labelled:", nrow(rna_pos_labelled), "rows\n")
cat("RNA Neg labelled:", nrow(rna_neg_labelled), "rows\n")
 
## ---- Model training with bin-level holdout ----
 
classes_all <- c("Low", "Medium", "High")
 
run_split_sample_label <- function(df, seed, test_prop = 0.20,
                                   max_tries = 500) {
    set.seed(seed)
 
    all_bins <- unique(df$Bin)
    n_test   <- max(2, floor(length(all_bins) * test_prop))
 
    rare          <- names(sort(table(df$binding_class_sample)))[1]
    min_rare_test <- max(2, floor(sum(df$binding_class_sample == rare) *
        test_prop * 0.5))
 
    tries <- 0
    repeat {
        tries <- tries + 1
        test_bins <- sample(all_bins, n_test)
        te <- df %>% filter(Bin %in% test_bins)
        ok <- n_distinct(te$binding_class_sample) == length(classes_all) &&
              sum(te$binding_class_sample == rare) >= min_rare_test
        if (ok) break
        if (tries >= max_tries) { warning("relaxed at seed ", seed); break }
    }
 
    tr <- df %>% filter(!Bin %in% test_bins)
 
    X_train <- tr %>% dplyr::select(-Mouse, -Bin, -binding_class_sample)
    y_train <- droplevels(tr$binding_class_sample)
    X_test  <- te %>% dplyr::select(-Mouse, -Bin, -binding_class_sample)
    y_test  <- te$binding_class_sample
 
    weights <- as.numeric((1 / table(y_train))[as.character(y_train)])
 
    rf <- ranger(x = X_train, y = y_train, num.trees = 500, probability = TRUE,
                 importance = "impurity", case.weights = weights, num.threads = 1)
 
    imp <- tibble::tibble(product = names(rf$variable.importance),
                          importance = rf$variable.importance)
 
    list(y_test = y_test, prob_preds = predict(rf, data = X_test)$predictions,
         importance = imp)
}
 
## ---- AUC computation ----
 
compute_class_auc <- function(fold_result, class_name) {
    y_test <- fold_result$y_test
    probs  <- fold_result$prob_preds
    if (!(class_name %in% colnames(probs))) return(NA_real_)
    binary_actual <- as.integer(y_test == class_name)
    if (length(unique(binary_actual)) < 2) return(NA_real_)
    tryCatch({
        roc_obj <- pROC::roc(response = binary_actual,
            predictor = probs[, class_name], quiet = TRUE)
        as.numeric(pROC::auc(roc_obj))
    }, error = function(e) NA_real_)
}
 
compute_overall_auc <- function(fold_result) {
    aucs <- sapply(classes_all, function(cls)
        compute_class_auc(fold_result, cls))
    mean(aucs, na.rm = TRUE)
}
 
## ---- Run 100 bootstrap seeds ----
 
seeds <- sample.int(1e6, 100)
test_prop_use <- 0.20
 
plan(multisession, workers = 3)
handlers(global = TRUE); handlers("txtprogressbar")
 
## ---- DNA ----
with_progress({p <- progressor(along = seeds)
    dna_full_prob <- future_map(seeds, function(s){
        r <- run_split_sample_label(dna_labelled, s, test_prop_use); p(); r
    }, .options = furrr_options(seed = TRUE))
})
 
## ---- RNA Positive ----
with_progress({p <- progressor(along = seeds)
    rna_pos_full_prob <- future_map(seeds, function(s){
        r <- run_split_sample_label(rna_pos_labelled, s, test_prop_use); p(); r
    }, .options = furrr_options(seed = TRUE))
})
 
## ---- RNA Negative ----
with_progress({p <- progressor(along = seeds)
    rna_neg_full_prob <- future_map(seeds, function(s){
        r <- run_split_sample_label(rna_neg_labelled, s, test_prop_use); p(); r
    }, .options = furrr_options(seed = TRUE))
})
 
## ---- Summary statistics ----
 
build_auc_long <- function(fold_list, group_label) {
    tibble(
        `Overall AUC` = purrr::map_dbl(fold_list, compute_overall_auc),
        `Low AUC` = purrr::map_dbl(fold_list, ~compute_class_auc(.x, "Low")),
        `Medium AUC` = purrr::map_dbl(fold_list, ~compute_class_auc(.x, "Medium")),
        `High AUC` = purrr::map_dbl(fold_list, ~compute_class_auc(.x, "High"))
    ) %>%
        pivot_longer(everything(), names_to = "Metric", values_to = "AUC") %>%
        mutate(Group = group_label)
}
 
auc_plot_data <- bind_rows(
    build_auc_long(dna_full_prob, "DNA"),
    build_auc_long(rna_pos_full_prob, "RNA (Positive)"),
    build_auc_long(rna_neg_full_prob, "RNA (Negative)")
) %>%
    mutate(
        Group = factor(Group, levels = c("DNA", "RNA (Positive)", "RNA (Negative)")),
        Metric = factor(Metric, levels = c("Overall AUC", "High AUC",
            "Low AUC", "Medium AUC"))
    )
 
## ---- Summary with 95% CI ----
auc_summary <- auc_plot_data %>%
    group_by(Metric, Group) %>%
    summarise(
        mean    = mean(AUC, na.rm = TRUE),
        sd      = sd(AUC, na.rm = TRUE),
        median  = median(AUC, na.rm = TRUE),
        ci_low  = quantile(AUC, 0.025, na.rm = TRUE),
        ci_high = quantile(AUC, 0.975, na.rm = TRUE),
        .groups = "drop") %>%
    arrange(Metric, Group)
print(auc_summary, n = Inf, width = Inf)
 
## ---- Pairwise comparison between models ----
cat("\nPairwise Wilcoxon p-values (between models, Overall AUC):\n")
dna_overall     <- purrr::map_dbl(dna_full_prob, compute_overall_auc)
rna_pos_overall <- purrr::map_dbl(rna_pos_full_prob, compute_overall_auc)
rna_neg_overall <- purrr::map_dbl(rna_neg_full_prob, compute_overall_auc)
 
cat("DNA vs RNA (Pos):", wilcox.test(dna_overall, rna_pos_overall)$p.value, "\n")
cat("DNA vs RNA (Neg):", wilcox.test(dna_overall, rna_neg_overall)$p.value, "\n")
cat("RNA (Pos) vs RNA (Neg):", wilcox.test(rna_pos_overall, rna_neg_overall)$p.value, "\n")
 
## ---- Plotting ----
 
## ---- 10a. AUC Boxplot (Overall only) ----
auc_plot_data_overall <- auc_plot_data %>%
    dplyr::filter(Metric == "Overall AUC")
 
auc_boxplot <- ggplot(auc_plot_data_overall,
        aes(x = Group, y = AUC, fill = Group)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85,
        linewidth = lw_sm) +
    geom_jitter(width = 0.1, height = 0, size = pt_small * 0.4, alpha = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
        colour = "grey50", linewidth = lw_sm) +
    scale_fill_manual(values = c("DNA"            = "#117733",
                                 "RNA (Positive)" = "#EE7733",
                                 "RNA (Negative)" = "#882255")) +
    scale_x_discrete(labels = c("DNA" = "DNA",
                                "RNA (Positive)" = "RNA\nPos",
                                "RNA (Negative)" = "RNA\nNeg")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "Overall AUC") +
    theme_for_figures(show_labels = TRUE) +
    theme(legend.position = "none",
          plot.margin = ggplot2::margin(1, 1, 1, 1, "mm"))
auc_boxplot
 
## ---- 10b. Macro-Averaged ROC Curves ----
build_macro_roc <- function(fold_list, model_label) {
    class_rocs <- purrr::map(classes_all, function(cls) {
        all_actual <- purrr::map(fold_list,
            ~as.integer(.x$y_test == cls)) %>% unlist()
        all_probs  <- purrr::map(fold_list,
            ~.x$prob_preds[, cls]) %>% unlist()
        roc_obj <- pROC::roc(response = all_actual,
            predictor = all_probs, quiet = TRUE)
 
        fpr_grid <- seq(0, 1, length.out = 500)
        tpr_interp <- approx(1 - roc_obj$specificities, roc_obj$sensitivities,
                             xout = fpr_grid, rule = 2, ties = mean)$y
        tibble(FPR = fpr_grid, TPR = tpr_interp, Class = cls)
    })
 
    bind_rows(class_rocs) %>%
        group_by(FPR) %>%
        summarise(TPR = mean(TPR), .groups = "drop") %>%
        mutate(
            AUC = round(pracma::trapz(FPR, TPR), 2),
            Model = model_label)
}
 
roc_data <- bind_rows(
    build_macro_roc(dna_full_prob, "DNA"),
    build_macro_roc(rna_pos_full_prob, "RNA (Positive)"),
    build_macro_roc(rna_neg_full_prob, "RNA (Negative)")
)
 
roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR, colour = Model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
        colour = "grey50", linewidth = lw_sm) +
    geom_line(linewidth = lw_main) +
    scale_colour_manual(values = c("DNA"            = "#117733",
                                   "RNA (Positive)" = "#EE7733",
                                   "RNA (Negative)" = "#882255")) +
    labs(x = "False positive rate", y = "True positive rate", colour = NULL) +
    theme_for_figures(show_labels = TRUE) +
    theme(legend.position = "bottom",
          legend.key.size = unit(2.5, "mm"),
          legend.text = element_text(size = base_sz),
          plot.margin = ggplot2::margin(1, 1, 1, 1, "mm"))
roc_plot
 
cat("DNA macro AUC:", unique(roc_data$AUC[roc_data$Model == "DNA"]), "\n")
cat("RNA Pos macro AUC:", unique(roc_data$AUC[roc_data$Model == "RNA (Positive)"]), "\n")
cat("RNA Neg macro AUC:", unique(roc_data$AUC[roc_data$Model == "RNA (Negative)"]), "\n")
 
## ---- 10c. IgA Binding Score Distribution ----
# Raw per-Mouse x Bin scores, colored by binding threshold
iga_plot_data <- beta_long %>%
    left_join(Bin_map %>%
        mutate(Bin = paste0("Bin.", Bin)) %>%
        dplyr::select(Identifier, Bin),
        by = c("Taxon" = "Identifier")) %>%
    filter(!is.na(Bin)) %>%
    mutate(
        IgA_score_rescaled = (IgA_score * 2) - 1,
        binding_color = case_when(
            IgA_score > 0.75 ~ "High",
            IgA_score < 0.25 ~ "Low",
            TRUE             ~ "Medium"
        ),
        binding_color = factor(binding_color,
            levels = c("Low", "Medium", "High"))
    ) %>%
    left_join(bin_binding_class %>%
        dplyr::select(Bin, binding_class_beta), by = "Bin")
 
# Order: group by binding class, then by mean score within class
bin_order <- iga_plot_data %>%
    group_by(Bin, binding_class_beta) %>%
    summarise(mean_score = mean(IgA_score, na.rm = TRUE),
        .groups = "drop") %>%
    mutate(binding_class_beta = factor(binding_class_beta,
        levels = c("Preferentially Unbound", "No Preference",
                   "Preferentially Bound"))) %>%
    arrange(binding_class_beta, mean_score)
 
iga_plot_data$Bin <- factor(iga_plot_data$Bin, levels = bin_order$Bin)
 
IgA_scores_plot <- ggplot(iga_plot_data,
    aes(y = Bin, x = IgA_score)) +
    geom_jitter(aes(color = binding_color),
        height = 0.2, width = 0, size = pt_small, alpha = 0.6) +
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = lw_sm) +
    geom_vline(xintercept = 0.75, linetype = "dotted",
        color = "grey50", linewidth = lw_sm) +
    geom_vline(xintercept = 0.25, linetype = "dotted",
        color = "grey50", linewidth = lw_sm) +
    scale_color_manual(values = c("High" = "#D62728",
        "Low" = "#1F77B4", "Medium" = "grey60")) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(y = NULL, x = "Binding probability") +
    theme_for_figures(show_labels = TRUE) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = ggplot2::margin(1, 1, 1, 1, "mm"))
IgA_scores_plot
 
## ---- Save ----
ggsave(file.path(fig_dir, "auc_boxplot.svg"),
    auc_boxplot + theme(legend.position = "none"),
    device = svglite::svglite, width = 37.5, height = 37.5, units = "mm")
 
ggsave(file.path(fig_dir, "auc_boxplot_legend.svg"),
    auc_boxplot + theme(legend.position = "bottom"),
    device = svglite::svglite, width = 37.5, height = 37.5, units = "mm")
 
ggsave(file.path(fig_dir, "roc_curve.svg"),
    roc_plot + theme(legend.position = "none"),
    device = svglite::svglite, width = 37.5, height = 37.5, units = "mm")
 
ggsave(file.path(fig_dir, "roc_curve_legend.svg"), roc_plot,
    device = svglite::svglite, width = 37.5, height = 37.5, units = "mm")
 
ggsave(file.path(fig_dir, "IgA_scores.svg"),
    IgA_scores_plot + theme(legend.position = "none"),
    device = svglite::svglite, width = 37.5, height = 75, units = "mm")
 
ggsave(file.path(fig_dir, "IgA_scores_legend.svg"), IgA_scores_plot,
    device = svglite::svglite, width = 37.5, height = 75, units = "mm")
 
