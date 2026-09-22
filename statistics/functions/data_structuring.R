# =============================================================================
# Script:       data_structuring.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Loads all data files, coerces types, creates phyloseq objects,
#               and builds core transformed data frames (pseudobulk, CLR,
#               gene-level, family-level) used across analysis notebooks.
# Input:        subsampled_dna_abundance.csv (https://doi.org/10.48420/33951214)
#               subsampled_rna_activity.csv (https://doi.org/10.48420/33951214)
#               metadata_30-09-2025.csv (https://doi.org/10.48420/33951214)
#               Taxonomy_30-09-2025.csv (https://doi.org/10.48420/33951214)
#               Bin_number_to_taxonomy.csv (https://doi.org/10.48420/33951214)
#               DNA_features_controlled.csv (https://doi.org/10.48420/33951214)
#               RNA_features_controlled.csv (https://doi.org/10.48420/33951214)
#               Pos_rel_abund_30-09-2025.csv (https://doi.org/10.48420/33951214)
#               Neg_rel_abund_30-09-2025.csv (https://doi.org/10.48420/33951214)
# Output:       Data, rna_data (phyloseq objects)
#               pseudobulk_data (gene-level) 
#               pseudobulk_filtered (product-level)
#               clr_genes_filtered (gene-level CLR)
#               clr_family_data (family-level CLR)
#               ratio_data (RNA/DNA ratios)
#               sample_metadata, sample_meta
# =============================================================================

library(phyloseq)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(compositions)

## ---- Paths (update for your system) ----
data_dir <- "/path/to/data"

## ---- Load data ----
counts          <- read.csv(file.path(data_dir, "subsampled_dna_abundance.csv"), row.names = 1)
metadata        <- read.csv(file.path(data_dir, "metadata_30-09-2025.csv"), row.names = 1)
Taxa            <- read.csv(file.path(data_dir, "Taxonomy_30-09-2025.csv"), row.names = 1)
Bin_map         <- read.csv(file.path(data_dir, "Bin_number_to_taxonomy.csv"))
gene_counts     <- read.csv(file.path(data_dir, "DNA_features_controlled.csv"), row.names = 1)
gene_expression <- read.csv(file.path(data_dir, "RNA_features_controlled.csv"), row.names = 1)
rna_counts      <- read.csv(file.path(data_dir, "subsampled_rna_activity.csv"), row.names = 1)
pos_taxa        <- read.csv(file.path(data_dir, "Pos_rel_abund_30-09-2025.csv"), row.names = 1)
neg_taxa        <- read.csv(file.path(data_dir, "Neg_rel_abund_30-09-2025.csv"), row.names = 1)
pos_binding     <- c(0.485, 0.325, 0.2620, 0.231, 0.248, 0.285, 0.379, 0.351)
neg_binding     <- c((1 - 0.485), (1 - 0.325), (1 - 0.2620), (1 - 0.231),
                     (1 - 0.248), (1 - 0.285), (1 - 0.379), (1 - 0.351))

## ---- Factor/Character coercion ----
Bin_map$Bin <- as.character(Bin_map$Bin)
metadata$Fraction <- factor(metadata$Fraction,
    levels = c("Positive", "Negative", "Native"))
metadata$Mouse <- as.character(metadata$Mouse)
metadata$Batch <- as.factor(metadata$Batch)
metadata$Sex <- as.factor(metadata$Sex)
metadata$Cage <- as.factor(metadata$Cage)

## ---- Phyloseq objects ----
OTU <- otu_table(counts, taxa_are_rows = TRUE)

tax_mat <- as.matrix(Taxa)
rownames(tax_mat) <- paste0("Bin.", as.integer(
    Bin_map$Bin[match(rownames(tax_mat), Bin_map$Identifier)]))
tax_mat <- tax_mat[!is.na(rownames(tax_mat)), ]
tax_mat <- tax_mat[taxa_names(OTU), ]
TAX <- tax_table(tax_mat)

phyloseq_metadata <- sample_data(metadata)

Data <- phyloseq(OTU, TAX, phyloseq_metadata)

rna_OTU <- otu_table(rna_counts, taxa_are_rows = TRUE)
rna_data <- phyloseq(rna_OTU, TAX, phyloseq_metadata)

## ---- Sample metadata ----
sample_metadata <- as.data.frame(metadata)
sample_metadata$Sample <- rownames(sample_metadata)

sample_meta <- metadata %>%
    rownames_to_column("Sample") %>%
    dplyr::select(Sample, Mouse, Fraction)

## ---- Gene-level merge and initial filter ----
expression_long <- gene_expression %>%
    rownames_to_column("Geneid") %>%
    pivot_longer(cols = -c(Geneid, Chr, Length, product),
        names_to = "Sample", values_to = "tpm_rna")

counts_long <- gene_counts %>%
    rownames_to_column("Geneid") %>%
    pivot_longer(cols = -c(Geneid, Chr, Length, product),
        names_to = "Sample", values_to = "tpm_dna")

merged <- expression_long %>%
    left_join(counts_long %>% dplyr::select(Geneid, Sample, tpm_dna),
        by = c("Geneid", "Sample"))

meta_merged <- merged %>%
    left_join(sample_metadata, by = "Sample")

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

## ---- Product-level CLR ----
pseudobulk_clr <- pseudobulk_data %>%
    filter(tpm_rna > 0, tpm_dna > 0) %>%
    group_by(Mouse, Fraction) %>%
    mutate(clr_rna = log2(tpm_rna) - mean(log2(tpm_rna)),
        clr_dna = log2(tpm_dna) - mean(log2(tpm_dna))) %>%
    ungroup()

pseudobulk_filtered <- pseudobulk_clr %>%
    filter(!grepl("hypothetical protein", product, ignore.case = TRUE)) %>%
    filter(!grepl("ribosomal rna|rrna|16s|23s|5s", product, ignore.case = TRUE))

## ---- Gene-level CLR ----
post_clr_gene_data <- intermediate_data_filtered %>%
    filter(tpm_rna > 0, tpm_dna > 0) %>%
    group_by(Mouse, Fraction) %>%
    mutate(clr_rna = log2(tpm_rna) - mean(log2(tpm_rna)),
        clr_dna = log2(tpm_dna) - mean(log2(tpm_dna))) %>%
    ungroup()

clr_genes_intermediate_filtered <- post_clr_gene_data %>%
    filter(tpm_rna > 0, tpm_dna > 0) %>%
    group_by(Geneid, Fraction) %>%
    filter(n() >= 2) %>%
    group_by(Geneid) %>%
    filter(n_distinct(Fraction) >= 2) %>%
    ungroup()

clr_genes_filtered <- clr_genes_intermediate_filtered %>%
    filter(!grepl("hypothetical protein", product, ignore.case = TRUE)) %>%
    filter(!grepl("ribosomal RNA|rRNA|16S|23S|5S", product, ignore.case = TRUE)) %>%
    mutate(product = tolower(as.character(product)))

cat("Genes after CLR and filtering:",
    n_distinct(clr_genes_filtered$Geneid), "\n")

## ---- Abundance tables ----
dna_abund <- as.data.frame(otu_table(Data)) %>%
    rownames_to_column("Identifier")

rna_abund <- as.data.frame(otu_table(rna_data)) %>%
    rownames_to_column("Identifier")

## ---- RNA/DNA ratios ----
dna_long <- dna_abund %>%
    pivot_longer(cols = -Identifier,
        names_to = "Sample", values_to = "dna_tpm")

rna_long <- rna_abund %>%
    pivot_longer(cols = -Identifier,
        names_to = "Sample", values_to = "rna_tpm")

ratio_data <- dna_long %>%
    left_join(rna_long, by = c("Identifier", "Sample")) %>%
    filter(dna_tpm > 0, rna_tpm > 0) %>%
    mutate(log2_ratio = log2(rna_tpm / dna_tpm)) %>%
    left_join(sample_meta, by = "Sample") %>%
    filter(!is.na(Fraction)) %>%
    left_join(Bin_map %>%
        mutate(Bin = paste0("Bin.", Bin)) %>%
        dplyr::select(Identifier, Bin), by = "Identifier") %>%
    filter(!is.na(Bin))

## ---- Family-level CLR ----
family_lookup <- Taxa %>%
    rownames_to_column("Identifier") %>%
    dplyr::select(Identifier, Family) %>%
    mutate(Family = gsub("f__", "", Family))

bin_to_identifier <- Bin_map %>%
    mutate(Bin_key = paste0("Bin.", as.integer(Bin))) %>%
    dplyr::select(Bin_key, Identifier)

dna_family <- dna_abund %>%
    pivot_longer(cols = -Identifier, names_to = "Sample",
        values_to = "dna_tpm") %>%
    left_join(bin_to_identifier, by = c("Identifier" = "Bin_key")) %>%
    left_join(family_lookup, by = c("Identifier.y" = "Identifier")) %>%
    filter(!is.na(Family)) %>%
    group_by(Sample, Family) %>%
    summarise(dna_tpm = sum(dna_tpm), .groups = "drop")

rna_family <- rna_abund %>%
    pivot_longer(cols = -Identifier, names_to = "Sample",
        values_to = "rna_tpm") %>%
    left_join(bin_to_identifier, by = c("Identifier" = "Bin_key")) %>%
    left_join(family_lookup, by = c("Identifier.y" = "Identifier")) %>%
    filter(!is.na(Family)) %>%
    group_by(Sample, Family) %>%
    summarise(rna_tpm = sum(rna_tpm), .groups = "drop")

clr_family <- function(long_df, value_col, out_col) {
    wide <- long_df %>%
        pivot_wider(names_from = Family, values_from = !!sym(value_col),
            values_fill = 0)
    mat <- wide %>% column_to_rownames("Sample") %>% as.matrix()
    clr_mat <- compositions::clr(mat)
    as.data.frame(clr_mat) %>%
        rownames_to_column("Sample") %>%
        pivot_longer(cols = -Sample, names_to = "Family", values_to = out_col)
}

dna_clr_family <- clr_family(dna_family, value_col = "dna_tpm",
    out_col = "clr_dna")
rna_clr_family <- clr_family(rna_family, value_col = "rna_tpm",
    out_col = "clr_rna")

clr_family_data <- dna_clr_family %>%
    left_join(rna_clr_family, by = c("Sample", "Family")) %>%
    left_join(sample_meta, by = "Sample") %>%
    filter(!is.na(Fraction))
