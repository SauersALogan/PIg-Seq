# =============================================================================
# Script:       gene_sets.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Builds functional gene set modules from Bakta annotations
#               and KEGG pathway search terms. Exports 
#               functional_modules_filtered for use in GSEA analyses.
# Dependencies: helpers.R must be sourced first
# Input:        parsed_gff.txt (https://doi.org/10.48420/33951214)
#               KEGG_pathway_search_terms.xlsx (https://doi.org/10.48420/33951214)
# Output:       functional_modules_filtered (list of locus tag vectors)
#               kegg_extraction_validation.csv
# =============================================================================

library(readr)
library(readxl)
library(dplyr)
library(tibble)

source("statistics/functions/helpers.R")

## ---- Paths (update for your system) ----
data_dir <- "/path/to/data"

## ---- Load data ----
Annotations <- read_delim(file.path(data_dir, "parsed_gff.txt"), delim = ",",
    col_types = cols(Locus_Tag = col_character(), Gene = col_character(),
    Product = col_character(), BlastRules = col_character(),
    COG = col_character(), EC = col_character(), GO = col_character(),
    IS = col_character(), KEGG = col_character(), NCBIFam = col_character(),
    NCBIProtein = col_character(), PFAM = col_character(),
    RFAM = col_character(), RefSeq = col_character(),
    SO = col_character(), UniParc = col_character(), 
    UniRef = col_character(), VFDB = col_character()))

kegg_pathway_search_terms <- read_excel(file.path(data_dir, "KEGG_pathway_search_terms.xlsx"))

## ---- Module patterns (regex-based) ----
module_patterns <- list(
    Dissimilatory_sulfite_reduction = list(include = paste(
        "dissimilatory.*sulfite reductase|sulfite reductase.*dissimilatory",
        "\\bdsrA\\b|\\bdsrB\\b|\\bdsrC\\b|\\bdsrD\\b|\\bdsrE\\b|\\bdsrF\\b",
        "\\bdsrM\\b|\\bdsrK\\b|\\bdsrJ\\b|\\bdsrO\\b|\\bdsrP\\b",
        "adenylyl-sulfate reductase|adenosine-5-phosphosulfate reductase",
        "sulfate adenylyltransferase|\\bsat\\b.*sulfate",
        sep = "|"),
        exclude = "assimilatory|domain-containing|putative"),
  
    Sporulation = "\\bspo0[A-Z]|\\bspo[II]+[A-Z]|\\bspo[IV]+[A-Z]|
        \\bspo[V]+[A-Z]|\\bsporulation\\b|\\bsigE\\b.*sporulation|
        \\bsigF\\b.*sporulation|\\bsigG\\b.*sporulation|
        \\bsigK\\b.*sporulation|spore coat|spore cortex")

## ---- KEGG pathway search terms ----
kegg_pathway_search_terms <- read_excel("KEGG_pathway_search_terms.xlsx")

## ---- Extract modules ----
kegg_matches <- extract_kegg_modules(
    kegg_pathway_data = kegg_pathway_search_terms,
    user_annotations = Annotations,
    product_col = "Product",
    kegg_col = "KEGG")

kegg_functional_modules <- kegg_matches %>%
    dplyr::group_by(Pathway) %>%
    dplyr::summarise(Locus_Tags = list(unique(Locus_Tag)), .groups = "drop") %>%
    {setNames(.$Locus_Tags, .$Pathway)}

pattern_modules <- extract_pattern_modules(
    module_patterns = module_patterns,
    user_annotations = Annotations,
    product_col = "Product")

## ---- Merge and filter ----
functional_modules <- c(kegg_functional_modules, pattern_modules)

min_genes <- 5
max_genes <- 5000
functional_modules_filtered <- functional_modules[
    sapply(functional_modules, length) >= min_genes & 
    sapply(functional_modules, length) <= max_genes
]

cat("\nTotal modules:", length(functional_modules_filtered), "\n")
cat("Total genes (non-unique):", 
    sum(sapply(functional_modules_filtered, length)), "\n")

## ---- Validate KEGG extraction ----
validation_results <- validate_kegg_extraction(
    kegg_pathway_data = kegg_pathway_search_terms)

write.csv(validation_results, "kegg_extraction_validation.csv", 
    row.names = FALSE)

## ---- Test module overlaps ----
overlap_df <- make_overlap_df(functional_modules_filtered)
