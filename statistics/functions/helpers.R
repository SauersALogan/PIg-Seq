# =============================================================================
# Script:       helpers.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Shared analysis functions sourced by all analysis notebooks.
#               Includes parallelized lmer, GSEA, module extraction,
#               validation, and plotting utilities.
# =============================================================================

library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(lme4)
library(emmeans)
library(MuMIn)
library(fgsea)
library(parallel)
library(future)
library(future.apply)
library(progressr)
library(ggplot2)
library(ggrastr)

## ---- GSEA helper functions ----
run_fgsea_perbin <- function(bin_id, DEGs, selected_contrast, gene_set_data,
    min_genes = 3, min_set_size = 3, max_set_size = 4000) {
  
    DEGs_bin <- DEGs %>%
        dplyr::filter(Bin == bin_id, contrast == selected_contrast)
  
    if (nrow(DEGs_bin) < min_genes) return(NULL)
  
    gene_ranks <- setNames(DEGs_bin$clr_t, DEGs_bin$Geneid)
    gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  
    bin_prefix <- paste0("^", bin_id, "\\.")
    pathway_list <- lapply(gene_set_data, function(gene_list) {
        gene_list[grepl(bin_prefix, gene_list)]
    })
  
    pathway_list <- pathway_list[sapply(pathway_list, length) >= min_set_size]
  
    if (length(pathway_list) == 0) return(NULL)
  
    fgsea_results <- fgsea(pathways = pathway_list,
        stats = gene_ranks, minSize = min_set_size,
        maxSize = max_set_size, nproc = 1
    )
  
    fgsea_results$Bin <- bin_id
    return(as.data.frame(fgsea_results))
}

run_fgsea_parallel <- function(bin_list, DEGs, selected_contrast, gene_set_data,
    analysis_name = "GSEA", min_genes = 100, min_set_size = 5, 
    max_set_size = 500, n_cores = parallel::detectCores() - 1) {
  
    cat("\n=== Starting parallelized", analysis_name, "===\n")
    cat("Bins to analyze:", length(bin_list), "\n")
    cat("Cores:", n_cores, "\n\n")
  
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("DEGs", "selected_contrast", "gene_set_data",
        "min_genes", "min_set_size", "max_set_size", "run_fgsea_perbin"),
        envir = environment())
    
    clusterEvalQ(cl, { 
        library(dplyr); library(fgsea); library(tidyr); library(tibble) 
    })
  
    all_results <- parLapply(cl, bin_list, function(bin_id) {
        run_fgsea_perbin(bin_id, DEGs, selected_contrast, gene_set_data,
            min_genes, min_set_size, max_set_size)
    })
  
    stopCluster(cl)
  
    gsea_results <- dplyr::bind_rows(all_results)
  
    cat("\n===", analysis_name, "Complete ===\n")
    cat("Total results:", nrow(gsea_results), "\n")
    cat("Unique bins:", length(unique(gsea_results$Bin)), "\n")
    cat("Unique pathways/modules:", length(unique(gsea_results$pathway)), "\n")
    cat("Significant (padj < 0.05):", sum(gsea_results$padj < 0.05, 
        na.rm = TRUE), "\n")
    cat("Significant (padj < 0.25):", sum(gsea_results$padj < 0.25, 
        na.rm = TRUE), "\n\n")
  
    return(gsea_results)
}

## ---- Module extraction helper function ----
extract_pattern_modules <- function(module_patterns, 
    user_annotations, product_col = "Product") {
     
    cat("\n=== Extracting pattern-based gene sets ===\n")

    standard_exclusions <- paste("domain-containing protein",
        "domain protein$", "-type domain", "family protein$", "C4-type",
        "-like protein", "putative", "hypothetical", "\\bfamily\\b",
        sep = "|")
  
    functional_modules <- list()
  
    for(module_name in names(module_patterns)) {
    
        cat("Processing:", module_name, "... ")
    
        pattern_spec <- module_patterns[[module_name]]
    
        if(is.character(pattern_spec)) {
            genes <- user_annotations %>%
                dplyr::filter(grepl(pattern_spec, !!sym(product_col), 
            ignore.case = TRUE, perl = TRUE)) %>%
                dplyr::filter(!grepl(standard_exclusions, !!sym(product_col), 
            ignore.case = TRUE, perl = TRUE)) %>%
                dplyr::pull(Locus_Tag) %>%
            unique()
      
            functional_modules[[module_name]] <- genes
            cat(length(genes), "genes\n")
            next
        }
    
        if(is.list(pattern_spec)) {
      
            if("include" %in% names(pattern_spec)) {
        
                genes <- user_annotations %>%
                dplyr::filter(grepl(pattern_spec$include, !!sym(product_col), 
                    ignore.case = TRUE, perl = TRUE))
        
                if("exclude" %in% names(pattern_spec)) {
                    genes <- genes %>%
                        dplyr::filter(!grepl(pattern_spec$exclude, 
                            !!sym(product_col), 
                                ignore.case = TRUE, perl = TRUE))
                }
        
            genes <- genes %>%
                dplyr::filter(!grepl(standard_exclusions, !!sym(product_col), 
                    ignore.case = TRUE, perl = TRUE))
        
            functional_modules[[module_name]] <- unique(genes$Locus_Tag)
            cat(length(unique(genes$Locus_Tag)), "genes\n")
            next
            }
        }
    
        cat("WARNING: Unrecognized pattern format\n")
        functional_modules[[module_name]] <- character(0)
    }
  
    cat("\n")
    return(functional_modules)
}

extract_kegg_modules <- function(kegg_pathway_data,
    user_annotations, product_col = "Product", kegg_col = "KEGG") {
  
    exclude_from_name_matching <- c("regulatory protein", 
    "outer membrane protein", "immunity protein", "hemolysin", "flagellin", 
    "acyltransferase", "alcohol dehydrogenasa e", "aldehyde dehydrogenase", 
    "glucosyltransferase", "serine protease", "diguanylate cyclase", 
    "alkaline phosphatase", "cysteine desulfurase", "adenylate cyclase", 
    "c-di-GMP phosphodiesterase", "kinase", "transferase", "starch synthase")
  
    cat("Excluded from name matching:", 
        length(exclude_from_name_matching), "generic terms\n")
  
    pathway_ko_exclusions <- list("Biofilm" = c("K01657", "K01658", "K01912"),
        "Valine biosynthesis" = c("K01649"), 
        "Oxidative phosphorylation" = c("K02319"))
  
    n_pathway_exclusions <- sum(sapply(pathway_ko_exclusions, length))
    cat("Pathway-specific KO exclusions:", n_pathway_exclusions, "\n")
  
    standard_exclusion_pattern <- paste("domain-containing",
        "domain protein$", "-type domain", "family protein$", "C4-type",
        "-like protein", "putative", "hypothetical", "\\bfamily\\b", sep = "|")
  
    specific_exclusions <- list(
        list(ko_numbers = "K01768", pattern = "diadenylate|di-adenylate|DisA"),
        list(ko_numbers = "K00640", pattern = "homoserine"),
        list(ko_numbers = "K18101", pattern = "RNA helicase"),
        list(ko_numbers = "K02659", pattern = "PilT|assembly chaperone"),
        list(ko_numbers = "K10924", pattern = "MshJ|MshL"),
        list(ko_numbers = "K06204", pattern = "CarD"),
        list(ko_numbers = "K03666", pattern = "RNA$|binding RNA"),
        list(ko_numbers = "K10932", pattern = "Conjugative transposon"),
        list(ko_numbers = "K03557", pattern = "PAS.*ATPase"),
        list(ko_numbers = "K03469", pattern = "\\bHII\\b|\\bHIII\\b"),
        list(ko_numbers = "K03470", pattern = "\\bHIII\\b"),
        list(ko_numbers = "K02340", pattern = "delta'|delta prime"),
        list(ko_numbers = "K02316", pattern = "TraC|conjugative"),
        list(ko_numbers = "K02314", pattern = "loader|DnaC"),
        list(ko_numbers = "K22391", 
            pattern = "cyclohydrolase II|cyclohydrolase 2"),
        list(ko_numbers = "K00011", pattern = "phosphonoacetaldehyde"),
        list(ko_numbers = "K03752", 
            pattern = "conjugal transfer|conjugative|plasmid"),
        list(ko_numbers = "K13507", 
            pattern = "tRNA|aminotransferase|amidotransferase"),
        list(ko_numbers = c("K00937", "K22468"), 
            pattern = "RNA degradosome|RNA processing"),
        list(ko_numbers = "K00343", pattern = "Nuo[A-MO-Z]"),
        list(ko_numbers = "K08738", 
            pattern = paste("biogenesis|biosynthesis|assembly|maturation",
                "oxidase|reductase|dehydrogenase",
                "nitrite|nitrate|peroxidase|TMAO",
                "CcsA|CcsB|CcmA|CcmE|CcmF|CcmH|CcdA|CcoG|CcoQ|ResB|CycH",
                "NrfA|NrfH|NrfG|NapC|QrcA|DsrJ", sep = "|")),
        list(ko_numbers = "K01659", 
            pattern = paste(
                "[Pp]rotein phosphatase|[Ss]erine/threonine|[Pp]hosphatase",
                sep = "|")),
        list(ko_numbers = "K01734", pattern = "[Rr]ecombination|DNA repair"),
        list(ko_numbers = "K18471", pattern = "[Uu]ncharacterized"),
        list(ko_numbers = "K00016", 
            pattern = "[Rr]egulatory|[Rr]egulator|operon regulatory"),
        list(ko_numbers = c("K00625", "K13788"), 
            pattern = "glucosamine|GlmU|uridyltransferase"),
        list(ko_numbers = "K00248", pattern = "3-hydroxybutyryl"),
        list(ko_numbers = c("K00024", "K00025", "K00026"), 
            pattern = paste(
                "isopropylmalate|isopropyl-malate|leucine", sep="|")),
        list(ko_numbers = "K01638", 
            pattern = paste(
                "isopropylmalate|citramalate|glucosaminyl|leucine|isoleucine",
                sep = "|")),
        list(ko_numbers = "K01958", 
            pattern = "[Pp]hosphoenolpyruvate carboxylase|PEP carboxylase"),
        list(ko_numbers = "K01573", pattern = "transporter|related small"),
        list(ko_numbers = "K03076", pattern = "interacting protein|Syd"),
        list(ko_numbers = "K12295", 
            pattern = "ComE operon protein [0-9]|operon protein"),
        list(ko_numbers = "K01728", 
            pattern = "Pel transporter|PelG|exopolysaccharide Pel"),
        list(ko_numbers = "K20272", 
            pattern = "YacL|TRAM-domain|recognition site|Deoxyribonuclease"),
        list(ko_numbers = "K20380", pattern = "relaxase"),
        list(ko_numbers = "K02250", pattern = "Control of|YlbF|YmcA"),
        list(ko_numbers = "K01998", pattern = "LivH"),
        list(ko_numbers = "K15581", pattern = "OppC"),
        list(ko_numbers = "K15582", pattern = "OppB"),
        list(ko_numbers = "K01078", 
            pattern = paste(
                "[Pp]hosphatidic acid|lipid|[Hh]istidine acid",
                "[Pp]urple acid|phosphate-irrepressible", sep ="|")),
        list(ko_numbers = "K03148", pattern = "UBA/THIF-type(?! thiamine)|generic"),
        list(ko_numbers = "K03153", 
             pattern = "[Dd]isulfide|thioredoxin|Thio:disulfide"),
        list(ko_numbers = NA,
             pathway = "Phage tail",
             pattern = "adhesin|collagen-binding|shufflon|PilV|capsid assembly"),
        list(ko_numbers = NA,
             pathway = "Phage lysis",
             pattern = paste("choline|betaine|carnitine|nicotinamid|pyrazinamid",
                 "omega-amidase|formamidase|acetamidase|carbamoylputrescine",
                 "mycothiol|glucosylceramidase|CheD|deamidase",
                 "ribosomal.*protease|tellurite.*protease|YraA|Prp",
                 "phosphatidylcholine|phosphocholine|cytidylyltransferase",
                 sep = "|")))
  
    cat("Specific false positive patterns:", 
    length(specific_exclusions), "\n\n")
  
    all_matches <- data.frame()
  
    for(i in seq_len(nrow(kegg_pathway_data))) {
    
        ko <- kegg_pathway_data$ko_number[i]
        abbr <- kegg_pathway_data$abbreviation[i]
        name <- kegg_pathway_data$name[i]
        pathway <- kegg_pathway_data$pathway[i]
    
        if(pathway %in% names(pathway_ko_exclusions)) {
            if(ko %in% pathway_ko_exclusions[[pathway]]) {
            next
            }
        }
    
        matches <- data.frame()
    
        if(!is.na(abbr) && abbr != "" && abbr != "NA") {
            abbr_pattern <- paste0("\\b", abbr, "\\b")
            abbr_matches <- user_annotations %>%
            dplyr::filter(grepl(abbr_pattern, !!sym(product_col), 
            ignore.case = TRUE, perl = TRUE)) %>%
            dplyr::mutate(Search_method = "Abbreviation")
      
            if(nrow(abbr_matches) > 0) {
                matches <- rbind(matches, abbr_matches)
            }
        }
    
        if(!is.na(name) && name != "" && !(name %in% exclude_from_name_matching)) {
            if(name == "DNA polymerase I") {
                name_matches <- user_annotations %>%
                dplyr::filter(grepl("DNA polymerase I", !!sym(product_col), 
                    ignore.case = TRUE)) %>%
                dplyr::filter(!grepl("DNA polymerase I(I|V)", 
                    !!sym(product_col), ignore.case = TRUE)) %>%
                dplyr::mutate(Search_method = "Name")
            } else {
                name_escaped <- gsub("([\\[\\]\\(\\)\\{\\}\\+\\*\\?\\.])", 
                "\\\\\\1", name)
                name_matches <- user_annotations %>%
                dplyr::filter(grepl(name_escaped, !!sym(product_col), 
                    ignore.case = TRUE, perl = TRUE)) %>%
                dplyr::mutate(Search_method = "Name")
            }
      
            if(nrow(name_matches) > 0) {
                matches <- rbind(matches, name_matches)
            }
        }
    
        if(nrow(matches) > 0) {
            matches <- matches %>%
            dplyr::group_by(Locus_Tag) %>%
            dplyr::mutate(Search_method = paste(unique(Search_method), 
                collapse = "+")) %>%
            dplyr::ungroup() %>%
            dplyr::distinct(Locus_Tag, .keep_all = TRUE)
      
            matches <- matches %>%
            dplyr::filter(!grepl(standard_exclusion_pattern, 
                !!sym(product_col), ignore.case = TRUE, perl = TRUE))
      
            for(exclusion in specific_exclusions) {
                if(!is.null(exclusion$pathway)) {
                    if(pathway == exclusion$pathway) {
                        matches <- matches %>%
                        dplyr::filter(!grepl(exclusion$pattern, 
                            !!sym(product_col), ignore.case = TRUE, perl = TRUE))
                    }
                } else {
                    if(ko %in% exclusion$ko_numbers) {
                        matches <- matches %>%
                        dplyr::filter(!grepl(exclusion$pattern, 
                            !!sym(product_col), ignore.case = TRUE, perl = TRUE))
                    }
                }
            }
      
            if(nrow(matches) > 0) {
                matches <- matches %>%
                dplyr::mutate(Pathway = pathway,
                    KO_Number = ko,
                    Abbreviation = ifelse(is.na(abbr), "", abbr),
                    KEGG_Name = ifelse(is.na(name), "", name)
                )
        
                all_matches <- rbind(all_matches, matches)
            }
        }
    }
  
    cat("\nTotal genes matched:", nrow(all_matches), "\n\n")
  
    return(all_matches)
}

## ---- Custom module validation function ----
validate_kegg_extraction <- function(kegg_pathway_data,
    user_product_col = "name", user_kegg_col = "ko_number") {

    cat("\n VALIDATING KEGG EXTRACTION \n")

    reference_as_annotations <- kegg_pathway_data %>%
    dplyr::filter(!is.na(name), name != "") %>%
    dplyr::rename(Locus_Tag = ko_number, Product = name)
        
    kegg_matches_validation <- extract_kegg_modules(
        kegg_pathway_data = kegg_pathway_data,
        user_annotations  = reference_as_annotations,
        product_col = "Product", kegg_col = user_kegg_col)

    results_summary <- data.frame()

    for(pathway_name in unique(kegg_pathway_data$pathway)) {

        cat("\n PATHWAY:", pathway_name, "\n")

        reference_products <- kegg_pathway_data %>%
            dplyr::filter(pathway == pathway_name,
                !is.na(name), name != "") %>%
            dplyr::select(ko_number, name) %>%
            dplyr::rename(Product = name)

        n_reference_total <- nrow(reference_products)
        cat("Reference genes (from KEGG names):", n_reference_total, "\n")

        if(n_reference_total == 0) {
            cat("No gene names in kegg_pathway_data for this pathway\n")
            next
        }

        pathway_matches <- kegg_matches_validation %>%
            dplyr::filter(Pathway == pathway_name)

        reference_products <- reference_products %>%
            dplyr::mutate(
                captured = ko_number %in% pathway_matches$KO_Number,
                search_method = dplyr::case_when(
                    ko_number %in% pathway_matches$KO_Number[
                    pathway_matches$Search_method == "Abbreviation"] ~ "Abbreviation",
                    ko_number %in% pathway_matches$KO_Number[
                    pathway_matches$Search_method == "Name"] ~ "Name",
                    ko_number %in% pathway_matches$KO_Number[
                    pathway_matches$Search_method == "Abbreviation+Name"] ~ "Both",
                    TRUE ~ NA_character_))

        n_captured <- sum(reference_products$captured)
        n_missed <- n_reference_total - n_captured
        pct_captured <- round(100 * n_captured / n_reference_total, 1)
        n_abbr_only <- sum(reference_products$search_method == "Abbreviation", 
            na.rm = TRUE)
        n_name_only <- sum(reference_products$search_method == "Name", 
            na.rm = TRUE)
        n_both <- sum(reference_products$search_method == "Both", 
            na.rm = TRUE)

        cat(sprintf("\nCAPTURED: %d / %d (%.1f%%)\n", 
            n_captured, n_reference_total, pct_captured))
        cat(sprintf("  - Abbreviation only: %d\n", n_abbr_only))
        cat(sprintf("  - Name only: %d\n", n_name_only))
        cat(sprintf("  - Both: %d\n", n_both))
        cat(sprintf("\nMISSED: %d (%.1f%%)\n", n_missed, 100 - pct_captured))

        if(n_missed > 0) {
            cat("\nMissed genes:\n")
            missed_genes <- reference_products %>% dplyr::filter(!captured)
            for(i in 1:min(30, nrow(missed_genes))) {
                cat(sprintf("  %2d. [%s] %s\n", i,
                    missed_genes$ko_number[i],
                    missed_genes$Product[i]))
            }
            if(n_missed > 30) cat(sprintf("  ... and %d more\n", n_missed - 30))
        }

        results_summary <- rbind(results_summary, data.frame(
            Pathway = pathway_name, Reference_genes = n_reference_total,
            Captured = n_captured, Missed = n_missed,
            Pct_captured = pct_captured, Abbr_only = n_abbr_only,
            Name_only = n_name_only, Both = n_both, stringsAsFactors = FALSE))
    }

    cat("\nVALIDATION SUMMARY\n")
    print(results_summary)

    cat("\nOverall statistics:\n")
    cat(sprintf("Total reference genes: %d\n", 
        sum(results_summary$Reference_genes)))
    cat(sprintf("Total captured: %d (%.1f%%)\n",
        sum(results_summary$Captured),
        100 * sum(results_summary$Captured) /
        sum(results_summary$Reference_genes)))
    cat(sprintf("Total missed: %d (%.1f%%)\n",
        sum(results_summary$Missed),
        100 * sum(results_summary$Missed) /
        sum(results_summary$Reference_genes)))

    return(results_summary)
}

## ---- Differential expression function (can use for abundance also) ----
run_lmer_parallel <- function(data_list, formula, emmeans_term,
    fdr_group, singular_action = c("null", "flag_keep", "flag_exclude"),
    output_cols = list(), inner_split = NULL) {
    
    singular_action <- match.arg(singular_action)
    
    cat("\n=== Running parallel lmer models ===\n")
    cat("Items to model:", length(data_list), "\n")
    cat("Formula:", deparse(formula), "\n")
    cat("Singular action:", singular_action, "\n")
    
    start_time <- Sys.time()
    
    with_progress({
        p <- progressor(steps = length(data_list))
        
        results <- future_lapply(names(data_list), function(item_name) {
            d <- data_list[[item_name]]
            
            sub_list <- if (!is.null(inner_split)){
                split(d, d[[inner_split]])
            } else { 
               list(d) 
            }
            
            sub_results <- lapply(sub_list, function(d) {
            
                tryCatch({
                    m <- lmer(formula, data = d,
                        control = lmerControl(optimizer = "bobyqa"))
                
                    singular_flag <- isSingular(m)
                
                    if (singular_flag && singular_action == "null")           
                        return(NULL)
                    if (singular_flag && singular_action == "flag_exclude")   
                        return(NULL)
                
                    r2_vals <- r.squaredGLMM(m)
                    residuals_m <- residuals(m)
                    shapiro_p <- if (length(residuals_m) >= 3 &&
                    length(residuals_m) <= 5000) {
                        shapiro.test(residuals_m)$p.value
                    } else NA
                
                    em   <- emmeans(m, emmeans_term)
                    cf   <- contrast(em, method = "pairwise", adjust = "none")
                    summ <- as.data.frame(summary(cf))
                
                    summ <- summ %>%
                    dplyr::mutate(r_squared = r2_vals[1, "R2m"],
                        adj_r_squared = r2_vals[1, "R2c"],shapiro_p = shapiro_p,
                        singular_fit = singular_flag,n_obs = nrow(d))
                
                
                    for (col_name in names(output_cols)) {
                        summ[[col_name]] <- eval(output_cols[[col_name]], 
                            envir = d)
                    }
                
                    summ
                
                }, error = function(e) {
                    message("Error in ", item_name, ": ", e$message)
                    return(NULL)
                })
            })
            
        p()
        return(dplyr::bind_rows(sub_results))
            
        }, future.seed = TRUE)
    })
    
    results_df <- dplyr::bind_rows(results)
    
    cat("Time:", round(difftime(Sys.time(), 
        start_time, units = "mins"), 2), "minutes\n")
    cat("Items with results:", 
        n_distinct(results_df[[names(output_cols)[1]]]), "\n")
    
    results_df <- results_df %>%
    dplyr::group_by(across(all_of(fdr_group))) %>%
    dplyr::mutate(p_adjust = p.adjust(p.value, method = "BH"),
        significant = p_adjust < 0.05) %>%
    dplyr::ungroup()
    
    cat("Significant results (FDR < 0.05):", 
        sum(results_df$significant, na.rm = TRUE), "\n")
    
    return(results_df)
}

## ---- Family volcano plots ----
plot_family_product_volcano <- function(data, family_name) {

    fam_data <- data %>% dplyr::filter(Family == family_name)

    cat(family_name, "- total products:", nrow(fam_data),
        "| significant:", sum(fam_data$significant, na.rm = TRUE), "\n")

    ggplot(fam_data, aes(x = estimate, y = neg_log10_p)) +
        geom_point_rast(aes(color = direction),
            alpha = 0.5, size = pt_small, raster.dpi = 600) +
        geom_vline(xintercept = 0, linetype = "dashed",
            color = "grey50", linewidth = lw_sm) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed",
            color = "grey50", linewidth = lw_sm) +
        scale_color_manual(name = NULL,
            values = c("Up in Positive"   = col_pos,
                       "Down in Positive" = col_neg,
                       "Not significant"  = "grey70")) +
        scale_x_continuous(breaks = c(-8, -6, -4, -2, 0, 2, 4, 6), 
                          limits = c(-8,6)) +
        theme_for_figures() +
        theme(plot.margin = ggplot2::margin(0, 0, 0, 0)) +
        labs(x = "Δ CLR(RNA)",
             y = "-log10(FDR)")
}

## ---- Custom module overlap function ----
make_overlap_df <- function(functional_modules) {
    module_names <- names(functional_modules)
    purrr::map_dfr(module_names, function(m1) {
        purrr::map_dfr(module_names, function(m2) {
            if(m1 >= m2) return(NULL)
            shared <- length(intersect(functional_modules[[m1]], 
                                       functional_modules[[m2]]))
            tibble::tibble(
                Module1       = m1,
                Module2       = m2,
                n_shared      = shared,
                pct_of_m1     = round(100 * shared / length(functional_modules[[m1]]), 1),
                pct_of_m2     = round(100 * shared / length(functional_modules[[m2]]), 1),
                jaccard       = round(100 * shared / length(union(
                    functional_modules[[m1]], functional_modules[[m2]])), 1)
            )
        })
    })
}
