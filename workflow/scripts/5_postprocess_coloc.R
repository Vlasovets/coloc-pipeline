#!/usr/bin/env Rscript

#' Post-processing Module for Coloc Results
#'
#' Aggregates and filters colocalization results
#' Adapted from src/5_postprocess_coloc.R for Snakemake workflow

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Get parameters from Snakemake
input_files <- snakemake@input
output_all <- snakemake@output$all_results
output_sig <- snakemake@output$sig_results
pp4_threshold <- snakemake@params$pp4_threshold
pp4_pp3_ratio <- snakemake@params$pp4_pp3_ratio

cat(sprintf("[%s] Starting aggregation of coloc results\n", Sys.time()))

#############################################################################
# 1. Read all coloc ABF results
#############################################################################
cat(sprintf("[%s] Reading %d coloc result files\n", Sys.time(), length(input_files)))

coloc_res <- rbindlist(lapply(input_files, function(f) {
  if(file.exists(f) && file.size(f) > 0) {
    dt <- fread(f)
    # Extract tissue and GWAS_ID from filename if not in data
    if(!"tissue" %in% colnames(dt)) {
      tissue <- gsub(".*colocABF_df_([^_]+)_.*", "\\1", basename(f))
      dt$tissue <- tissue
    }
    if(!"GWAS_ID" %in% colnames(dt)) {
      gwas_id <- gsub(".*results_([^/]+)/.*", "\\1", f)
      dt$GWAS_ID <- gwas_id
    }
    return(dt)
  } else {
    return(data.table())
  }
}), fill = TRUE)

if(nrow(coloc_res) == 0) {
  cat(sprintf("[%s] Warning: No coloc results found\n", Sys.time()))
  # Create empty output files
  fwrite(data.table(), output_all)
  fwrite(data.table(), output_sig)
  quit(save="no", status=0)
}

#############################################################################
# 2. Calculate summary statistics
#############################################################################
cat(sprintf("[%s] Total coloc tests: %d\n", Sys.time(), nrow(coloc_res)))
cat(sprintf("[%s] Unique GWAS signals tested: %d\n", Sys.time(), 
            length(unique(coloc_res$gwas_signal))))
cat(sprintf("[%s] Unique genes tested: %d\n", Sys.time(), 
            length(unique(coloc_res$cpg))))

# Summary by tissue and GWAS_ID
summary_dt <- coloc_res %>%
  group_by(GWAS_ID, tissue) %>%
  summarise(
    n_tests = n(),
    n_genes = n_distinct(cpg),
    n_signals = n_distinct(gwas_signal),
    mean_PP4 = mean(PP4, na.rm=TRUE),
    median_PP4 = median(PP4, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(GWAS_ID, tissue)

cat("\nSummary by tissue and trait:\n")
print(summary_dt)

#############################################################################
# 3. Save all results
#############################################################################
cat(sprintf("\n[%s] Saving all results to %s\n", Sys.time(), output_all))
fwrite(coloc_res, output_all, sep="\t")

#############################################################################
# 4. Filter for significant results
#############################################################################
cat(sprintf("[%s] Filtering for significant colocalizations\n", Sys.time()))
cat(sprintf("[%s] Criteria: PP4 >= %.2f OR (PP4 > 0.6 AND PP4/PP3 > %.1f)\n", 
            Sys.time(), pp4_threshold, pp4_pp3_ratio))

coloc_res_sign <- coloc_res[
  PP4 >= pp4_threshold | 
  (PP4 > 0.6 & PP4 < pp4_threshold & PP4/PP3 > pp4_pp3_ratio)
]

cat(sprintf("[%s] Significant colocalizations: %d\n", Sys.time(), nrow(coloc_res_sign)))

if(nrow(coloc_res_sign) > 0) {
  cat(sprintf("[%s] Unique significant GWAS signals: %d\n", Sys.time(), 
              length(unique(coloc_res_sign$gwas_signal))))
  cat(sprintf("[%s] Unique significant genes: %d\n", Sys.time(), 
              length(unique(coloc_res_sign$cpg))))
  
  # Summary of significant results by tissue
  sig_summary <- coloc_res_sign %>%
    group_by(GWAS_ID, tissue) %>%
    summarise(
      n_colocs = n(),
      n_genes = n_distinct(cpg),
      mean_PP4 = mean(PP4, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_colocs))
  
  cat("\nSignificant colocalizations by tissue and trait:\n")
  print(sig_summary)
}

#############################################################################
# 5. Save significant results
#############################################################################
cat(sprintf("\n[%s] Saving significant results to %s\n", Sys.time(), output_sig))
fwrite(coloc_res_sign, output_sig, sep="\t")

cat(sprintf("[%s] Aggregation complete\n", Sys.time()))
