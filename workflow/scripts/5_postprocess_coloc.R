#!/usr/bin/env Rscript

#' Post-processing Module for Coloc Results
#'
#' Aggregates and filters colocalization results
#' Refactored from src/5_postprocess_coloc.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(optparse)
})

#' Aggregate coloc results across tissues and traits
#'
#' @param results_dir Directory containing coloc result files
#' @param pp4_threshold Numeric: PP4 threshold for significance
#' @return data.table with aggregated results
aggregate_coloc_results <- function(results_dir, pp4_threshold = 0.8) {
  log_message("Aggregating coloc results...")
  
  result_files <- list.files(results_dir, pattern = "colocABF.*\\.txt$", 
                             full.names = TRUE, recursive = TRUE)
  
  all_results <- rbindlist(lapply(result_files, fread), fill = TRUE)
  
  # Filter significant results
  significant <- all_results[PP4 >= pp4_threshold | 
                             (PP4 > 0.6 & PP4 < 0.8 & PP4/PP3 > 2)]
  
  log_message(sprintf("Found %d significant colocalizations", nrow(significant)))
  
  return(significant)
}

#' Create summary statistics
#'
#' @param results data.table with coloc results
#' @return data.table with summary by trait and tissue
summarize_results <- function(results) {
  summary <- results %>%
    group_by(GWAS_ID, tissue) %>%
    summarise(
      n_colocs = n(),
      n_genes = n_distinct(cpg),
      mean_PP4 = mean(PP4),
      .groups = "drop"
    ) %>%
    arrange(desc(n_colocs))
  
  return(as.data.table(summary))
}

# CLI interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-i", "--input"), type = "character",
                help = "Input directory with coloc results"),
    make_option(c("-o", "--output"), type = "character",
                help = "Output file for aggregated results"),
    make_option(c("--pp4"), type = "double", default = 0.8,
                help = "PP4 threshold for significance")
  )
  
  opt <- parse_args(OptionParser(option_list = option_list))
  
  results <- aggregate_coloc_results(opt$input, opt$pp4)
  summary <- summarize_results(results)
  
  fwrite(results, opt$output)
  fwrite(summary, gsub("\\.txt$", "_summary.txt", opt$output))
  
  log_message("Post-processing complete")
}

log_message <- function(msg, level = "INFO") {
  cat(sprintf("[%s] %s: %s\n", Sys.time(), level, msg))
}
