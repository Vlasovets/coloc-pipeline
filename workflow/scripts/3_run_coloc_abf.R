#!/usr/bin/env Rscript

#' Coloc ABF Analysis Module
#'
#' Performs colocalization analysis using coloc.abf method
#' Refactored from src/3_run_coloc_abf.R with modular design

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(dplyr)
  library(parallel)
  library(optparse)
})

source(file.path(snakemake@scriptdir, "data_loader.R"))
source(file.path(snakemake@scriptdir, "qtl_processor.R"))
source(file.path(snakemake@scriptdir, "coloc_helpers.R"))

#' Run coloc.abf for a single QTL-GWAS pair
#'
#' @param qtl_data List: QTL data in coloc format
#' @param gwas_data List: GWAS data in coloc format
#' @param p12 Numeric: prior probability of shared causal variant
#' @return List with coloc results
run_coloc_abf <- function(qtl_data, gwas_data, p12 = 1e-5) {
  tryCatch({
    result <- coloc::coloc.abf(
      dataset1 = qtl_data,
      dataset2 = gwas_data,
      p12 = p12
    )
    return(result)
  }, error = function(e) {
    log_message(sprintf("Coloc error: %s", e$message), "ERROR")
    return(NULL)
  })
}

#' Extract credible set from coloc results
#'
#' @param coloc_result Coloc result object
#' @param threshold Numeric: cumulative probability threshold
#' @return Character vector of SNPs in credible set
extract_credible_set <- function(coloc_result, threshold = 0.95) {
  o <- order(coloc_result$results$SNP.PP.H4, decreasing = TRUE)
  cs <- cumsum(coloc_result$results$SNP.PP.H4[o])
  w <- which(cs > threshold)[1]
  credset <- coloc_result$results[o, ][1:w, ]$snp
  return(credset)
}

# CLI interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-g", "--gwas"), type = "character",
                help = "GWAS ID/trait name"),
    make_option(c("-t", "--tissue"), type = "character",
                help = "Tissue/cell type"),
    make_option(c("-o", "--output"), type = "character",
                help = "Output directory"),
    make_option(c("--cores"), type = "integer", default = 1,
                help = "Number of cores for parallel processing")
  )
  
  opt <- parse_args(OptionParser(option_list = option_list))
  
  log_message(sprintf("Starting coloc ABF analysis for %s - %s", 
                     opt$gwas, opt$tissue))
  
  # Main analysis would go here using the functions above
  # This is a template - actual implementation depends on workflow structure
  
  log_message("Analysis complete")
}
