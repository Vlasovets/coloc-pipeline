#!/usr/bin/env Rscript

#' Coloc SuSiE Analysis Module
#'
#' Performs colocalization with SuSiE for multiple causal variants
#' Refactored from src/4_run_coloc_susie.R

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(dplyr)
  library(optparse)
})

source("workflow/scripts/coloc_helpers.R")

#' Run coloc.susie analysis
#'
#' @param qtl_data List: QTL data with LD matrix
#' @param gwas_data List: GWAS data with LD matrix
#' @return List with SuSiE coloc results
run_coloc_susie <- function(qtl_data, gwas_data) {
  tryCatch({
    result <- coloc.susie(qtl_data, gwas_data)
    return(result)
  }, error = function(e) {
    log_message(sprintf("Coloc SuSiE error: %s", e$message), "ERROR")
    return(NULL)
  })
}

# CLI interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-g", "--gwas"), type = "character"),
    make_option(c("-t", "--tissue"), type = "character"),
    make_option(c("-o", "--output"), type = "character")
  )
  
  opt <- parse_args(OptionParser(option_list = option_list))
  log_message(sprintf("Starting coloc SuSiE for %s - %s", opt$gwas, opt$tissue))
}
