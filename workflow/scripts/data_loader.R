#!/usr/bin/env Rscript

#' Data Loader Module for Colocalization Pipeline
#'
#' This module provides functions for loading and validating input data:
#' - GWAS VCF files
#' - QTL summary statistics
#' - Variant annotations
#' - Independent signals
#'
#' Refactored from src/ files with improvements:
#' - Modular design with clear function boundaries
#' - Comprehensive error handling and validation
#' - Logging for debugging and provenance
#' - CLI argument support
#' - Type checking and documentation

# =============================================================================
# DEPENDENCIES
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(VariantAnnotation)
  library(gwasvcf)
  library(GenomicRanges)
  library(dplyr)
  library(rtracklayer)
})

# =============================================================================
# LOGGING UTILITIES
# =============================================================================

#' Log message with timestamp
#'
#' @param message Character string to log
#' @param level Character: "INFO", "WARNING", "ERROR"
#' @return NULL (prints to console)
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
}

# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

#' Load and validate variant annotation file
#'
#' @param variant_ann_file Path to variant annotation file (gzipped TSV)
#' @param filter_biallelic Logical: filter to biallelic SNPs only
#' @return data.frame with variant annotations
#' @export
load_variant_annotation <- function(variant_ann_file, filter_biallelic = TRUE) {
  log_message(sprintf("Loading variant annotation from: %s", variant_ann_file))
  
  if (!file.exists(variant_ann_file)) {
    stop(sprintf("Variant annotation file not found: %s", variant_ann_file))
  }
  
  tryCatch({
    variant_ann <- fread(variant_ann_file, data.table = FALSE)
    log_message(sprintf("Loaded %d variants from annotation file", nrow(variant_ann)))
    
    # Validate required columns
    required_cols <- c("chr", "pos", "rsid", "ref", "alt")
    missing_cols <- setdiff(required_cols, colnames(variant_ann))
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    
    # Convert chromosome to integer
    variant_ann$chr <- as.integer(variant_ann$chr)
    original_count <- nrow(variant_ann)
    
    # Filter for valid chromosomes and biallelic SNPs if requested
    if (filter_biallelic) {
      variant_ann <- variant_ann[!is.na(variant_ann$chr) & 
                                  nchar(variant_ann$ref) == 1 & 
                                  nchar(variant_ann$alt) == 1, ]
      log_message(sprintf("Filtered to %d biallelic SNPs (removed %d variants)", 
                         nrow(variant_ann), original_count - nrow(variant_ann)))
    }
    
    return(variant_ann)
    
  }, error = function(e) {
    log_message(sprintf("Error loading variant annotation: %s", e$message), "ERROR")
    stop(e)
  })
}

#' Load GWAS summary statistics from VCF file
#'
#' @param vcf_file Path to GWAS VCF file (bgzipped and indexed)
#' @param regions GRanges object with regions to extract
#' @param genome_build Character: genome build ("hg19" or "hg38")
#' @return data.frame with GWAS summary statistics
#' @export
load_gwas_vcf <- function(vcf_file, regions, genome_build = "hg38") {
  log_message(sprintf("Loading GWAS data from VCF: %s", basename(vcf_file)))
  
  if (!file.exists(vcf_file)) {
    stop(sprintf("GWAS VCF file not found: %s", vcf_file))
  }
  
  if (!file.exists(paste0(vcf_file, ".tbi"))) {
    stop(sprintf("GWAS VCF index not found: %s.tbi", vcf_file))
  }
  
  tryCatch({
    param <- ScanVcfParam(which = regions)
    vcf <- readVcf(vcf_file, genome_build, param)
    
    gwas_data <- gwasvcf::vcf_to_granges(vcf) %>% 
      dplyr::as_tibble() %>%
      dplyr::select(
        chr = seqnames,
        position = start,
        rsid = ID,
        ea = ALT,
        nea = REF,
        eaf = AF,
        beta = ES,
        se = SE,
        p = LP,
        n = SS
      ) %>%
      dplyr::mutate(
        p = 10^(-p),
        chr = as.integer(as.character(chr))
      )
    
    log_message(sprintf("Extracted %d variants from GWAS VCF", nrow(gwas_data)))
    
    # Create variant ID for colocalization
    gwas_data$varid_for_coloc <- paste(gwas_data$chr, gwas_data$position, 
                                        toupper(gwas_data$ea), sep = "_")
    
    return(gwas_data)
    
  }, error = function(e) {
    log_message(sprintf("Error loading GWAS VCF: %s", e$message), "ERROR")
    stop(e)
  })
}

#' Extract GWAS summary statistics for multiple windows
#'
#' @param vcf_file Path to GWAS VCF file
#' @param signals data.frame with columns: rsid, chr, position
#' @param window_length Numeric: window size around each signal (bp)
#' @param genome_build Character: genome build
#' @return data.frame with GWAS summary statistics
#' @export
extract_gwas_windows <- function(vcf_file, signals, window_length = 1e6, 
                                  genome_build = "hg38") {
  log_message(sprintf("Extracting GWAS data for %d signals with window size %d bp", 
                     nrow(signals), window_length))
  
  # Create genomic ranges for each signal
  signals$start <- ifelse(signals$position - window_length < 0, 
                          1, signals$position - window_length)
  signals$end <- signals$position + window_length
  
  gwas_loci <- GRanges(
    seqnames = paste0("chr", signals$chr),
    ranges = IRanges(start = signals$start, end = signals$end)
  )
  
  # Extract data for all regions
  gwas_associations <- load_gwas_vcf(vcf_file, gwas_loci, genome_build)
  
  # Remove duplicates
  gwas_associations <- unique(gwas_associations)
  
  log_message(sprintf("Extracted %d unique GWAS variants", 
                     nrow(gwas_associations)))
  
  return(gwas_associations)
}

#' Load QTL permutation results (eGenes/eQTLs)
#'
#' @param perm_file Path to QTL permutation results file
#' @param qval_threshold Numeric: q-value threshold for significance
#' @return data.table with significant QTLs
#' @export
load_qtl_permutations <- function(perm_file, qval_threshold = 0.05) {
  log_message(sprintf("Loading QTL permutation results from: %s", basename(perm_file)))
  
  if (!file.exists(perm_file)) {
    stop(sprintf("QTL permutation file not found: %s", perm_file))
  }
  
  tryCatch({
    mqtl <- fread(perm_file)
    log_message(sprintf("Loaded %d QTL results", nrow(mqtl)))
    
    # Parse variant ID
    mqtl[, c("chr", "pos", "REF", "ALT", "rsid") := tstrsplit(variant_id, "_")]
    mqtl[, pos := as.integer(pos)]
    mqtl[, chr := as.integer(sub("chr", "", chr))]
    mqtl[, gene_id := phenotype_id]
    
    # Filter by q-value
    mqtl <- mqtl[qval < qval_threshold]
    log_message(sprintf("Filtered to %d significant QTLs (qval < %.2f)", 
                       nrow(mqtl), qval_threshold))
    
    return(mqtl)
    
  }, error = function(e) {
    log_message(sprintf("Error loading QTL permutations: %s", e$message), "ERROR")
    stop(e)
  })
}

#' Load QTL nominal pass results for specific genes
#'
#' @param nominal_dir Directory containing nominal pass parquet files
#' @param gene_ids Character vector of gene IDs to load
#' @param chromosomes Integer vector of chromosomes to load
#' @param file_prefix Character: file prefix pattern
#' @return data.table with QTL nominal results
#' @export
load_qtl_nominal <- function(nominal_dir, gene_ids, chromosomes, 
                              file_prefix = "cis_qtl_pairs.chr") {
  log_message(sprintf("Loading QTL nominal data for %d genes across %d chromosomes", 
                     length(gene_ids), length(chromosomes)))
  
  if (!dir.exists(nominal_dir)) {
    stop(sprintf("QTL nominal directory not found: %s", nominal_dir))
  }
  
  all_qtl <- data.table()
  
  for (chr in chromosomes) {
    parquet_file <- file.path(nominal_dir, 
                              paste0(file_prefix, chr, ".parquet"))
    
    if (!file.exists(parquet_file)) {
      log_message(sprintf("Skipping chr%d: file not found", chr), "WARNING")
      next
    }
    
    log_message(sprintf("Loading chr%d...", chr))
    
    tryCatch({
      qtl_chr <- as.data.table(arrow::read_parquet(parquet_file))
      qtl_chr <- qtl_chr[phenotype_id %in% gene_ids]
      
      if (nrow(qtl_chr) > 0) {
        qtl_chr[, c("chr", "pos", "REF", "ALT", "rsid") := tstrsplit(variant_id, "_")]
        qtl_chr[, maf := ifelse(af < 0.5, af, 1 - af)]
        qtl_chr[, pos := as.integer(pos)]
        qtl_chr[, chr := as.integer(sub("chr", "", chr))]
        qtl_chr[, gene_id := phenotype_id]
        qtl_chr <- qtl_chr[rsid != "."]
        
        all_qtl <- rbind(all_qtl, qtl_chr)
        log_message(sprintf("  Loaded %d QTL associations for chr%d", 
                           nrow(qtl_chr), chr))
      }
    }, error = function(e) {
      log_message(sprintf("Error loading chr%d: %s", chr, e$message), "WARNING")
    })
  }
  
  log_message(sprintf("Total QTL associations loaded: %d", nrow(all_qtl)))
  
  return(all_qtl)
}

#' Load independent GWAS signals
#'
#' @param signals_file Path to file with independent signals
#' @param phenotype Character: phenotype to filter for (optional)
#' @return data.frame with signals
#' @export
load_gwas_signals <- function(signals_file, phenotype = NULL) {
  log_message(sprintf("Loading GWAS signals from: %s", basename(signals_file)))
  
  if (!file.exists(signals_file)) {
    stop(sprintf("Signals file not found: %s", signals_file))
  }
  
  tryCatch({
    signals <- fread(signals_file)
    log_message(sprintf("Loaded %d signals", nrow(signals)))
    
    if (!is.null(phenotype)) {
      signals <- signals[signals$phenotype == phenotype, ]
      log_message(sprintf("Filtered to %d signals for phenotype: %s", 
                         nrow(signals), phenotype))
    }
    
    # Ensure required columns exist
    if (!"chr" %in% colnames(signals) && "chromosome" %in% colnames(signals)) {
      signals$chr <- signals$chromosome
    }
    
    return(as.data.frame(signals))
    
  }, error = function(e) {
    log_message(sprintf("Error loading signals: %s", e$message), "ERROR")
    stop(e)
  })
}

# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

#' Validate GWAS data structure
#'
#' @param gwas_data data.frame with GWAS results
#' @return Logical: TRUE if valid, stops with error otherwise
validate_gwas_data <- function(gwas_data) {
  required_cols <- c("chr", "position", "rsid", "ea", "nea", "beta", "se", "p")
  missing_cols <- setdiff(required_cols, colnames(gwas_data))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("GWAS data missing required columns: %s", 
                paste(missing_cols, collapse = ", ")))
  }
  
  if (nrow(gwas_data) == 0) {
    stop("GWAS data is empty")
  }
  
  log_message("GWAS data validation passed")
  return(TRUE)
}

#' Validate QTL data structure
#'
#' @param qtl_data data.table with QTL results
#' @return Logical: TRUE if valid, stops with error otherwise
validate_qtl_data <- function(qtl_data) {
  required_cols <- c("gene_id", "chr", "pos", "slope", "slope_se", "pval_nominal")
  missing_cols <- setdiff(required_cols, colnames(qtl_data))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("QTL data missing required columns: %s", 
                paste(missing_cols, collapse = ", ")))
  }
  
  if (nrow(qtl_data) == 0) {
    stop("QTL data is empty")
  }
  
  log_message("QTL data validation passed")
  return(TRUE)
}

# =============================================================================
# CLI INTERFACE (optional - for standalone execution)
# =============================================================================

if (sys.nframe() == 0) {
  # This code runs only when script is executed directly
  log_message("Data loader module loaded successfully")
  log_message("Available functions:")
  log_message("  - load_variant_annotation()")
  log_message("  - load_gwas_vcf()")
  log_message("  - extract_gwas_windows()")
  log_message("  - load_qtl_permutations()")
  log_message("  - load_qtl_nominal()")
  log_message("  - load_gwas_signals()")
}
