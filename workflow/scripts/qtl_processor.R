#!/usr/bin/env Rscript

#' QTL Processor Module for Colocalization Pipeline
#'
#' This module handles QTL-GWAS overlap detection and data preparation:
#' - Create genomic windows around QTL lead variants
#' - Find overlaps between QTL and GWAS regions
#' - Harmonize alleles and create variant IDs
#' - Filter and prepare data for colocalization
#'
#' Refactored from src/3_run_coloc_abf.R with improvements:
#' - Cleaner separation of concerns
#' - Better memory management
#' - Progress tracking and logging
#' - Reproducible window creation logic

# =============================================================================
# DEPENDENCIES
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
})

# =============================================================================
# WINDOW CREATION FUNCTIONS
# =============================================================================

#' Create genomic windows around GWAS signals
#'
#' @param signals data.frame with columns: rsid, chr, position
#' @param window_size Numeric: window size in base pairs (default 1Mb)
#' @return list with signals data.frame and GenomicRanges object
#' @export
create_gwas_windows <- function(signals, window_size = 1e6) {
  log_message(sprintf("Creating genomic windows of %d bp around %d GWAS signals", 
                     window_size, nrow(signals)))
  
  # Calculate window coordinates
  signals$start <- ifelse(signals$position - window_size >= 1, 
                          signals$position - window_size, 1)
  signals$end <- signals$position + window_size
  
  # Create GenomicRanges object
  signals_gr <- makeGRangesFromDataFrame(
    signals,
    keep.extra.columns = FALSE,
    ignore.strand = FALSE,
    seqinfo = NULL,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    starts.in.df.are.0based = FALSE
  )
  
  log_message(sprintf("Created %d genomic windows", length(signals_gr)))
  
  return(list(signals = signals, signals_gr = signals_gr))
}

#' Create genomic windows around QTL lead variants
#'
#' @param qtl_data data.table with QTL permutation results
#' @param window_size Numeric: window size in base pairs
#' @return data.table with window coordinates added
#' @export
create_qtl_windows <- function(qtl_data, window_size = 1e6) {
  log_message(sprintf("Creating %d bp windows around %d QTL lead variants", 
                     window_size, nrow(qtl_data)))
  
  qtl_data[, mqtl_start_coord := pmax(1, pos - window_size)]
  qtl_data[, mqtl_end_coord := pos + window_size]
  
  log_message("QTL windows created")
  
  return(qtl_data)
}

# =============================================================================
# OVERLAP DETECTION
# =============================================================================

#' Find overlapping regions between QTL and GWAS windows
#'
#' @param qtl_windows data.table with QTL windows (from create_qtl_windows)
#' @param gwas_signals data.frame with GWAS signals and windows
#' @return data.table with overlapping QTL-GWAS pairs
#' @export
find_qtl_gwas_overlaps <- function(qtl_windows, gwas_signals) {
  log_message("Finding overlaps between QTL and GWAS windows...")
  
  # Convert to data.table if needed
  setDT(qtl_windows)
  setDT(gwas_signals)
  
  # Ensure chromosome is numeric
  gwas_signals[, chr := as.numeric(as.character(chr))]
  qtl_windows[, chr := as.numeric(as.character(chr))]
  
  # Create position2 column for foverlaps
  gwas_signals[, position2 := position]
  
  # Set keys for fast overlap join
  setkey(qtl_windows, chr, mqtl_start_coord, mqtl_end_coord)
  setkey(gwas_signals, chr, position, position2)
  
  # Find overlaps
  overlaps <- foverlaps(qtl_windows, gwas_signals, type = "any", 
                        nomatch = NULL, mult = "all")
  
  log_message(sprintf("Found %d QTL-GWAS overlapping regions", nrow(overlaps)))
  log_message(sprintf("  Unique QTL genes: %d", length(unique(overlaps$gene_id))))
  log_message(sprintf("  Unique GWAS signals: %d", length(unique(overlaps$rsid))))
  
  return(overlaps)
}

# =============================================================================
# DATA HARMONIZATION
# =============================================================================

#' Create variant ID for colocalization (chr_pos_allele format)
#'
#' @param chr Numeric or character vector of chromosomes
#' @param pos Numeric vector of positions
#' @param allele Character vector of alleles
#' @return Character vector of variant IDs
#' @export
create_variant_id <- function(chr, pos, allele) {
  paste(chr, pos, toupper(allele), sep = "_")
}

#' Harmonize QTL and GWAS data for colocalization
#'
#' @param qtl_data data.table with QTL summary statistics
#' @param gwas_data data.frame with GWAS summary statistics
#' @param min_variants Numeric: minimum number of overlapping variants
#' @return list with harmonized QTL and GWAS data
#' @export
harmonize_qtl_gwas <- function(qtl_data, gwas_data, min_variants = 100) {
  log_message("Harmonizing QTL and GWAS data...")
  
  # Create variant IDs for matching
  if (!"varid_for_coloc" %in% colnames(qtl_data)) {
    qtl_data[, varid_for_coloc := paste(chr, pos, ALT, sep = "_")]
  }
  
  if (!"varid_for_coloc" %in% colnames(gwas_data)) {
    gwas_data$varid_for_coloc <- paste(gwas_data$chr, gwas_data$position, 
                                        toupper(gwas_data$ea), sep = "_")
  }
  
  # Filter QTL data to only variants in GWAS
  qtl_filtered <- qtl_data[varid_for_coloc %in% gwas_data$varid_for_coloc]
  
  log_message(sprintf("  Overlapping variants: %d", nrow(qtl_filtered)))
  
  if (nrow(qtl_filtered) < min_variants) {
    log_message(sprintf("  Insufficient overlap (%d < %d), skipping", 
                       nrow(qtl_filtered), min_variants), "WARNING")
    return(NULL)
  }
  
  # Keep only SNPs (biallelic single nucleotide variants)
  if (all(c("REF", "ALT") %in% colnames(qtl_filtered))) {
    qtl_filtered <- qtl_filtered[REF %in% c("A", "T", "C", "G") & 
                                  ALT %in% c("A", "T", "C", "G")]
    log_message(sprintf("  After filtering to SNPs: %d variants", 
                       nrow(qtl_filtered)))
  }
  
  # Filter GWAS to match QTL variants
  gwas_filtered <- gwas_data[gwas_data$varid_for_coloc %in% 
                              qtl_filtered$varid_for_coloc, ]
  
  log_message(sprintf("Harmonization complete: %d shared variants", 
                     nrow(gwas_filtered)))
  
  return(list(
    qtl = as.data.frame(qtl_filtered),
    gwas = gwas_filtered,
    n_variants = nrow(gwas_filtered)
  ))
}

# =============================================================================
# DATA PREPARATION FOR COLOC
# =============================================================================

#' Prepare QTL data for coloc analysis
#'
#' @param qtl_data data.frame with QTL summary statistics
#' @param sample_size Numeric: QTL sample size
#' @return list in coloc format
#' @export
prepare_coloc_qtl <- function(qtl_data, sample_size) {
  
  # Ensure required columns exist
  required_cols <- c("varid_for_coloc", "maf", "slope", "slope_se", "pval_nominal")
  missing <- setdiff(required_cols, colnames(qtl_data))
  if (length(missing) > 0) {
    stop(sprintf("QTL data missing columns: %s", paste(missing, collapse = ", ")))
  }
  
  # Create coloc input list
  coloc_data <- list(
    N = sample_size,
    sdY = 1,
    MAF = qtl_data$maf,
    beta = qtl_data$slope,
    varbeta = (qtl_data$slope_se)^2,
    type = "quant",
    snp = qtl_data$varid_for_coloc
  )
  
  # Ensure MAF is <= 0.5
  coloc_data$MAF <- ifelse(coloc_data$MAF <= 0.5, coloc_data$MAF, 
                           1 - coloc_data$MAF)
  
  return(coloc_data)
}

#' Prepare GWAS data for coloc analysis
#'
#' @param gwas_data data.frame with GWAS summary statistics
#' @param sample_size Numeric or vector: GWAS sample size (total or c(cases, controls))
#' @param gwas_type Character: "quant" or "cc" (case-control)
#' @return list in coloc format
#' @export
prepare_coloc_gwas <- function(gwas_data, sample_size, gwas_type = "cc") {
  
  # Ensure required columns exist
  required_cols <- c("varid_for_coloc", "eaf", "beta", "se", "p")
  missing <- setdiff(required_cols, colnames(gwas_data))
  if (length(missing) > 0) {
    stop(sprintf("GWAS data missing columns: %s", paste(missing, collapse = ", ")))
  }
  
  # Create coloc input list
  coloc_data <- list(
    MAF = gwas_data$eaf,
    sdY = 1,
    beta = gwas_data$beta,
    varbeta = (gwas_data$se)^2,
    type = gwas_type,
    snp = gwas_data$varid_for_coloc
  )
  
  # Add sample size info
  if (gwas_type == "cc" && length(sample_size) == 2) {
    coloc_data$N <- sum(sample_size)
    coloc_data$s <- sample_size[1] / sum(sample_size)  # Proportion of cases
  } else {
    coloc_data$N <- sample_size
  }
  
  # Ensure MAF is <= 0.5
  coloc_data$MAF <- ifelse(coloc_data$MAF <= 0.5, coloc_data$MAF, 
                           1 - coloc_data$MAF)
  
  return(coloc_data)
}

# =============================================================================
# OVERLAP SUMMARY FUNCTIONS
# =============================================================================

#' Create overlap dataframe for coloc analysis
#'
#' @param gwas_signals data.frame with GWAS signals
#' @param qtl_genes data.frame with QTL genes
#' @return data.frame with overlap information
#' @export
create_overlap_df <- function(gwas_signals, qtl_genes) {
  
  # This function creates a structured overlap dataframe
  # that maps GWAS signals to QTL genes for systematic testing
  
  overlap_df <- expand.grid(
    signal_gwas.rsid = gwas_signals$rsid,
    signal_gwas.chr = gwas_signals$chr,
    signal_gwas.position = gwas_signals$position,
    cpg.id = qtl_genes$gene_id,
    cpg.pos = qtl_genes$pos,
    stringsAsFactors = FALSE
  )
  
  log_message(sprintf("Created overlap dataframe with %d tests", 
                     nrow(overlap_df)))
  
  return(overlap_df)
}

# =============================================================================
# UTILITY FUNCTIONS
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

#' Filter QTL-GWAS pairs by genomic distance
#'
#' @param overlaps data.table with QTL-GWAS overlaps
#' @param max_distance Numeric: maximum distance in bp
#' @return data.table filtered by distance
#' @export
filter_by_distance <- function(overlaps, max_distance = 1e6) {
  overlaps[, distance := abs(pos - position)]
  filtered <- overlaps[distance <= max_distance]
  
  log_message(sprintf("Filtered to %d pairs within %d bp", 
                     nrow(filtered), max_distance))
  
  return(filtered)
}

# =============================================================================
# CLI INTERFACE
# =============================================================================

if (sys.nframe() == 0) {
  log_message("QTL processor module loaded successfully")
  log_message("Available functions:")
  log_message("  - create_gwas_windows()")
  log_message("  - create_qtl_windows()")
  log_message("  - find_qtl_gwas_overlaps()")
  log_message("  - harmonize_qtl_gwas()")
  log_message("  - prepare_coloc_qtl()")
  log_message("  - prepare_coloc_gwas()")
}
