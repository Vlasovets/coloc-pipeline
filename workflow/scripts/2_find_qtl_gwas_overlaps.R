#!/usr/bin/env Rscript
# Stage 2: Find QTL-GWAS Overlaps
# This script identifies genomic regions where QTL and GWAS signals overlap

# Early error wrapper to capture issues before sink()
options(error = function() {
  cat("\n=== R ERROR TRACE ===\n", file=stderr())
  traceback(2, file=stderr())
  if (!interactive()) quit(status = 2)
})

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

# Get script directory - Snakemake should provide this
if (exists("snakemake")) {
  script_dir <- snakemake@scriptdir
} else {
  # Fallback for testing
  script_dir <- "/home/itg/oleg.vlasovets/projects/coloc-pipeline/workflow/scripts"
}

cat("Script directory:", script_dir, "\n", file=stderr())

# Source modular functions with error handling
tryCatch({
  cat("Sourcing coloc_helpers.R...\n", file=stderr())
  source(file.path(script_dir, "coloc_helpers.R"))
  cat("Sourcing data_loader.R...\n", file=stderr())
  source(file.path(script_dir, "data_loader.R"))
  cat("Sourcing qtl_processor.R...\n", file=stderr())
  source(file.path(script_dir, "qtl_processor.R"))
  cat("All helpers loaded successfully\n", file=stderr())
}, error = function(e) {
  cat("FATAL: Failed to source helper scripts\n", file=stderr())
  cat("Error:", conditionMessage(e), "\n", file=stderr())
  quit(status = 1)
})

# Get parameters from Snakemake
trait <- snakemake@wildcards[["trait"]]
tissue <- snakemake@wildcards[["tissue"]]
output_overlaps <- snakemake@output[["overlaps"]]
output_qtl_subset <- snakemake@output[["qtl_subset"]]
log_file <- snakemake@log[[1]]

# Redirect output to log with early error handling
tryCatch({
  sink(log_file, append = FALSE, split = TRUE)
  cat("Stage 2: Find QTL-GWAS Overlaps\n")
  cat("Script directory used:", script_dir, "\n")
  cat("Trait:", trait, "\n")
  cat("Tissue:", tissue, "\n")
  cat("Started at:", as.character(Sys.time()), "\n\n")
}, error = function(e) {
  # If even sink fails, write to stderr
  message("FATAL: Cannot open log file: ", log_file)
  message("Error: ", conditionMessage(e))
  quit(status = 1)
})

tryCatch({
  
  # Load GWAS signals
  cat("Loading GWAS signals...\n")
  gwas_signals <- load_gwas_signals(
    signals_file = snakemake@config[["gwas_signals_file"]],
    phenotype = trait
  )
  cat("Loaded", nrow(gwas_signals), "GWAS signals\n\n")
  
  # Load QTL permutations
  cat("Loading QTL permutation results...\n")
  qtl_perm_file <- snakemake@config[["qtl_permutation_files"]][[tissue]]
  qtl_perms <- load_qtl_permutations(
    perm_file = qtl_perm_file,
    qval_threshold = snakemake@config[["qval_threshold"]]
  )
  cat("Loaded", nrow(qtl_perms), "significant QTLs (qval <", 
      snakemake@config[["qval_threshold"]], ")\n\n")
  
  # Create windows around signals
  cat("Creating genomic windows...\n")
  window_size <- snakemake@config[["window_size"]]
  gwas_windows <- create_gwas_windows(gwas_signals, window_size = window_size)
  qtl_windows <- create_qtl_windows(qtl_perms, window_size = window_size)
  cat("GWAS windows:", nrow(gwas_windows), "\n")
  cat("QTL windows:", nrow(qtl_windows), "\n\n")
  
  # Find overlaps
  cat("Finding QTL-GWAS overlaps...\n")
  overlaps <- find_qtl_gwas_overlaps(qtl_windows, gwas_windows)
  cat("Found", nrow(overlaps), "overlapping regions\n\n")
  
  # Save results
  cat("Saving results...\n")
  save(overlaps, gwas_windows, qtl_windows, file = output_overlaps)
  cat("Saved overlaps to:", output_overlaps, "\n")
  
  # Save list of genes for next stage
  gene_list <- unique(overlaps$gene_id)
  fwrite(data.table(gene_id = gene_list), output_qtl_subset, sep = "\t")
  cat("Saved", length(gene_list), "genes to:", output_qtl_subset, "\n\n")
  
  cat("Stage 2 completed successfully at:", as.character(Sys.time()), "\n")
  
}, error = function(e) {
  cat("\nERROR in Stage 2:\n")
  cat(conditionMessage(e), "\n")
  cat(paste(conditionCall(e), collapse = "\n"), "\n")
  quit(status = 1)
})

sink()
