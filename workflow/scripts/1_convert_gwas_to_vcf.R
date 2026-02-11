#' GWAS VCF Conversion Script for Snakemake Pipeline
#'
#' This script converts GWAS summary statistics for osteoarthritis (OA) traits
#' into standardized VCF format with liftover from hg19 to hg38.
#'
#' This version is designed to work with Snakemake and receives parameters
#' via the snakemake object.

# =============================================================================
# LOAD DEPENDENCIES
# =============================================================================

# Install gwasvcf from GitHub if not available
if (!requireNamespace("gwasvcf", quietly = TRUE)) {
  message("Installing gwasvcf from GitHub...")
  remotes::install_github("MRCIEU/gwasvcf", upgrade = "never")
}

library(data.table)
library(VariantAnnotation)
library(gwasvcf)
library(magrittr)
library(GenomicRanges)
library(rtracklayer)
library(coloc)
library(dplyr)

# =============================================================================
# GET SNAKEMAKE PARAMETERS
# =============================================================================

# Set temporary directory to writeable location
tmpdir <- file.path(dirname(snakemake@output[["vcf"]]), "..", "tmp")
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmpdir)
Sys.setenv(R_TEMPDIR = tmpdir)
options(tempdir = tmpdir)

# Access Snakemake inputs, outputs, and parameters
gwas_file <- snakemake@input[["gwas"]]
variant_annotation_file <- snakemake@input[["variant_ann"]]
output_vcf <- snakemake@output[["vcf"]]
n_cases <- snakemake@params[["n_cases"]]
n_controls <- snakemake@params[["n_controls"]]
helper_script <- snakemake@params[["helper_script"]]
genome_build <- snakemake@params[["genome_build"]]
liftover <- snakemake@params[["liftover"]]
log_file <- snakemake@log[[1]]

# Redirect output to log file
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

# =============================================================================
# LOAD HELPER FUNCTIONS
# =============================================================================

if (!file.exists(helper_script)) {
  stop("ERROR: Helper functions script not found: ", helper_script)
}
source(helper_script, local = TRUE)
message("✓ Helper functions script loaded")

# =============================================================================
# LOAD REFERENCE ANNOTATION
# =============================================================================

message("\n========================================")
message("GWAS VCF Conversion")
message("========================================")
message("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Input: ", basename(gwas_file))
message("Output: ", basename(output_vcf))
message("========================================\n")

# Load variant annotation with error handling
variant_ann <- tryCatch({
  message("Loading variant annotation file...")
  fread(variant_annotation_file)
}, error = function(e) {
  stop("ERROR: Failed to load variant annotation: ", e$message)
})

message("Filtering variants...")
original_count <- nrow(variant_ann)
variant_ann$chr <- as.integer(variant_ann$chr)
variant_ann <- variant_ann[!is.na(variant_ann$chr) & nchar(variant_ann$ref) == 1 & nchar(variant_ann$alt) == 1, ]
filtered_count <- nrow(variant_ann)
message("✓ Retained ", filtered_count, " of ", original_count, " variants after filtering\n")

# =============================================================================
# CONVERT GWAS TO VCF
# =============================================================================

gwas_n <- c(n_cases, n_controls)

message("Sample size:")
message("  - Cases: ", formatC(n_cases, format = "d", big.mark = ","))
message("  - Controls: ", formatC(n_controls, format = "d", big.mark = ","))
message("  - Total: ", formatC(sum(gwas_n), format = "d", big.mark = ","), "\n")

message("Converting to VCF format...")

# Convert GWAS summary statistics to VCF format
make_vcf(
  GWASfile = gwas_file,
  chrom = "CHR",
  pos = "POS",
  nea = "NEA",
  ea = "EA",
  snp = "CPTID",
  ea_af = "EAF",
  effect = "BETA",
  se = "SE",
  pval = "P",
  WantToLiftOver = liftover,
  GWAS_n = gwas_n,
  variant_ann = variant_ann,
  output = gsub("\\.gz$", "", output_vcf)  # Remove .gz if present, make_vcf handles compression
)

message("\n========================================")
message("✓ Conversion completed successfully!")
message("========================================")

# Close log file
sink(type = "message")
sink(type = "output")
close(log_con)
