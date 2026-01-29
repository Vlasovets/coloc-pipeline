#' GWAS VCF Conversion Pipeline for GO2 OA Traits
#'
#' This script converts GWAS summary statistics for osteoarthritis (OA) traits
#' into standardized VCF format with liftover from hg19 to hg38.
#'
#' Input: GO2 GWAS summary statistics (KNEE, TKR, ALLOA traits)
#' Output: Tabix-indexed VCF files in hg38 coordinates
#' Dependencies: data.table, VariantAnnotation, gwasvcf, GenomicRanges,
#'               rtracklayer, coloc, dplyr, parallel

# TODO: Change allele2 to be also a SNP (%in% c("A", "T", "C", "G"))
# TODO: Use CPTID instead of SNP

# =============================================================================
# CONFIGURATION
# =============================================================================

# Load configuration from external file
# Users should copy config.R.example to config.R and update paths
config_file <- "config.R"
if (!file.exists(config_file)) {
  stop("Configuration file not found. Please copy config.R.example to config.R and update the paths.")
}
source(config_file)

# =============================================================================
# LOAD DEPENDENCIES
# =============================================================================

library(data.table)
library(VariantAnnotation)
library(gwasvcf)
library(magrittr)
library(GenomicRanges)
library(rtracklayer)
library(coloc)
library(dplyr)
library(parallel)

# =============================================================================
# MAIN PIPELINE
# =============================================================================

message("\n========================================")
message("GWAS VCF Conversion Pipeline")
message("========================================")
message("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Traits to process: KNEE, TKR, ALLOA")
message("========================================\n")


# =============================================================================
# INPUT VALIDATION
# =============================================================================

#' Validate required files and directories exist before processing
#'
#' @description Checks that all input files and output directories are accessible
#' @return NULL (stops execution with error message if validation fails)

validate_inputs <- function(variant_annotation_file,
                            gwas_sumstats_dir,
                            output_path,
                            helper_script = "scripts/Coloc_helper_functions.R") {
  
  # Check variant annotation file exists
  if (!file.exists(variant_annotation_file)) {
    stop("ERROR: Variant annotation file not found: ", variant_annotation_file)
  }
  message("✓ Variant annotation file found")
  
  # Check GWAS summary statistics directory exists
  if (!dir.exists(gwas_sumstats_dir)) {
    stop("ERROR: GWAS summary statistics directory not found: ", gwas_sumstats_dir)
  }
  message("✓ GWAS summary statistics directory found")
  
  # Check output directory exists (create if not)
  if (!dir.exists(output_path)) {
    message("Creating output directory: ", output_path)
    dir.create(output_path, recursive = TRUE)
  }
  message("✓ Output directory ready")
  
  # Check helper functions script exists, and source it
  if (!file.exists(helper_script)) {
    stop("ERROR: Helper functions script not found: ", helper_script)
  }
  source(helper_script, local = TRUE)
  message("✓ Helper functions script loaded")
  
  # Validate required packages are installed
  required_packages <- c("data.table", "VariantAnnotation", "gwasvcf",
                         "GenomicRanges", "rtracklayer", "coloc", "dplyr", "parallel")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("ERROR: Missing required packages: ", paste(missing_packages, collapse = ", "))
  }
  message("✓ All required packages available")
  
  message("\n--- All validation checks passed ---\n")
}

# Run validation before processing
validate_inputs(
  variant_annotation_file = VARIANT_ANNOTATION_FILE,
  gwas_sumstats_dir = GWAS_SUMSTATS_DIR,
  output_path = OUTPUT_PATH,
  helper_script = HELPER_SCRIPT
)

# =============================================================================
# LOAD REFERENCE ANNOTATION
# =============================================================================

# Load variant annotation with error handling
variant_ann <- tryCatch({
  message("Loading variant annotation file...")
  fread(VARIANT_ANNOTATION_FILE)
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
# PROCESS TRAITS
# =============================================================================

for (trait in c("KNEE", "TKR", "ALLOA")) {
  
  message("========================================")
  message("Processing trait: ", trait)
  message("========================================")
  
  # Define sample sizes for each trait
  n_cases <- ifelse(trait == "KNEE", 172256, ifelse(trait == "TKR", 48161, 489952))
  n_controls <- ifelse(trait == "KNEE", 1144244, ifelse(trait == "TKR", 958463, 1471094))
  gwas_n <- c(n_cases, n_controls)
  
  message("Sample size:")
  message("  - Cases: ", formatC(n_cases, format = "d", big.mark = ","))
  message("  - Controls: ", formatC(n_controls, format = "d", big.mark = ","))
  message("  - Total: ", formatC(sum(gwas_n), format = "d", big.mark = ","), "\n")
  
  # Construct input file path
  gwas_file <- file.path(
    GWAS_SUMSTATS_DIR,
    paste0("GO.FILTER.GW.final.meta.results.ALL.", trait, ".FULL.MAFless0.01.Nupdated.txt.gz")
  )
  
  message("Input file: ", basename(gwas_file))
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
    hg = "hg19",
    lift_down = FALSE,
    WantToLiftOver = TRUE,
    GWAS_n = gwas_n,
    variant_ann = variant_ann,
    output = file.path(OUTPUT_PATH, paste0("GO2_b38_", trait, "_ana_new.vcf"))
  )
  
  message("✓ Completed conversion for ", trait, "\n")
}

# =============================================================================
# PIPELINE COMPLETION
# =============================================================================

message("========================================")
message("Pipeline completed successfully!")
message("========================================")
message("All 3 traits processed")
message("Output directory: ", OUTPUT_PATH)
message("Finished at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("========================================")
