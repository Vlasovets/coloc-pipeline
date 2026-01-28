# Change allele2 to be also a SNP (%in% c("A", "T", "C", "G"))
# Use CPTID instead of SNP

.libPaths(c("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/R", .libPaths()))
project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc"
data.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/"
setwd(project.path)

source("scripts/Coloc_helper_functions.R")

library(data.table)
library(VariantAnnotation)
library(gwasvcf)
library(magrittr)
library("GenomicRanges")
library(rtracklayer)
library(coloc)
library(dplyr)
library(parallel)

cat("\n")
cat("========================================\n")
cat("GWAS VCF Conversion Pipeline\n")
cat("========================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Traits to process: KNEE, TKR, ALLOA\n")
cat("========================================\n\n")

cat("Step 1: Loading variant annotation file...\n")
variant_ann <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/variant_ann_hg38.txt.gz")
cat("✓ Loaded", nrow(variant_ann), "variants\n\n")

cat("Step 2: Filtering variants...\n")
variant_ann$chr <- as.integer(variant_ann$chr)
original_count <- nrow(variant_ann)
variant_ann <- variant_ann[!is.na(variant_ann$chr) & nchar(variant_ann$ref)==1 & nchar(variant_ann$alt)==1,]
filtered_count <- nrow(variant_ann)
cat("✓ Retained", filtered_count, "of", original_count, "variants after filtering\n\n")


#' ============================================================================
#' INPUT VALIDATION
#' ============================================================================

#' Validate required files and directories exist before processing
#' 
#' @description Checks that all input files and output directories are accessible
#' @return NULL (stops execution with error message if validation fails)

validate_inputs <- function(VARIANT_ANNOTATION_FILE ,
                            GWAS_SUMSTATS_DIR,
                            OUTPUT_PATH,
                            HELPER_SCRIPT = "scripts/Coloc_helper_functions.R") {

  # Check variant annotation file exists
  if (!file.exists(VARIANT_ANNOTATION_FILE)) {
    stop("ERROR: Variant annotation file not found: ", VARIANT_ANNOTATION_FILE)
  }
  message("✓ Variant annotation file found")

  # Check GWAS summary statistics directory exists
  if (!dir.exists(GWAS_SUMSTATS_DIR)) {
    stop("ERROR: GWAS summary statistics directory not found: ", GWAS_SUMSTATS_DIR)
  }
  message("✓ GWAS summary statistics directory found")

  # Check output directory exists (create if not)
  if (!dir.exists(OUTPUT_PATH)) {
    message("Creating output directory: ", OUTPUT_PATH)
    dir.create(OUTPUT_PATH, recursive = TRUE)
  }
  message("✓ Output directory ready")

  # Check helper functions script exists, and source it
  if (!file.exists(HELPER_SCRIPT)) {
    stop("ERROR: Helper functions script not found: ", HELPER_SCRIPT)
  }
  source(HELPER_SCRIPT)
  message("✓ Helper functions script found")

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
validate_inputs(VARIANT_ANNOTATION_FILE  ,
                            GWAS_SUMSTATS_DIR,
                            OUTPUT_PATH,
                            HELPER_SCRIPT = "scripts/Coloc_helper_functions.R")

# =================================================
# Load the reference annotation
# =================================================

# Load variant annotation with error handling
variant_ann <- tryCatch({
  message("Loading variant annotation file...")
  fread(VARIANT_ANNOTATION_FILE)
}, error = function(e) {
  stop("ERROR: Failed to load variant annotation: ", e$message)
})
variant_ann$chr <- as.integer(variant_ann$chr)
variant_ann <- variant_ann[!is.na(variant_ann$chr) & nchar(variant_ann$ref)==1 & nchar(variant_ann$alt)==1,]


for(trait in c("KNEE", "TKR", "ALLOA")){
  
  cat("========================================\n")
  cat("Processing trait:", trait, "\n")
  cat("========================================\n")
  
  ncases <- ifelse(trait=="KNEE", 172256, ifelse(trait=="TKR", 48161, 489952))
  ncontrols <- ifelse(trait=="KNEE", 1144244, ifelse(trait=="TKR", 958463, 1471094))
  GWAS_n <- c(ncases, ncontrols)
  
  cat("Sample size:\n")
  cat("  - Cases:", formatC(ncases, format="d", big.mark=","), "\n")
  cat("  - Controls:", formatC(ncontrols, format="d", big.mark=","), "\n")
  cat("  - Total:", formatC(sum(GWAS_n), format="d", big.mark=","), "\n\n")
  
  GWASfile <- paste0("/lustre/groups/itg/teams/zeggini/projects/GO2/GO2SummStats/ALL.MAIN/GO.FILTER.GW.final.meta.results.ALL.", trait, ".FULL.MAFless0.01.Nupdated.txt.gz")
  
  cat("Input file:", basename(GWASfile), "\n")
  
  cat("Converting to VCF format...\n")
  make_vcf(GWASfile=GWASfile,
           chrom="CHR",
           pos="POS",
           nea="NEA",
           ea="EA",
           snp="CPTID",
           ea_af="EAF",
           effect="BETA",
           se="SE",
           pval="P",
           hg="hg19", # hg build of summary stats !! For GO2 it's hg19 so change here and then change the want to lift over here (not sure where the chain file has to be)
           lift_down=FALSE, # if lift_down=TRUE then hg38 is converted to hg19. If FALSE, hg19 is converted to hg38. Ignored if WantToLiftOver is set ot FALSE
           WantToLiftOver=TRUE, # Whether or not you want to lift over the coordinates
           GWAS_n=GWAS_n, # a vector of one or two elements. If quantitative trait, total sample size, if case control, number of cases and controls
           variant_ann=variant_ann, # reference file to map the missing values
           # output=paste0("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/GO2_sumstats/GO2_b38_", trait, "_ody_cptid.vcf"))
           output=paste0("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/GO2_sumstats/GO2_b38_", trait, "_ana_new.vcf"))
  
  cat("✓ Completed conversion for", trait, "\n\n")
}

cat("========================================\n")
cat("Pipeline completed successfully!\n")
cat("========================================\n")
cat("All 3 traits processed\n")
cat("Output directory: /lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/GO2_sumstats/\n")
cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n")
