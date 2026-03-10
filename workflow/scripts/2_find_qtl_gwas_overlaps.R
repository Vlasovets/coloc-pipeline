#!/usr/bin/env Rscript
# Stage 2: Find QTL-GWAS Overlaps
# This script identifies genomic regions where QTL and GWAS signals overlap

# Early error wrapper to capture issues before sink()
options(error = function() {
  cat("\n=== R ERROR TRACE ===\n", file=stderr())
  traceback(2)
  if (!interactive()) quit(status = 2)
})

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

# Configure local scratch for temporary files (faster I/O)
local_scratch <- Sys.getenv("TMPDIR", "/tmp")
if(dir.exists(local_scratch)) {
  cat("Using local scratch:", local_scratch, "\n", file=stderr())
  # Ensure R uses this for temporary files
  Sys.setenv(TMPDIR = local_scratch)
  tempdir(check = TRUE)  # Force R to recognize new TMPDIR
} else {
  cat("Warning: Local scratch not available, using default temp\n", file=stderr())
}

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
  
  # Extract and prepare GWAS data for Stage 3
  cat("Extracting GWAS data in overlapping windows for Stage 3...\n")
  suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(gwasvcf)
  })
  
  # Load variant annotation for allele flipping
  cat("Loading variant annotation file...\n")
  variant_ann <- fread(snakemake@params[["variant_ann"]], data.table=FALSE)
  cat("  Loaded", nrow(variant_ann), "variants\n")
  cat("  Columns:", paste(colnames(variant_ann), collapse=", "), "\n")
  
  variant_ann$chr <- as.integer(variant_ann$chr)
  variant_ann <- variant_ann[!is.na(variant_ann$chr) & nchar(variant_ann$ref)==1 & nchar(variant_ann$alt)==1,]
  cat("  After filtering to biallelic SNPs:", nrow(variant_ann), "variants\n")
  
  # Check if variant_id_chrpos column exists
  if(!"variant_id_chrpos" %in% colnames(variant_ann)) {
    # Create it if missing
    cat("  Creating variant_id_chrpos column\n")
    variant_ann$variant_id_chrpos <- paste(variant_ann$chr, variant_ann$position, variant_ann$alt, sep="_")
  }
  cat("  Variant annotation ready with", nrow(variant_ann), "variants\n")
  
  # Read GWAS VCF for overlapping regions
  gwas_vcf_file <- snakemake@input[["gwas_vcf"]]
  
  # Copy GWAS VCF to local scratch for faster I/O (5-10x speedup)
  if(dir.exists(local_scratch) && !startsWith(gwas_vcf_file, local_scratch)) {
    cat("Copying GWAS VCF to local scratch for faster access...\n")
    local_vcf <- file.path(local_scratch, basename(gwas_vcf_file))
    local_vcf_index <- paste0(local_vcf, ".tbi")
    
    # Copy VCF and index
    file.copy(gwas_vcf_file, local_vcf, overwrite = TRUE)
    file.copy(paste0(gwas_vcf_file, ".tbi"), local_vcf_index, overwrite = TRUE)
    
    cat("  VCF copied to:", local_vcf, "\n")
    gwas_vcf_file <- local_vcf  # Use local copy
  }
  
  # Checkpoint file for resuming (saves 2 hours of VCF reading)
  checkpoint_file <- paste0(snakemake@output[["gwas_data"]], ".checkpoint_raw.rds")
  
  if(file.exists(checkpoint_file)) {
    cat("*** RESUMING: Found checkpoint file, loading pre-extracted GWAS data ***\n")
    GWAS_associations <- readRDS(checkpoint_file)
    cat("Loaded", nrow(GWAS_associations), "GWAS variants from checkpoint\n")
  } else {
    cat("No checkpoint found, extracting GWAS data from VCF...\n")
    GWAS_associations <- data.table()
    
    total_signals <- nrow(gwas_signals)
    
    # TEST MODE: Limit signals for faster debugging
    test_mode <- Sys.getenv("COLOC_TEST_MODE", "false")
    if(test_mode == "true") {
      test_limit <- min(10, total_signals)
      cat("*** TEST MODE: Processing only first", test_limit, "signals (set COLOC_TEST_MODE=false for full run) ***\n")
      total_signals <- test_limit
    }
    
    cat("Processing", total_signals, "GWAS signals...\n")
  
  for(i in 1:total_signals) {
    # Progress indicator
    if(i %% 50 == 0 || i == total_signals) {
      cat("  Progress:", i, "of", total_signals, "signals processed\n")
    }
    
    chr <- gwas_signals$chromosome[i]
    pos <- gwas_signals$position[i]
    window_start <- max(1, pos - window_size)
    window_end <- pos + window_size
    
    # Extract region from VCF (without 'chr' prefix for this VCF)
    region_gr <- GRanges(seqnames = as.character(chr),
                        ranges = IRanges(start = window_start, end = window_end))
    
    tryCatch({
      vcf <- readVcf(gwas_vcf_file, "hg38", param = region_gr)
      cat("    Region chr", chr, ":", window_start, "-", window_end, " - Found", length(vcf), "variants\n")
      
      if(length(vcf) > 0) {
        # Extract summary statistics from FORMAT fields (not INFO)
        geno_fields <- geno(vcf)
        fixed_fields <- fixed(vcf)
        
        cat("      DEBUG: Available GENO (FORMAT) fields:", paste(names(geno_fields), collapse=", "), "\n")
        
        # Handle ALT field (might be DNAStringSetList)
        alt_alleles <- sapply(fixed_fields$ALT, function(x) as.character(x)[1])
        
        region_data <- data.table(
          rsid = names(vcf),
          chr = as.integer(as.character(seqnames(vcf))),
          position = start(vcf),
          ea = alt_alleles,
          nea = as.character(fixed_fields$REF),
          beta = as.numeric(geno_fields$ES[,1]),
          se = as.numeric(geno_fields$SE[,1]),
          p = as.numeric(geno_fields$LP[,1]),  # -log10(p)
          eaf = as.numeric(geno_fields$AF[,1])
        )
        
        cat("      DEBUG: nrow(region_data):", nrow(region_data), "\n")
        
        # Convert log10p back to p-value
        region_data$p <- 10^(-region_data$p)
        
        cat("      Extracted", nrow(region_data), "variants, adding to dataset\n")
        GWAS_associations <- rbind(GWAS_associations, region_data, fill=TRUE)
        cat("      Total now:", nrow(GWAS_associations), "\n")
      } else {
        cat("      No variants in this region\n")
      }
    }, error = function(e) {
      cat("    ERROR extracting region chr", chr, ":", window_start, "-", window_end, "\n")
      cat("    Error message:", conditionMessage(e), "\n")
      cat("    Error call:", paste(deparse(conditionCall(e)), collapse=" "), "\n")
    })
  }
  
  cat("Extracted", nrow(GWAS_associations), "GWAS variants from VCF\n")
  
  # Save checkpoint before allele flipping (so we can resume if that step fails)
  if(!file.exists(checkpoint_file)) {
    cat("Saving checkpoint...\n")
    saveRDS(GWAS_associations, checkpoint_file)
    cat("Checkpoint saved to:", checkpoint_file, "\n")
  }
}
  
  # Check if we have data
  if(nrow(GWAS_associations) == 0) {
    cat("WARNING: No GWAS data extracted. Creating empty dataset.\n")
    # Save empty dataset
    output_gwas_data <- snakemake@output[["gwas_data"]]
    save(GWAS_associations, file = output_gwas_data)
    cat("Saved empty GWAS data to:", output_gwas_data, "\n\n")
  } else {
    # Perform allele flipping
    cat("Performing allele flipping...\n")
    cat("  GWAS variants:", nrow(GWAS_associations), "\n")
    cat("  Reference annotation variants:", nrow(variant_ann), "\n")
    
    # Ensure numeric columns are actually numeric
    cat("  Converting columns to numeric...\n")
    setDT(GWAS_associations)
    GWAS_associations[, c("beta", "se", "p", "eaf") := lapply(.SD, as.numeric), .SDcols = c("beta", "se", "p", "eaf")]
    
    GWAS_associations$variant_id_chrpos = paste(GWAS_associations$chr, GWAS_associations$position, GWAS_associations$ea, sep="_")
    GWAS_associations$variant_id_rev = paste(GWAS_associations$chr, GWAS_associations$position, GWAS_associations$nea, sep="_")
    
    cat("  Checking allele orientation...\n")
    tmp = ifelse(GWAS_associations$variant_id_chrpos %in% variant_ann$variant_id_chrpos, "orig",
                 ifelse(GWAS_associations$variant_id_rev %in% variant_ann$variant_id_chrpos, "rev", NA))
    
    cat("  Orientation: orig=", sum(tmp == "orig", na.rm=TRUE), 
        ", rev=", sum(tmp == "rev", na.rm=TRUE),
        ", NA=", sum(is.na(tmp)), "\n")
    
    cat("  Creating variant IDs for coloc...\n")
    # Use data.table operations instead of nested ifelse (more robust for large datasets)
    setDT(GWAS_associations)
    GWAS_associations[, varid_for_coloc := NA_character_]
    GWAS_associations[tmp == "orig" & !is.na(tmp), varid_for_coloc := variant_id_chrpos]
    GWAS_associations[tmp == "rev" & !is.na(tmp), varid_for_coloc := variant_id_rev]
    
    cat("  Flipping betas...\n")
    # Flip betas for reversed alleles (modify in place)
    GWAS_associations[tmp == "rev" & !is.na(tmp), beta := -1 * beta]
    
    cat("  Assigning effect alleles...\n")
    GWAS_associations[, ea_coloc := NA_character_]
    GWAS_associations[, nea_coloc := NA_character_]
    GWAS_associations[tmp == "orig" & !is.na(tmp), `:=`(ea_coloc = ea, nea_coloc = nea)]
    GWAS_associations[tmp == "rev" & !is.na(tmp), `:=`(ea_coloc = nea, nea_coloc = ea)]
    
    GWAS_associations[, ea := ea_coloc]
    GWAS_associations[, nea := nea_coloc]
    GWAS_associations <- GWAS_associations[!is.na(varid_for_coloc)]
    
    cat("  After filtering to matched variants:", nrow(GWAS_associations), "\n")
    
    # Save GWAS data
    output_gwas_data <- snakemake@output[["gwas_data"]]
    save(GWAS_associations, file = output_gwas_data)
    cat("Saved GWAS data (", nrow(GWAS_associations), " SNPs) to:", output_gwas_data, "\n\n")
    
    # Clean up checkpoint file on success
    if(file.exists(checkpoint_file)) {
      file.remove(checkpoint_file)
      cat("Removed checkpoint file\n")
    }
  }
  
  cat("Stage 2 completed successfully at:", as.character(Sys.time()), "\n")
  
}, error = function(e) {
  cat("\nERROR in Stage 2:\n")
  cat(conditionMessage(e), "\n")
  cat(paste(conditionCall(e), collapse = "\n"), "\n")
  quit(status = 1)
})

sink()
