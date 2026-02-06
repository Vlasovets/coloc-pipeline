#!/usr/bin/env Rscript
# Minimal test to verify R environment and sourcing

cat("=== Minimal R Environment Test ===\n")

# Test 1: Load libraries
cat("\n1. Testing library loading...\n")
tryCatch({
  suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
  })
  cat("✓ Libraries loaded successfully\n")
}, error = function(e) {
  cat("✗ Error loading libraries:", conditionMessage(e), "\n")
  quit(status = 1)
})

# Test 2: Source helper scripts
cat("\n2. Testing script sourcing...\n")
script_dir <- "/home/itg/oleg.vlasovets/projects/coloc-pipeline/workflow/scripts"

tryCatch({
  cat("   Sourcing coloc_helpers.R...\n")
  source(file.path(script_dir, "coloc_helpers.R"))
  cat("   ✓ coloc_helpers.R loaded\n")
  
  cat("   Sourcing data_loader.R...\n")
  source(file.path(script_dir, "data_loader.R"))
  cat("   ✓ data_loader.R loaded\n")
  
  cat("   Sourcing qtl_processor.R...\n")
  source(file.path(script_dir, "qtl_processor.R"))
  cat("   ✓ qtl_processor.R loaded\n")
}, error = function(e) {
  cat("✗ Error sourcing scripts:", conditionMessage(e), "\n")
  traceback()
  quit(status = 1)
})

# Test 3: Check functions exist
cat("\n3. Checking if functions exist...\n")
functions_to_check <- c("load_gwas_signals", "load_qtl_permutations", 
                        "create_gwas_windows", "create_qtl_windows", 
                        "find_qtl_gwas_overlaps")
for (func in functions_to_check) {
  if (exists(func)) {
    cat("   ✓", func, "exists\n")
  } else {
    cat("   ✗", func, "NOT FOUND\n")
    quit(status = 1)
  }
}

cat("\n=== ALL TESTS PASSED ===\n")
cat("The R environment and scripts are working correctly.\n")
cat("Ready to test with actual data.\n")
