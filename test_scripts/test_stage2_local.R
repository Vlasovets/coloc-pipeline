#!/usr/bin/env Rscript
# Test Stage 2 locally

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

# Use absolute path for testing
script_dir <- "/home/itg/oleg.vlasovets/projects/coloc-pipeline/workflow/scripts"

cat("=== Testing helper script sourcing ===\n")
cat("Sourcing coloc_helpers.R...\n")
source(file.path(script_dir, "coloc_helpers.R"))
cat("✓ coloc_helpers.R loaded\n")

cat("Sourcing data_loader.R...\n")
source(file.path(script_dir, "data_loader.R"))
cat("✓ data_loader.R loaded\n")

cat("Sourcing qtl_processor.R...\n")
source(file.path(script_dir, "qtl_processor.R"))
cat("✓ qtl_processor.R loaded\n")

# Test loading data
cat("\n=== Testing data loading ===\n")
trait <- "KNEE"
tissue <- "fat_pad"

gwas_signals_file <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/GO2_index_signals_b38.csv"
qtl_perm_file <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/eQTL_results/fat_pad/cov_df_fat_pad_K06/cis_permutation/egenes_df.txt"

cat("Loading GWAS signals...\n")
gwas_signals <- load_gwas_signals(
  signals_file = gwas_signals_file,
  phenotype = trait
)
cat("✓ Loaded", nrow(gwas_signals), "GWAS signals\n")

cat("Loading QTL permutations...\n")
qtl_perms <- load_qtl_permutations(
  perm_file = qtl_perm_file,
  qval_threshold = 0.05
)
cat("✓ Loaded", nrow(qtl_perms), "QTL permutations\n")

# Test window creation
cat("\n=== Testing window creation ===\n")
gwas_windows <- create_gwas_windows(gwas_signals, window_size = 1e6)
cat("✓ Created", nrow(gwas_windows), "GWAS windows\n")

qtl_windows <- create_qtl_windows(qtl_perms, window_size = 1e6)
cat("✓ Created", nrow(qtl_windows), "QTL windows\n")

# Test overlap detection
cat("\n=== Testing overlap detection ===\n")
overlaps <- find_qtl_gwas_overlaps(qtl_windows, gwas_windows)
cat("✓ Found", nrow(overlaps), "overlaps\n")

cat("\n=== ALL TESTS PASSED ===\n")
