#!/usr/bin/env Rscript
# Test Stage 4: Aggregate Coloc Results
# Tests 5_postprocess_coloc.R logic with:
#   a) synthetic mock data (always runs)
#   b) real Stage 3 output if available (skipped if not present)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ── helpers ──────────────────────────────────────────────────────────────────
pass <- function(msg) cat(sprintf("[PASS] %s\n", msg))
fail <- function(msg) { cat(sprintf("[FAIL] %s\n", msg)); quit(save="no", status=1) }
skip <- function(msg) cat(sprintf("[SKIP] %s\n", msg))
section <- function(msg) cat(sprintf("\n── %s ──\n", msg))

OUTPUT_DIR  <- "/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"
GWAS_ID     <- "KNEE"
PP4_THRESH  <- 0.8
PP4_PP3_RATIO <- 2.0

cat("══════════════════════════════════════════════\n")
cat(" Stage 4 Aggregate Results — Test Suite\n")
cat("══════════════════════════════════════════════\n")

# ── helper: run the aggregation logic ────────────────────────────────────────
run_aggregation <- function(input_files, pp4_threshold, pp4_pp3_ratio) {
  coloc_res <- rbindlist(lapply(input_files, function(f) {
    if (file.exists(f) && file.size(f) > 0) {
      dt <- fread(f)
      if (!"tissue" %in% colnames(dt)) {
        dt$tissue <- sub("^[^.]+\\.(.+)\\.colocABF_results\\.txt$", "\\1", basename(f))
      }
      if (!"GWAS_ID" %in% colnames(dt)) {
        dt$GWAS_ID <- sub("^([^.]+)\\..*", "\\1", basename(f))
      }
      return(dt)
    } else {
      return(data.table())
    }
  }), fill = TRUE)
  return(coloc_res)
}

# ══════════════════════════════════════════════
# Part A: Mock data tests (always run)
# ══════════════════════════════════════════════
section("Part A: Mock data tests")

# Build a synthetic coloc result table
mock_coloc <- data.table(
  cpg            = c("ENSG00000001", "ENSG00000002", "ENSG00000003", "ENSG00000004", "ENSG00000005"),
  cpg_pos        = c(1000000L, 2000000L, 3000000L, 4000000L, 5000000L),
  gwas_signal    = c("rs001", "rs002", "rs003", "rs004", "rs005"),
  gwas_lead_snp  = c("1:1000000", "2:2000000", "3:3000000", "4:4000000", "5:5000000"),
  PP0 = c(0.001, 0.001, 0.001, 0.001, 0.001),
  PP1 = c(0.100, 0.050, 0.010, 0.005, 0.002),
  PP2 = c(0.050, 0.030, 0.010, 0.005, 0.002),
  PP3 = c(0.049, 0.119, 0.029, 0.089, 0.395),
  PP4 = c(0.800, 0.800, 0.950, 0.900, 0.600),
  nvariants      = c(500L, 600L, 700L, 800L, 900L),
  tissue         = "high_grade_cartilage",
  GWAS_ID        = "KNEE"
)

# Write mock file
mock_dir  <- file.path(OUTPUT_DIR, "test_stage4")
mock_file <- file.path(mock_dir, "KNEE.high_grade_cartilage.colocABF_results.txt")
dir.create(mock_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(mock_coloc, mock_file, sep = "\t")
pass(sprintf("Created mock coloc file: %d rows", nrow(mock_coloc)))

# Test 1: Read and combine
section("A1. Read input files")
coloc_res <- run_aggregation(mock_file, PP4_THRESH, PP4_PP3_RATIO)
if (nrow(coloc_res) == 0) fail("Combined result is empty")
pass(sprintf("Read %d coloc rows from mock file", nrow(coloc_res)))

# Test 2: Required columns present
section("A2. Required columns")
required_cols <- c("cpg", "PP4", "PP3", "gwas_signal", "tissue", "GWAS_ID")
missing <- setdiff(required_cols, colnames(coloc_res))
if (length(missing) > 0) fail(paste("Missing columns:", paste(missing, collapse=", ")))
pass(paste("All required columns present:", paste(required_cols, collapse=", ")))

# Test 3: PP values in valid range
section("A3. PP value range")
stopifnot(all(coloc_res$PP4 >= 0 & coloc_res$PP4 <= 1))
stopifnot(all(coloc_res$PP3 >= 0 & coloc_res$PP3 <= 1))
pass("PP3 and PP4 values in [0, 1]")

# Test 4: Significant filtering (PP4 >= 0.8 OR PP4 > 0.6 AND PP4/PP3 > 2)
section("A4. Significant coloc filtering")
coloc_sig <- coloc_res[PP4 >= PP4_THRESH | (PP4 > 0.6 & PP4 < PP4_THRESH & PP4/PP3 > PP4_PP3_RATIO)]
# Expected: rows 1 (PP4=0.80), 2 (PP4=0.80), 3 (PP4=0.95), 4 (PP4=0.90) => 4 hits
# Row 5: PP4=0.60, PP4/PP3=0.6/0.395=1.52 < 2 => not significant
expected_sig <- 4L
if (nrow(coloc_sig) != expected_sig) {
  fail(sprintf("Expected %d significant colocalizations, got %d", expected_sig, nrow(coloc_sig)))
}
pass(sprintf("Filtering correct: %d/%d rows pass threshold", nrow(coloc_sig), nrow(coloc_res)))

# Test 5: Summary statistics
section("A5. Summary statistics")
summary_dt <- coloc_res %>%
  group_by(GWAS_ID, tissue) %>%
  summarise(
    n_tests   = n(),
    n_genes   = n_distinct(cpg),
    n_signals = n_distinct(gwas_signal),
    mean_PP4  = mean(PP4, na.rm = TRUE),
    .groups   = "drop"
  )
if (nrow(summary_dt) == 0) fail("Summary is empty")
if (summary_dt$n_tests[1] != nrow(mock_coloc)) fail("Summary n_tests mismatch")
pass(sprintf("Summary correct: %d tests, mean PP4=%.3f", summary_dt$n_tests[1], summary_dt$mean_PP4[1]))

# Test 6: Empty file handling
section("A6. Empty file handling")
empty_file <- file.path(mock_dir, "KNEE.synovium.colocABF_results.txt")
fwrite(data.table(), empty_file, sep = "\t")
combined <- run_aggregation(c(mock_file, empty_file), PP4_THRESH, PP4_PP3_RATIO)
if (nrow(combined) != nrow(mock_coloc)) fail("Empty file handling failed: wrong row count")
pass("Empty input file handled gracefully")

# Test 7: Output files can be written
section("A7. Write output files")
out_all <- file.path(mock_dir, paste0(GWAS_ID, "_all_coloc_results.txt"))
out_sig <- file.path(mock_dir, paste0(GWAS_ID, "_significant_coloc_results.txt"))
fwrite(coloc_res, out_all, sep = "\t")
fwrite(coloc_sig, out_sig, sep = "\t")
if (!file.exists(out_all) || file.size(out_all) == 0) fail("all_coloc_results.txt not written")
if (!file.exists(out_sig)) fail("significant_coloc_results.txt not written")
pass(sprintf("Output written: all=%d rows, significant=%d rows", nrow(coloc_res), nrow(coloc_sig)))

cat("\n[SUMMARY] Part A (mock data): ALL TESTS PASSED\n")

# ══════════════════════════════════════════════
# Part B: Real Stage 3 output (if available)
# ══════════════════════════════════════════════
section("Part B: Real Stage 3 output")

real_coloc_file <- file.path(
  OUTPUT_DIR, "coloc_abf",
  paste0(GWAS_ID, ".high_grade_cartilage.colocABF_results.txt")
)

if (!file.exists(real_coloc_file)) {
  skip(sprintf("Stage 3 output not yet available: %s", real_coloc_file))
  skip("Skipping Part B — run after Stage 3 completes")
} else {
  real_size <- file.size(real_coloc_file)
  if (real_size == 0) {
    skip("Stage 3 output file is empty — skipping Part B")
  } else {
    # B1: Read real output
    section("B1. Read real Stage 3 output")
    real_res <- run_aggregation(real_coloc_file, PP4_THRESH, PP4_PP3_RATIO)
    pass(sprintf("Read %d coloc rows from real Stage 3 output (%s)",
                 nrow(real_res), format(real_size, big.mark=",")))

    # B2: Check columns
    section("B2. Column check")
    missing_b <- setdiff(c("cpg", "PP4", "PP3", "gwas_signal"), colnames(real_res))
    if (length(missing_b) > 0) fail(paste("Missing real columns:", paste(missing_b, collapse=", ")))
    pass("Real output has required columns")

    # B3: PP ranges
    section("B3. PP value validity")
    stopifnot(all(real_res$PP4 >= 0 & real_res$PP4 <= 1))
    pass("Real PP4 values in [0, 1]")

    # B4: Filter significant
    section("B4. Significant filter on real data")
    real_sig <- real_res[PP4 >= PP4_THRESH | (PP4 > 0.6 & PP4 < PP4_THRESH & PP4/PP3 > PP4_PP3_RATIO)]
    pass(sprintf("Significant colocalizations: %d / %d total", nrow(real_sig), nrow(real_res)))

    # B5: Write to results dir
    section("B5. Write real aggregated output")
    results_dir <- file.path(OUTPUT_DIR, "results")
    dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
    real_out_all <- file.path(results_dir, paste0(GWAS_ID, "_all_coloc_results.txt"))
    real_out_sig <- file.path(results_dir, paste0(GWAS_ID, "_significant_coloc_results.txt"))
    fwrite(real_res, real_out_all, sep = "\t")
    fwrite(real_sig, real_out_sig, sep = "\t")
    pass(sprintf("Written: %s (%d rows)", basename(real_out_all), nrow(real_res)))
    pass(sprintf("Written: %s (%d rows)", basename(real_out_sig), nrow(real_sig)))

    cat("\n[SUMMARY] Part B (real data): ALL TESTS PASSED\n")
  }
}

# ── Final summary ─────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════\n")
cat(" Stage 4 test suite COMPLETE\n")
cat("══════════════════════════════════════════════\n")
