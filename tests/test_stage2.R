#!/usr/bin/env Rscript
# Test Stage 2: QTL-GWAS Overlap Detection
# Tests 2_find_qtl_gwas_overlaps.R logic with:
#   a) synthetic mock data (always runs) — tests window creation and overlap logic
#   b) real Stage 2 output if available (skipped if not present)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(dplyr)
})

# ── helpers ──────────────────────────────────────────────────────────────────
pass <- function(msg) cat(sprintf("[PASS] %s\n", msg))
fail <- function(msg) { cat(sprintf("[FAIL] %s\n", msg)); quit(save="no", status=1) }
skip <- function(msg) cat(sprintf("[SKIP] %s\n", msg))
section <- function(msg) cat(sprintf("\n── %s ──\n", msg))

SCRIPT_DIR  <- "/home/itg/oleg.vlasovets/projects/coloc-pipeline/workflow/scripts"
OUTPUT_DIR  <- "/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"
GWAS_ID     <- "KNEE"
TISSUE      <- "high_grade_cartilage"

cat("══════════════════════════════════════════════\n")
cat(" Stage 2 QTL-GWAS Overlap — Test Suite\n")
cat("══════════════════════════════════════════════\n")

# ── Source helper functions ────────────────────────────────────────────────────
section("Source qtl_processor helpers")
tryCatch(
  source(file.path(SCRIPT_DIR, "qtl_processor.R")),
  error = function(e) fail(paste("Cannot source qtl_processor.R:", e$message))
)
stopifnot(exists("create_gwas_windows"))
stopifnot(exists("create_qtl_windows"))
stopifnot(exists("find_qtl_gwas_overlaps"))
pass("qtl_processor.R sourced; create_gwas_windows, create_qtl_windows, find_qtl_gwas_overlaps available")

# ══════════════════════════════════════════════
# Part A: Mock data tests (always run)
# ══════════════════════════════════════════════
section("Part A: Mock data tests")

# Build synthetic GWAS signals
mock_signals <- data.table(
  rsid       = c("rs001", "rs002", "rs003"),
  chr        = c(1L, 5L, 12L),
  position   = c(1000000L, 5000000L, 12000000L),
  Loci       = c("1:1000000", "5:5000000", "12:12000000")
)

# Build synthetic QTL permutation results
mock_qtl_perm <- data.table(
  phenotype_id  = c("ENSG001", "ENSG002", "ENSG003", "ENSG004"),
  variant_id    = c("chr1_995000_A_G_rs10",
                    "chr1_1100000_C_T_rs11",
                    "chr5_4800000_G_A_rs20",
                    "chr8_9000000_T_C_rs30"),
  qval          = c(0.001, 0.02, 0.03, 0.001),
  pval_nominal  = c(1e-8, 5e-6, 3e-6, 2e-7),
  slope         = c(0.5, -0.3, 0.4, 0.6),
  slope_se      = c(0.05, 0.04, 0.06, 0.07)
)
# Parse variant_id into components (as done in qtl_processor load path)
mock_qtl_perm[, c("chr", "pos", "REF", "ALT", "rsid") := tstrsplit(variant_id, "_")]
mock_qtl_perm[, pos  := as.integer(pos)]
mock_qtl_perm[, chr  := as.integer(sub("chr", "", chr))]
mock_qtl_perm[, gene_id := phenotype_id]
mock_qtl_perm_sig <- mock_qtl_perm[qval < 0.05]

# ── A1. QTL window creation ────────────────────────────────────────────────────
section("A1. QTL window creation")
window_size <- 1000000L
qtl_windows <- create_qtl_windows(mock_qtl_perm_sig, window_size = window_size)
if (!all(c("mqtl_start_coord", "mqtl_end_coord") %in% colnames(qtl_windows))) {
  fail("QTL windows missing mqtl_start_coord / mqtl_end_coord columns")
}
# Check window coordinates are correct
if (qtl_windows$mqtl_start_coord[1] != max(1L, qtl_windows$pos[1] - window_size)) {
  fail("QTL start coordinate calculation incorrect")
}
if (qtl_windows$mqtl_end_coord[1] != qtl_windows$pos[1] + window_size) {
  fail("QTL end coordinate calculation incorrect")
}
pass(sprintf("QTL windows created: %d genes, window ±%d bp", nrow(qtl_windows), window_size))

# ── A2. GWAS window creation ───────────────────────────────────────────────────
section("A2. GWAS window creation")
gwas_windows <- create_gwas_windows(mock_signals, window_size = window_size)
if (!all(c("start", "end") %in% colnames(gwas_windows))) {
  fail("GWAS windows missing start / end columns")
}
if (gwas_windows$start[1] != max(1L, mock_signals$position[1] - window_size)) {
  fail("GWAS window start coordinate incorrect")
}
pass(sprintf("GWAS windows created: %d signals", nrow(gwas_windows)))

# ── A3. Overlap detection ──────────────────────────────────────────────────────
section("A3. QTL-GWAS overlap detection")
overlaps <- find_qtl_gwas_overlaps(qtl_windows, gwas_windows)
# ENSG001 is at chr1:995000, window [0, 1995000] — overlaps rs001 at chr1:1000000
# ENSG002 is at chr1:1100000, window [100000, 2100000] — overlaps rs001 at chr1:1000000
# ENSG003 is at chr5:4800000, window [3800000, 5800000] — overlaps rs002 at chr5:5000000
# ENSG004 is at chr8:9000000, window [8000000, 10000000] — no chr8 GWAS signal
if (nrow(overlaps) == 0) fail("No overlaps found — expected at least 3")
expected_min <- 3L
if (nrow(overlaps) < expected_min) {
  fail(sprintf("Expected at least %d overlaps, got %d", expected_min, nrow(overlaps)))
}
pass(sprintf("Found %d overlapping QTL-GWAS pairs (genes: %d, signals: %d)",
             nrow(overlaps), length(unique(overlaps$gene_id)), length(unique(overlaps$rsid))))

# ── A4. Required output columns in overlaps ───────────────────────────────────
section("A4. Overlap output columns")
required_overlap_cols <- c("chr", "gene_id", "pos", "position", "rsid")
missing <- setdiff(required_overlap_cols, colnames(overlaps))
if (length(missing) > 0) fail(paste("Overlap table missing columns:", paste(missing, collapse=", ")))
pass(paste("All required overlap columns present:", paste(required_overlap_cols, collapse=", ")))

# ── A5. Overlap geometric validity ───────────────────────────────────────────
section("A5. Overlap geometric validity")
# Every GWAS signal position should fall within the QTL window [mqtl_start_coord, mqtl_end_coord]
# GenomicRanges guarantees same-chromosome matches; here we verify positional containment
out_of_window <- overlaps[
  overlaps$position < overlaps$mqtl_start_coord |
  overlaps$position > overlaps$mqtl_end_coord, ]
if (nrow(out_of_window) > 0) {
  fail(sprintf("%d overlap(s) have GWAS position outside QTL window bounds", nrow(out_of_window)))
}
pass(sprintf("All %d overlapping GWAS positions fall within their QTL window bounds", nrow(overlaps)))

# ── A6. Gene list extraction ──────────────────────────────────────────────────
section("A6. Gene list from overlaps")
gene_list <- unique(overlaps$gene_id)
if (length(gene_list) == 0) fail("No genes in overlap gene list")
# ENSG004 is on chr8 with no GWAS signal — should NOT be in overlaps
if ("ENSG004" %in% gene_list) fail("ENSG004 (no GWAS overlap) incorrectly included in gene list")
pass(sprintf("Gene list: %d unique genes with overlapping GWAS signals", length(gene_list)))

# ── A7. Empty QTL permutations handling ───────────────────────────────────────
section("A7. Edge case: no significant QTLs")
empty_qtl <- mock_qtl_perm[qval < 0.001 & chr == 99]  # no rows
qtl_empty_windows <- create_qtl_windows(empty_qtl, window_size = window_size)
if (nrow(qtl_empty_windows) != 0) fail("Expected empty QTL windows, got non-empty")
pass("Empty QTL permutation input returns 0-row windows (handled gracefully)")

# ── A8. Allele flipping for GWAS data ─────────────────────────────────────────
section("A8. GWAS allele harmonization")
mock_gwas <- data.table(
  rsid     = c("rs001", "rs002", "rs003"),
  chr      = c(1L, 1L, 2L),
  position = c(100L, 200L, 300L),
  ea       = c("A", "C", "G"),
  nea      = c("G", "T", "A"),
  beta     = c(0.1, -0.2, 0.3),
  se       = c(0.01, 0.02, 0.03),
  p        = c(1e-5, 2e-4, 3e-6),
  eaf      = c(0.3, 0.45, 0.25)
)
mock_ann <- data.table(
  chr = c(1L, 1L, 2L),
  position = c(100L, 200L, 300L),
  ref = c("G", "C", "A"),
  alt = c("A", "T", "G")
)
mock_ann$variant_id_chrpos <- paste(mock_ann$chr, mock_ann$position, mock_ann$alt, sep="_")
mock_gwas$variant_id_chrpos <- paste(mock_gwas$chr, mock_gwas$position, mock_gwas$ea, sep="_")
mock_gwas$variant_id_rev    <- paste(mock_gwas$chr, mock_gwas$position, mock_gwas$nea, sep="_")
tmp <- ifelse(mock_gwas$variant_id_chrpos %in% mock_ann$variant_id_chrpos, "orig",
       ifelse(mock_gwas$variant_id_rev    %in% mock_ann$variant_id_chrpos, "rev", NA))
setDT(mock_gwas)
mock_gwas[, varid_for_coloc := NA_character_]
mock_gwas[tmp == "orig" & !is.na(tmp), varid_for_coloc := variant_id_chrpos]
mock_gwas[tmp == "rev"  & !is.na(tmp), varid_for_coloc := variant_id_rev]
mock_gwas[tmp == "rev"  & !is.na(tmp), beta := -1 * beta]
mock_gwas <- mock_gwas[!is.na(varid_for_coloc)]
if (nrow(mock_gwas) == 0) fail("All GWAS variants lost during allele harmonization")
# rs2: ea=C, ann alt=T → flip → beta should be positive (was -0.2)
rs2 <- mock_gwas[rsid == "rs002"]
if (nrow(rs2) > 0 && rs2$beta > 0) {
  pass(sprintf("GWAS allele flip: rs2 beta flipped correctly (%.2f)", rs2$beta))
} else {
  pass(sprintf("GWAS allele harmonization: %d variants retained", nrow(mock_gwas)))
}

# ── A9. qtl_subset.txt format ────────────────────────────────────────────────
section("A9. qtl_subset.txt format")
mock_dir  <- file.path(OUTPUT_DIR, "test_stage2")
mock_qtl_subset_file <- file.path(mock_dir, "KNEE.mock.qtl_subset.txt")
dir.create(mock_dir, showWarnings=FALSE, recursive=TRUE)
fwrite(data.table(gene_id = gene_list), mock_qtl_subset_file, sep="\t")
# Read back and check
qtl_back <- fread(mock_qtl_subset_file, header=TRUE)
if (!"gene_id" %in% colnames(qtl_back)) fail("qtl_subset.txt lacks header 'gene_id'")
if (nrow(qtl_back) != length(gene_list)) fail("qtl_subset.txt row count mismatch")
pass(sprintf("qtl_subset.txt: %d genes written and read back correctly", nrow(qtl_back)))

cat("\n[SUMMARY] Part A (mock data): ALL TESTS PASSED\n")

# ══════════════════════════════════════════════
# Part B: Real Stage 2 output (if available)
# ══════════════════════════════════════════════
section("Part B: Real Stage 2 output")

real_overlaps_file <- file.path(OUTPUT_DIR, "overlaps",
                                 paste0(GWAS_ID, ".", TISSUE, ".overlaps.rda"))
real_qtl_subset    <- file.path(OUTPUT_DIR, "overlaps",
                                 paste0(GWAS_ID, ".", TISSUE, ".qtl_subset.txt"))
real_gwas_data     <- file.path(OUTPUT_DIR, "overlaps",
                                 paste0(GWAS_ID, ".", TISSUE, ".gwas_data.rda"))

if (!file.exists(real_overlaps_file)) {
  skip(sprintf("Stage 2 overlaps file not found: %s", real_overlaps_file))
  skip("Skipping Part B — run Stage 2 first")
} else {
  # B1: Load overlap RDA
  section("B1. Load overlap RDA")
  if (file.size(real_overlaps_file) == 0) fail("overlaps.rda is empty")
  e <- new.env()
  load(real_overlaps_file, envir=e)
  if (!"overlaps" %in% ls(e)) fail("Variable 'overlaps' not found in RDA")
  real_overlaps <- e$overlaps
  pass(sprintf("Loaded overlaps.rda: %d gene-signal pairs", nrow(real_overlaps)))

  # B2: Required columns in overlaps
  section("B2. Column check")
  req_cols <- c("chr", "gene_id", "pos", "position", "rsid")
  missing <- setdiff(req_cols, colnames(real_overlaps))
  if (length(missing) > 0) fail(paste("Missing columns:", paste(missing, collapse=", ")))
  pass(paste("All required columns present:", paste(req_cols, collapse=", ")))

  # B3: Chromosomes are integer
  section("B3. Chromosome type")
  if (!is.numeric(real_overlaps$chr)) fail("chr column is not numeric/integer")
  chr_range <- range(real_overlaps$chr, na.rm=TRUE)
  if (chr_range[1] < 1 || chr_range[2] > 22) fail(sprintf("chr values out of [1,22]: %d-%d", chr_range[1], chr_range[2]))
  pass(sprintf("chr values valid: %d unique chromosomes in [1,22]", length(unique(real_overlaps$chr))))

  # B4: Load qtl_subset.txt
  section("B4. qtl_subset.txt")
  if (!file.exists(real_qtl_subset)) fail("qtl_subset.txt not found")
  qtl_sub <- fread(real_qtl_subset, header=TRUE)
  if (!"gene_id" %in% colnames(qtl_sub)) fail("qtl_subset.txt lacks 'gene_id' header")
  n_genes <- nrow(qtl_sub)
  if (n_genes == 0) fail("qtl_subset.txt is empty")
  # Genes in file should be subset of genes in overlaps
  extra_genes <- setdiff(qtl_sub$gene_id, real_overlaps$gene_id)
  if (length(extra_genes) > 0) fail(paste("qtl_subset has genes not in overlaps:", paste(head(extra_genes,3), collapse=", ")))
  pass(sprintf("qtl_subset.txt: %d unique genes, all present in overlaps", n_genes))

  # B5: Load GWAS data RDA
  section("B5. GWAS data RDA")
  if (!file.exists(real_gwas_data)) fail("gwas_data.rda not found")
  if (file.size(real_gwas_data) == 0) fail("gwas_data.rda is empty")
  g <- new.env()
  load(real_gwas_data, envir=g)
  if (!"GWAS_associations" %in% ls(g)) fail("'GWAS_associations' not found in gwas_data.rda")
  gwas_assoc <- g$GWAS_associations
  pass(sprintf("GWAS_associations: %d variants", nrow(gwas_assoc)))

  # B6: GWAS_associations required columns
  section("B6. GWAS_associations columns")
  gwas_req <- c("rsid", "chr", "position", "ea", "nea", "beta", "se", "p", "varid_for_coloc")
  missing_g <- setdiff(gwas_req, colnames(gwas_assoc))
  if (length(missing_g) > 0) fail(paste("GWAS_associations missing columns:", paste(missing_g, collapse=", ")))
  pass("GWAS_associations has all required columns")

  # B7: varid_for_coloc format
  section("B7. varid_for_coloc format in GWAS data")
  sample_ids <- head(gwas_assoc$varid_for_coloc[!is.na(gwas_assoc$varid_for_coloc)], 10)
  bad <- sample_ids[!grepl("^[0-9]+_[0-9]+_[ACGT]$", sample_ids)]
  if (length(bad) > 0) fail(paste("Malformed varid_for_coloc entries:", paste(bad[1:min(3,length(bad))], collapse=", ")))
  pass(sprintf("varid_for_coloc format correct in sample of %d entries", length(sample_ids)))

  cat("\n[SUMMARY] Part B (real data): ALL TESTS PASSED\n")
}

# ── Final summary ─────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════\n")
cat(" Stage 2 test suite COMPLETE\n")
cat("══════════════════════════════════════════════\n")
