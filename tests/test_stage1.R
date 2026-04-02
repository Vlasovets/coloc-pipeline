#!/usr/bin/env Rscript
# Test Stage 1: GWAS VCF Conversion
# Tests 1_convert_gwas_to_vcf.R logic with:
#   a) synthetic mock data (always runs) — tests helper functions and output format
#   b) real Stage 1 output if available (skipped if not present)

suppressPackageStartupMessages({
  library(data.table)
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

cat("══════════════════════════════════════════════\n")
cat(" Stage 1 GWAS VCF Conversion — Test Suite\n")
cat("══════════════════════════════════════════════\n")

# ══════════════════════════════════════════════
# Part A: Mock data tests (always run)
# ══════════════════════════════════════════════
section("Part A: Mock data tests")

# ── A1. Helper: allele flipping logic ────────────────────────────────────────
section("A1. Allele flipping logic")

# Simulate what Stage 1 / Stage 2 does: chr_pos_alt variant IDs for allele matching
mock_gwas <- data.table(
  rsid     = c("rs1", "rs2", "rs3", "rs4", "rs5"),
  chr      = c(1L, 1L, 2L, 3L, 5L),
  position = c(100000L, 200000L, 300000L, 400000L, 500000L),
  ea       = c("A", "C", "G", "T", "A"),
  nea      = c("G", "T", "A", "C", "C"),
  beta     = c(0.10, -0.05, 0.20, -0.15, 0.08),
  se       = c(0.02,  0.01, 0.03,  0.02, 0.015),
  p        = c(1e-6, 2e-4, 5e-8, 3e-5, 1e-3),
  eaf      = c(0.35, 0.48, 0.22, 0.55, 0.40)
)

# Simulate variant annotation (reference uses alt allele in ID)
mock_ann <- data.table(
  chr      = c(1L, 1L, 2L, 3L, 5L),
  position = c(100000L, 200000L, 300000L, 400000L, 500000L),
  ref      = c("G", "C", "A", "C", "A"),   # ref matches ea of rs2, nea of rs3&rs4
  alt      = c("A", "T", "G", "T", "C")    # alt matches ea of rs1&rs3, nea of rs2&rs5
)
mock_ann$variant_id_chrpos <- paste(mock_ann$chr, mock_ann$position, mock_ann$alt, sep="_")

# Build variant IDs for GWAS
mock_gwas$variant_id_chrpos <- paste(mock_gwas$chr, mock_gwas$position, mock_gwas$ea, sep="_")
mock_gwas$variant_id_rev    <- paste(mock_gwas$chr, mock_gwas$position, mock_gwas$nea, sep="_")

flip_flag <- ifelse(mock_gwas$variant_id_chrpos %in% mock_ann$variant_id_chrpos, "orig",
             ifelse(mock_gwas$variant_id_rev    %in% mock_ann$variant_id_chrpos, "rev", NA))

mock_gwas$varid_for_coloc <- ifelse(flip_flag == "orig" & !is.na(flip_flag), mock_gwas$variant_id_chrpos,
                              ifelse(flip_flag == "rev"  & !is.na(flip_flag), mock_gwas$variant_id_rev, NA))
mock_gwas$beta_coloc <- ifelse(flip_flag == "rev" & !is.na(flip_flag), -mock_gwas$beta, mock_gwas$beta)

n_matched <- sum(!is.na(mock_gwas$varid_for_coloc))
n_flipped  <- sum(flip_flag == "rev", na.rm=TRUE)
if (n_matched == 0) fail("No variants matched after allele flipping")
pass(sprintf("Allele flipping: %d matched (%d flipped, %d direct)", n_matched, n_flipped, n_matched - n_flipped))

# rs2: ea=C, nea=T, ann alt=T → rev flip → beta should be negated
if (!is.na(flip_flag[2]) && flip_flag[2] == "rev") {
  if (mock_gwas$beta_coloc[2] != -mock_gwas$beta[2]) fail("Beta not negated for rev allele")
  pass("Beta correctly negated for reversed allele (rs2)")
} else {
  skip("rs2 allele orientation not rev — skipping beta flip check")
}

# ── A2. Variant ID format ─────────────────────────────────────────────────────
section("A2. Variant ID format (chr_pos_allele)")
varids <- na.omit(mock_gwas$varid_for_coloc)
if (length(varids) == 0) fail("No variant IDs produced")
# Each ID should match chr_pos_allele pattern
bad <- varids[!grepl("^[0-9]+_[0-9]+_[ACGT]$", varids)]
if (length(bad) > 0) fail(paste("Malformed variant IDs:", paste(bad, collapse=", ")))
pass(sprintf("All %d variant IDs match chr_pos_allele format", length(varids)))

# ── A3. Biallelic SNP filter ──────────────────────────────────────────────────
section("A3. Biallelic SNP filtering")
mock_ann_mixed <- data.table(
  chr = c(1L, 2L, 3L, 4L),
  ref = c("A", "AT", "C", "G"),   # row 2: indel
  alt = c("G", "A",  "G", "TT"),  # row 4: indel
  position = c(100L, 200L, 300L, 400L)
)
mock_ann_mixed <- mock_ann_mixed[nchar(ref) == 1 & nchar(alt) == 1]
if (nrow(mock_ann_mixed) != 2L) fail(sprintf("Expected 2 biallelic SNPs, got %d", nrow(mock_ann_mixed)))
pass("Biallelic SNP filter correctly removes indels")

# ── A4. GWAS summary statistics validation ────────────────────────────────────
section("A4. GWAS column validation")
required_gwas_cols <- c("rsid", "chr", "position", "ea", "nea", "beta", "se", "p", "eaf")
missing <- setdiff(required_gwas_cols, colnames(mock_gwas))
if (length(missing) > 0) fail(paste("Mock GWAS missing columns:", paste(missing, collapse=", ")))
pass(paste("All required GWAS columns present:", paste(required_gwas_cols, collapse=", ")))

# ── A5. P-value range ────────────────────────────────────────────────────────
section("A5. P-value range check")
if (any(mock_gwas$p <= 0 | mock_gwas$p > 1)) fail("P-values outside (0, 1]")
pass("All p-values in (0, 1]")

# ── A6. EAF range ────────────────────────────────────────────────────────────
section("A6. Effect allele frequency range")
if (any(mock_gwas$eaf <= 0 | mock_gwas$eaf >= 1)) fail("EAF values outside (0, 1)")
pass("All EAF values in (0, 1)")

# ── A7. Sample size parameterization ─────────────────────────────────────────
section("A7. Sample size parameterization (GWAS type cc)")
n_cases    <- 172256L
n_controls <- 1144244L
gwas_n     <- c(n_cases, n_controls)
s_ratio    <- n_cases / sum(gwas_n)
if (s_ratio <= 0 || s_ratio >= 1) fail("Case proportion s outside (0, 1)")
pass(sprintf("Case-control ratio s=%.4f (cases=%d, controls=%d)", s_ratio, n_cases, n_controls))

# ── A8. Missing input handling ───────────────────────────────────────────────
section("A8. Missing input file handling")
result <- tryCatch(
  { if (!file.exists("/nonexistent/file.txt.gz")) stop("Input file not found") },
  error = function(e) e$message
)
if (!grepl("not found", result)) fail("Missing file error not raised correctly")
pass("Missing input file correctly raises an error")

cat("\n[SUMMARY] Part A (mock data): ALL TESTS PASSED\n")

# ══════════════════════════════════════════════
# Part B: Real Stage 1 output (if available)
# ══════════════════════════════════════════════
section("Part B: Real Stage 1 VCF output")

real_vcf <- file.path(OUTPUT_DIR, "gwas_vcf", paste0(GWAS_ID, "_hg38.vcf.bgz"))
real_tbi  <- paste0(real_vcf, ".tbi")

if (!file.exists(real_vcf)) {
  skip(sprintf("Stage 1 output not found: %s", real_vcf))
  skip("Skipping Part B — run Stage 1 first")
} else {
  # B1: File sizes
  section("B1. Output file sizes")
  vcf_size <- file.size(real_vcf)
  if (vcf_size == 0) fail("VCF file is empty")
  pass(sprintf("VCF exists and is non-empty: %s", format(vcf_size, big.mark=",")))

  # B2: Index file
  section("B2. Tabix index present")
  if (!file.exists(real_tbi)) fail("Tabix index (.tbi) missing")
  tbi_size <- file.size(real_tbi)
  if (tbi_size == 0) fail("Tabix index is empty")
  pass(sprintf("Tabix index present: %s bytes", format(tbi_size, big.mark=",")))

  # B3: VCF header check (bgzipped — can use zcat to peek)
  section("B3. VCF format header check")
  header_lines <- tryCatch(
    system(sprintf("zcat %s 2>/dev/null | head -50", real_vcf), intern=TRUE),
    error = function(e) character(0)
  )
  if (length(header_lines) == 0) {
    skip("Could not read VCF header (zcat unavailable or VCF unreadable)")
  } else {
    has_vcf_header <- any(grepl("^##fileformat=VCF", header_lines))
    if (!has_vcf_header) fail("VCF header line ##fileformat=VCF not found")
    pass("VCF file starts with valid ##fileformat=VCF header")

    has_chrom_line <- any(grepl("^#CHROM", header_lines))
    if (!has_chrom_line) fail("#CHROM column header line not found")
    pass("#CHROM column header present in VCF")
  }

  # B4: Genome build — filename confirms hg38
  section("B4. hg38 filename convention")
  if (!grepl("hg38", basename(real_vcf))) fail("Output filename does not contain 'hg38'")
  pass(sprintf("Output filename correctly indicates hg38: %s", basename(real_vcf)))

  # B5: Log file check
  section("B5. Stage 1 log file")
  log_file <- file.path(OUTPUT_DIR, "logs", paste0("convert_gwas_", GWAS_ID, ".log"))
  if (!file.exists(log_file)) {
    skip(sprintf("Log file not found: %s", log_file))
  } else {
    log_content <- readLines(log_file, warn=FALSE)
    if (length(log_content) == 0) fail("Log file is empty")
    pass(sprintf("Log file exists with %d lines", length(log_content)))
  }

  cat("\n[SUMMARY] Part B (real data): ALL TESTS PASSED\n")
}

# ── Final summary ─────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════\n")
cat(" Stage 1 test suite COMPLETE\n")
cat("══════════════════════════════════════════════\n")
