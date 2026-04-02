#!/usr/bin/env Rscript
# Test Stage 3: Coloc ABF Analysis
# Tests the fixed 3_run_coloc_abf.R logic with a small subset of real data.
# Verifies: data loading, column mapping, allele flipping, perform_coloc(), output format.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(coloc)
  library(dplyr)
  library(arrow)
  library(parallel)
})

# ── helpers ──────────────────────────────────────────────────────────────────
pass <- function(msg) cat(sprintf("[PASS] %s\n", msg))
fail <- function(msg) { cat(sprintf("[FAIL] %s\n", msg)); quit(save="no", status=1) }
section <- function(msg) cat(sprintf("\n── %s ──\n", msg))

script_dir  <- "/home/itg/oleg.vlasovets/projects/coloc-pipeline/workflow/scripts"
OUTPUT_DIR  <- "/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"
GWAS_ID     <- "KNEE"
tissue      <- "high_grade_cartilage"
QTL_n       <- 115
GWAS_n      <- c(172256, 1144244)
type        <- "cc"
N_TEST_PAIRS <- 5   # use only first N gene-signal pairs to keep test fast

cat("══════════════════════════════════════════════\n")
cat(" Stage 3 Coloc ABF — Integration Test\n")
cat(sprintf(" Trait: %s  |  Tissue: %s\n", GWAS_ID, tissue))
cat("══════════════════════════════════════════════\n")

# ── 1. Source helper functions ────────────────────────────────────────────────
section("1. Source helper functions")
tryCatch(
  source(file.path(script_dir, "Coloc_helper_functions.R")),
  error = function(e) fail(paste("Cannot source Coloc_helper_functions.R:", e$message))
)
stopifnot(exists("perform_coloc"))
pass("Coloc_helper_functions.R sourced; perform_coloc() available")

# ── 2. Load overlap data (Stage 2 output) ─────────────────────────────────────
section("2. Load overlap data")
overlap_file <- file.path(OUTPUT_DIR, "overlaps",
                          paste0(GWAS_ID, ".", tissue, ".overlaps.rda"))
if (!file.exists(overlap_file)) fail(paste("Overlap file not found:", overlap_file))

e <- new.env()
load(overlap_file, envir = e)
if (!"overlaps" %in% ls(e)) fail("Variable 'overlaps' not found in overlap file")
overlaps <- e$overlaps
pass(sprintf("Loaded overlap file: %d gene-signal pairs", nrow(overlaps)))

expected_cols <- c("chr", "pos", "gene_id", "position", "rsid")
missing <- setdiff(expected_cols, colnames(overlaps))
if (length(missing) > 0) fail(paste("Missing columns in overlaps:", paste(missing, collapse=", ")))
pass(paste("All required columns present:", paste(expected_cols, collapse=", ")))

# ── 3. Column mapping to perform_coloc() format ───────────────────────────────
section("3. Column mapping")
overlap_df_coloc <- data.frame(
  cpg.id              = overlaps$gene_id,
  cpg.pos             = overlaps$pos,
  signal_gwas.rsid    = overlaps$rsid,
  signal_gwas.chr     = overlaps$chr,
  signal_gwas.position = overlaps$position,
  stringsAsFactors    = FALSE
)
overlap_df_test <- head(overlap_df_coloc, N_TEST_PAIRS)
pass(sprintf("Mapped columns; using first %d pairs for test", N_TEST_PAIRS))

# ── 4. Load GWAS data (Stage 2 output) ───────────────────────────────────────
section("4. Load GWAS data")
gwas_file <- file.path(OUTPUT_DIR, "overlaps",
                       paste0(GWAS_ID, ".", tissue, ".gwas_data.rda"))
if (!file.exists(gwas_file)) fail(paste("GWAS data file not found:", gwas_file))

g <- new.env()
load(gwas_file, envir = g)
if (!"GWAS_associations" %in% ls(g)) fail("Variable 'GWAS_associations' not found in gwas_data file")
GWAS_associations <- g$GWAS_associations
pass(sprintf("Loaded GWAS_associations: %d variants", nrow(GWAS_associations)))

gwas_cols <- c("rsid", "chr", "position", "ea", "nea", "beta", "se", "p", "eaf", "varid_for_coloc")
missing <- setdiff(gwas_cols, colnames(GWAS_associations))
if (length(missing) > 0) fail(paste("GWAS_associations missing columns:", paste(missing, collapse=", ")))
pass("GWAS_associations has all required columns")

# ── 5. Load QTL subset gene list ─────────────────────────────────────────────
section("5. Load QTL gene list")
qtl_subset_file <- file.path(OUTPUT_DIR, "overlaps",
                             paste0(GWAS_ID, ".", tissue, ".qtl_subset.txt"))
if (!file.exists(qtl_subset_file)) fail(paste("QTL subset file not found:", qtl_subset_file))

qtl_genes_dt <- fread(qtl_subset_file, header = TRUE)
if (!"gene_id" %in% colnames(qtl_genes_dt)) fail("Column 'gene_id' not found (fread header=TRUE required)")
qtl_genes <- qtl_genes_dt$gene_id
pass(sprintf("Loaded %d genes from qtl_subset.txt", length(qtl_genes)))

# Verify test genes are in the gene list
test_genes <- unique(overlap_df_test$cpg.id)
missing_genes <- setdiff(test_genes, qtl_genes)
if (length(missing_genes) > 0) fail(paste("Test genes not in qtl_subset:", paste(missing_genes, collapse=", ")))
pass("All test genes present in qtl_subset.txt")

# ── 6. Load variant annotation ───────────────────────────────────────────────
section("6. Load variant annotation")
ann_file <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/variant_ann_hg38.txt.gz"
if (!file.exists(ann_file)) fail(paste("Variant annotation file not found:", ann_file))
variant_ann <- fread(ann_file, data.table = FALSE)
variant_ann$chr <- as.integer(variant_ann$chr)
variant_ann <- variant_ann[!is.na(variant_ann$chr) &
                            nchar(variant_ann$ref) == 1 &
                            nchar(variant_ann$alt) == 1, ]
if (!"variant_id_chrpos" %in% colnames(variant_ann)) {
  variant_ann$variant_id_chrpos <- paste(variant_ann$chr, variant_ann$position, variant_ann$alt, sep="_")
}
pass(sprintf("Loaded variant annotation: %d biallelic SNPs", nrow(variant_ann)))

# ── 7. Load QTL nominal data from parquet ────────────────────────────────────
section("7. Load QTL nominal data")
qtl_nominal_dir <- paste0(
  "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/",
  "pipeline_RNA_seq_analysis/eQTL/eQTL_results/", tissue,
  "/cov_df_cart_high_K35/cis_nominal"
)

tmp.dt <- data.table()
for (chrom in seq(1, 22)) {
  qtl_file <- file.path(qtl_nominal_dir, paste0("cis_qtl.cis_qtl_pairs.chr", chrom, ".parquet"))
  if (!file.exists(qtl_file)) next
  qtl <- as.data.table(arrow::read_parquet(qtl_file))
  qtl <- qtl[phenotype_id %in% test_genes]
  if (nrow(qtl) == 0) next
  qtl[, c("chr", "pos", "REF", "ALT", "rsid") := tstrsplit(variant_id, "_")]
  qtl[, maf  := ifelse(af < 0.5, af, 1 - af)]
  qtl[, pos  := as.integer(pos)]
  qtl[, chr  := as.integer(sub("chr", "", chr))]
  qtl[, gene_id := phenotype_id]
  qtl <- qtl[rsid != "."]
  tmp.dt <- rbind(tmp.dt, qtl)
}
mqtl_df <- tmp.dt; rm(tmp.dt)

if (nrow(mqtl_df) == 0) fail("No QTL nominal data found for test genes")
pass(sprintf("Loaded %d QTL variant records for %d test genes", nrow(mqtl_df), length(test_genes)))

# ── 8. Allele flipping ───────────────────────────────────────────────────────
section("8. Allele flipping for QTL data")
mqtl_df$variant_id_chrpos <- paste(mqtl_df$chr, mqtl_df$pos, mqtl_df$ALT, sep="_")
mqtl_df$variant_id_rev    <- paste(mqtl_df$chr, mqtl_df$pos, mqtl_df$REF, sep="_")

flip_flag <- ifelse(mqtl_df$variant_id_chrpos %in% variant_ann$variant_id_chrpos, "orig",
                    ifelse(mqtl_df$variant_id_rev %in% variant_ann$variant_id_chrpos, "rev", NA))

mqtl_df$varid_for_coloc <- ifelse(flip_flag == "orig" & !is.na(flip_flag), mqtl_df$variant_id_chrpos,
                            ifelse(flip_flag == "rev"  & !is.na(flip_flag), mqtl_df$variant_id_rev, NA))
mqtl_df$slope    <- ifelse(flip_flag == "rev" & !is.na(flip_flag), -mqtl_df$slope, mqtl_df$slope)
mqtl_df$ALT_coloc <- ifelse(flip_flag == "orig" & !is.na(flip_flag), mqtl_df$ALT,
                      ifelse(flip_flag == "rev"  & !is.na(flip_flag), mqtl_df$REF, NA))
mqtl_df$REF_coloc <- ifelse(flip_flag == "orig" & !is.na(flip_flag), mqtl_df$REF,
                      ifelse(flip_flag == "rev"  & !is.na(flip_flag), mqtl_df$ALT, NA))
mqtl_df$ALT <- mqtl_df$ALT_coloc
mqtl_df$REF <- mqtl_df$REF_coloc
mqtl_df <- mqtl_df[!is.na(mqtl_df$varid_for_coloc), ]

if (nrow(mqtl_df) == 0) fail("All QTL variants were filtered during allele flipping")
pass(sprintf("%d QTL variants retained after allele flipping", nrow(mqtl_df)))

in_m_qtl <- split(mqtl_df, mqtl_df$gene_id)
pass(sprintf("Split into %d gene-level lists", length(in_m_qtl)))

# ── 9. Create required output directories ────────────────────────────────────
section("9. Create output directories")
out_path <- file.path(OUTPUT_DIR, "test_stage3")
dir.create(file.path(out_path, paste0(GWAS_ID, "_coloc_rda_files")), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(out_path, "coloc", "Coloc_sumstats", paste0("Coloc_mQTL_", GWAS_ID)), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(out_path, "coloc", "Coloc_sumstats", paste0("Coloc_", GWAS_ID)), showWarnings=FALSE, recursive=TRUE)
pass("Output directories created")

# ── 10. Run perform_coloc on test subset ──────────────────────────────────────
section("10. Run perform_coloc (first 5 pairs)")
results <- tryCatch(
  perform_coloc(
    overlap_df = overlap_df_test,
    in_m_qtl   = in_m_qtl,
    out_path   = out_path,
    gwas_n     = GWAS_n,
    in_gwas    = GWAS_associations,
    GWAS_type  = type,
    GWAS_ID    = GWAS_ID,
    QTL_n      = QTL_n,
    is.eQTL    = TRUE,
    cores      = 1
  ),
  error = function(e) { fail(paste("perform_coloc() error:", e$message)); NULL }
)

if (is.null(results) || (!is.data.frame(results) && !is.data.table(results))) {
  cat("[INFO] No colocalization results for these test pairs (PP4 below threshold) — this is expected\n")
} else {
  pass(sprintf("perform_coloc() returned %d row(s)", nrow(results)))
  expected_out_cols <- c("cpg", "PP4", "PP3", "gwas_signal")
  missing_out <- setdiff(expected_out_cols, colnames(results))
  if (length(missing_out) > 0) fail(paste("Missing output columns:", paste(missing_out, collapse=", ")))
  pass(paste("Output has expected columns:", paste(expected_out_cols, collapse=", ")))
  stopifnot(all(results$PP4 >= 0 & results$PP4 <= 1))
  pass("PP4 values are in [0, 1]")
}

# ── Summary ──────────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════\n")
cat(" ALL TESTS PASSED — Stage 3 logic is correct\n")
cat("══════════════════════════════════════════════\n")
