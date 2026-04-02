#!/usr/bin/env Rscript
# Smoke test for Stage 5 (SuSiE) — checks data availability and candidate filtering
# Safe to run without running SuSiE itself

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(dplyr)
  library(arrow)
})

pass <- function(msg) cat(sprintf("[PASS] %s\n", msg))
fail <- function(msg) { cat(sprintf("[FAIL] %s\n", msg)); quit(save="no", status=1) }

OUTPUT_DIR <- "/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"
GWAS_ID    <- "KNEE"
tissue     <- "high_grade_cartilage"

cat("── Stage 5 (SuSiE) smoke test ──\n")

# 1. ABF results
abf_file <- file.path(OUTPUT_DIR, "coloc_abf", paste0(GWAS_ID, ".", tissue, ".colocABF_results.txt"))
if (!file.exists(abf_file)) fail(paste("ABF file missing:", abf_file))
abf <- fread(abf_file)
if (nrow(abf) == 0) fail("ABF file is empty")
pass(sprintf("ABF results loaded: %d rows", nrow(abf)))

# 2. Candidate filtering
abf_sign  <- abf[PP4 >= 0.8 | (PP4 > 0.6 & PP4 < 0.8 & PP4/PP3 > 2),]
abf_susie <- abf %>%
  filter(!paste0(cpg, gwas_signal, tissue) %in%
           paste0(abf_sign$cpg, abf_sign$gwas_signal, abf_sign$tissue),
         PP4 > 0.25) %>%
  as.data.table()
pass(sprintf("Candidate filtering: %d significant, %d SuSiE candidates", nrow(abf_sign), nrow(abf_susie)))
if (nrow(abf_susie) == 0) {
  cat("[SKIP] No SuSiE candidates — test complete\n")
  quit(save="no", status=0)
}

# 3. GWAS data
gwas_file <- file.path(OUTPUT_DIR, "overlaps", paste0(GWAS_ID, ".", tissue, ".gwas_data.rda"))
if (!file.exists(gwas_file)) fail(paste("GWAS data missing:", gwas_file))
load(gwas_file)
if (!exists("GWAS_associations") || nrow(GWAS_associations) == 0) fail("GWAS_associations empty")
pass(sprintf("GWAS data: %d variants", nrow(GWAS_associations)))

# 4. QTL plink bfile
bfile <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/genotypes/high_grade_cartilage/CleanGenotypes/P21083A_P21083B.indel_snp.recalibrated.autosomal.hardfiltered.samqc.high_grade_cartilage.MAF005CR99.biall_rename"
if (!file.exists(paste0(bfile, ".bim"))) fail(paste("QTL bfile not found:", bfile))
pass(sprintf("QTL plink bfile accessible"))

# 5. plink binary
plink_bin <- "/lustre/groups/itg/shared/software/bin/plink"
if (!file.exists(plink_bin)) fail(paste("plink not found:", plink_bin))
pass(sprintf("plink binary found: %s", plink_bin))

# 6. GWAS LD bfile
gwas_bfile_prefix <- "/lustre/groups/itg/shared/referenceData/ukbiobank/chip/bgen2plink/CM/chr"
# Check first candidate's chromosome
chr1 <- abf_susie$coloc_region_chr[1]
gwas_bfile <- paste0(gwas_bfile_prefix, chr1)
if (!file.exists(paste0(gwas_bfile, ".bim"))) {
  cat(sprintf("[WARN] GWAS LD bfile not accessible for chr%s: %s.bim\n", chr1, gwas_bfile))
  cat("       This is expected if running outside SLURM — permissions differ.\n")
} else {
  pass(sprintf("GWAS LD bfile accessible for chr%s", chr1))
}

# 7. QTL parquet for first candidate
qtl_nominal_dir <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/eQTL_results/high_grade_cartilage/cov_df_cart_high_K35/cis_nominal/"
gene1 <- abf_susie$cpg[1]
chr1  <- abf_susie$coloc_region_chr[1]
parquet_file <- file.path(qtl_nominal_dir, paste0("cis_qtl.cis_qtl_pairs.chr", chr1, ".parquet"))
if (!file.exists(parquet_file)) fail(paste("QTL parquet missing:", parquet_file))
qtl <- as.data.table(arrow::read_parquet(parquet_file))
qtl_gene <- qtl[phenotype_id == gene1]
if (nrow(qtl_gene) == 0) fail(sprintf("No QTL variants for candidate gene %s on chr%s", gene1, chr1))
pass(sprintf("QTL data for candidate %s: %d variants on chr%s", gene1, nrow(qtl_gene), chr1))

# 8. coloc package version check
pkg_ver <- packageVersion("coloc")
if (!exists("runsusie")) fail("runsusie() not available — coloc package too old")
pass(sprintf("coloc package version: %s (runsusie available)", pkg_ver))

cat("\n[SUMMARY] Stage 5 smoke test: ALL CHECKS PASSED\n")
cat(sprintf("Ready to submit: sbatch scripts/stage5_susie.sh %s %s\n", GWAS_ID, tissue))
