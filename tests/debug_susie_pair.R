#!/usr/bin/env Rscript
# Debug a single SuSiE gene-signal pair interactively.
# Requires checkpoint saved by 4_run_coloc_susie.R (section 4).
#
# Usage (interactive node, ~5 min startup vs 2h):
#   srun --partition=cpu_p --qos=cpu_normal --mem=32G --time=1:00:00 \
#     /path/to/snakemake-conda-env/bin/Rscript tests/debug_susie_pair.R synovium 1
#
# Or in an interactive R session: source("tests/debug_susie_pair.R")

suppressPackageStartupMessages({
  library(coloc); library(data.table)
})

args   <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "synovium"
pair_i <- if (length(args) >= 2) as.integer(args[2]) else 1L

RESULTS_DIR <- "/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"
CHECKPOINT  <- file.path(RESULTS_DIR, "coloc_susie",
                         sprintf("KNEE.%s.colocSuSiE_results.txt.checkpoint.rds", tissue))
LD_DIR      <- file.path(RESULTS_DIR, "coloc_susie", "LD", sprintf("KNEE.%s", tissue))
PLINK_BIN   <- "/lustre/groups/itg/shared/software/bin/plink"
GWAS_BFILE  <- "/lustre/groups/itg/shared/referenceData/ukbiobank/chip/bgen2plink/CM/chr"
QTL_BFILE   <- c(
  high_grade_cartilage = "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/genotypes/high_grade_cartilage/CleanGenotypes/P21083A_P21083B.indel_snp.recalibrated.autosomal.hardfiltered.samqc.high_grade_cartilage.MAF005CR99.biall_rename",
  low_grade_cartilage  = "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/genotypes/low_grade_cartilage/CleanGenotypes/P21083A_P21083B.indel_snp.recalibrated.autosomal.hardfiltered.samqc.low_grade_cartilage.MAF005CR99.biall_rename",
  synovium             = "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/genotypes/synovium/CleanGenotypes/P21083A_P21083B.indel_snp.recalibrated.autosomal.hardfiltered.samqc.synovium.MAF005CR99.biall_rename"
)[tissue]
QTL_N  <- 133
GWAS_N <- c(cases = 24955, controls = 378169)  # KNEE

# ── Load checkpoint ───────────────────────────────────────────────────────────
stopifnot(file.exists(CHECKPOINT))
cat(sprintf("[%s] Loading checkpoint: %s\n", Sys.time(), CHECKPOINT))
ck <- readRDS(CHECKPOINT)
GWAS_associations <- ck$GWAS_associations
mqtl_df           <- ck$mqtl_df
abf_susie         <- ck$abf_susie
cat(sprintf("  %d GWAS, %d QTL, %d candidates\n",
            nrow(GWAS_associations), nrow(mqtl_df), nrow(abf_susie)))

susie_pairs <- unique(abf_susie[, .(cpg, gwas_signal, coloc_region_chr,
                                    coloc_region_start, coloc_region_end)])
cat(sprintf("  Pairs: %d total, debugging pair %d\n", nrow(susie_pairs), pair_i))
p <- susie_pairs[pair_i]
cat(sprintf("  Gene: %s  Signal: %s  Region: chr%s:%s-%s\n",
            p$cpg, p$gwas_signal, p$coloc_region_chr,
            p$coloc_region_start, p$coloc_region_end))

# ── Pull data for this pair ───────────────────────────────────────────────────
QTL <- mqtl_df[gene_id == p$cpg]
cat(sprintf("  QTL variants for gene: %d\n", nrow(QTL)))
QTL[, SNP := variant_id]
QTL[, A2  := REF]; QTL[, A1 := ALT]
QTL[, A1_freq := af]
QTL[, beta := slope]; QTL[, se := slope_se]; QTL[, pval := pval_nominal]

GWAS_win <- GWAS_associations[chr == p$coloc_region_chr &
                               position >= p$coloc_region_start &
                               position <= p$coloc_region_end]
cat(sprintf("  GWAS variants in window: %d\n", nrow(GWAS_win)))
GWAS_win[, A1_freq := eaf]
GWAS_win[, n := sum(GWAS_N)]
GWAS_win[, s := GWAS_N["cases"] / sum(GWAS_N)]

chr_n <- p$coloc_region_chr

# ── QTL plink LD ─────────────────────────────────────────────────────────────
dir.create(LD_DIR, showWarnings = FALSE, recursive = TRUE)
fname    <- sprintf("%s_%s_mQTL_%s_GWAS_KNEE.txt", p$cpg, p$gwas_signal, tissue)
base_qtl <- file.path(LD_DIR, sub("\\.[^.]*$", "", fname))
f_ref    <- paste0(base_qtl, ".SNP.REF.txt")
f_ext    <- paste0(base_qtl, ".SNP.extract")
f_ld_qtl <- paste0(base_qtl, ".LD")

qtl_for_ld <- unique(QTL[, .(SNP, A2, pos)])
write.table(qtl_for_ld[, .(SNP, A2)], f_ref,
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(qtl_for_ld[, .(SNP)],     f_ext,
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

cmd_qtl <- paste(PLINK_BIN, "--bfile", QTL_BFILE, "--const-fid",
                 "--extract", f_ext, "--threads 1",
                 "--update-ref-allele", f_ref, "2 1",
                 "--r square --out", f_ld_qtl, "--make-just-bim --silent")
cat(sprintf("[%s] Running QTL plink...\n", Sys.time()))
ret <- system(cmd_qtl, ignore.stdout = TRUE, ignore.stderr = TRUE)
cat(sprintf("  plink exit: %d\n", ret))

bim_qtl <- as.data.table(read.table(paste0(f_ld_qtl, ".bim"), header = FALSE,
             col.names = c("chr","rsid","cM","position","a1","a2")))
ld_qtl  <- as.matrix(read.table(paste0(f_ld_qtl, ".ld"), header = FALSE))
colnames(ld_qtl) <- rownames(ld_qtl) <- bim_qtl$rsid
cat(sprintf("  QTL LD: %dx%d\n", nrow(ld_qtl), ncol(ld_qtl)))

# ── GWAS plink LD (two-step) ──────────────────────────────────────────────────
base_gwas  <- file.path(LD_DIR, paste0("GWAS_", sub("\\.[^.]*$", "", fname)))
f_bim_tmp  <- paste0(base_gwas, ".region")
f_ext_gwas <- paste0(base_gwas, ".SNP.extract")
f_ld_gwas  <- paste0(base_gwas, ".LD")
gwas_bfile <- paste0(GWAS_BFILE, chr_n)

cmd_bim <- paste(PLINK_BIN, "--bfile", gwas_bfile, "--const-fid",
                 "--chr", chr_n, "--from-bp", p$coloc_region_start,
                 "--to-bp", p$coloc_region_end,
                 "--threads 1 --make-just-bim --out", f_bim_tmp, "--silent")
cat(sprintf("[%s] Running GWAS region bim...\n", Sys.time()))
system(cmd_bim, ignore.stdout = TRUE, ignore.stderr = TRUE)
bim_full <- as.data.table(read.table(paste0(f_bim_tmp, ".bim"), header = FALSE,
              col.names = c("chr","snp_id","cM","position","a1","a2")))
cat(sprintf("  Region bim: %d variants\n", nrow(bim_full)))

GWAS_win[, match_key_fwd := paste(position, toupper(ea),  toupper(nea), sep = "_")]
GWAS_win[, match_key_rev := paste(position, toupper(nea), toupper(ea),  sep = "_")]
bim_full[,  match_key    := paste(position, toupper(a1),  toupper(a2),  sep = "_")]
GWAS_win[, snp_id := bim_full$snp_id[match(match_key_fwd, bim_full$match_key)]]
GWAS_win[is.na(snp_id),
         snp_id := bim_full$snp_id[match(match_key_rev[is.na(snp_id)], bim_full$match_key)]]
GWAS_matched <- GWAS_win[!is.na(snp_id)]
cat(sprintf("  GWAS matched to bim: %d / %d variants\n", nrow(GWAS_matched), nrow(GWAS_win)))

write.table(GWAS_matched$snp_id, f_ext_gwas,
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
cmd_ld <- paste(PLINK_BIN, "--bfile", gwas_bfile, "--const-fid",
                "--extract", f_ext_gwas, "--threads 1",
                "--r square --out", f_ld_gwas, "--make-just-bim --silent")
cat(sprintf("[%s] Running GWAS LD...\n", Sys.time()))
ret2 <- system(cmd_ld, ignore.stdout = TRUE, ignore.stderr = TRUE)
cat(sprintf("  plink exit: %d\n", ret2))

bim_gwas <- as.data.table(read.table(paste0(f_ld_gwas, ".bim"), header = FALSE,
              col.names = c("chr","snp_id","cM","position","a1","a2")))
ld_gwas  <- as.matrix(read.table(paste0(f_ld_gwas, ".ld"), header = FALSE))
colnames(ld_gwas) <- rownames(ld_gwas) <- bim_gwas$snp_id
cat(sprintf("  GWAS LD: %dx%d\n", nrow(ld_gwas), ncol(ld_gwas)))

# ── Build D1 / D2 ────────────────────────────────────────────────────────────
qtl_snps <- rownames(ld_qtl)
QTL_sub  <- QTL[match(qtl_snps, SNP)]

GWAS_sub <- GWAS_matched[!is.na(A1_freq) & !duplicated(snp_id) &
                          snp_id %in% rownames(ld_gwas)]
ld_gwas  <- ld_gwas[GWAS_sub$snp_id, GWAS_sub$snp_id, drop = FALSE]

cat(sprintf("  D1 snps: %d  D2 snps: %d\n", nrow(QTL_sub), nrow(GWAS_sub)))

D1 <- list(N=QTL_N, MAF=QTL_sub$A1_freq, beta=QTL_sub$beta,
           varbeta=QTL_sub$se^2, type="quant", snp=QTL_sub$SNP,
           pvalue=QTL_sub$pval, position=QTL_sub$pos, LD=ld_qtl)
D1$MAF <- ifelse(D1$MAF <= 0.5, D1$MAF, 1 - D1$MAF)

D2 <- list(N=sum(GWAS_N), MAF=GWAS_sub$A1_freq, beta=GWAS_sub$beta,
           varbeta=GWAS_sub$se^2, type="cc", s=GWAS_N["cases"]/sum(GWAS_N),
           snp=GWAS_sub$snp_id, position=GWAS_sub$position, LD=ld_gwas)
D2$MAF <- ifelse(D2$MAF <= 0.5, D2$MAF, 1 - D2$MAF)

keep1 <- !is.na(D1$MAF) & !is.na(D1$beta) & !is.na(D1$varbeta) & !is.na(D1$position)
keep2 <- !is.na(D2$MAF) & !is.na(D2$beta) & !is.na(D2$varbeta) & !is.na(D2$position)
for (f in c("MAF","beta","varbeta","snp","pvalue","position"))
  if (!is.null(D1[[f]])) D1[[f]] <- D1[[f]][keep1]
D1$LD <- D1$LD[keep1, keep1, drop=FALSE]
for (f in c("MAF","beta","varbeta","snp","position"))
  if (!is.null(D2[[f]])) D2[[f]] <- D2[[f]][keep2]
D2$LD <- D2$LD[keep2, keep2, drop=FALSE]

cat(sprintf("  After NA filter — D1: %d  D2: %d\n", length(D1$snp), length(D2$snp)))
cat(sprintf("  LD dim — D1: %dx%d  D2: %dx%d\n",
            nrow(D1$LD), ncol(D1$LD), nrow(D2$LD), ncol(D2$LD)))
stopifnot(length(D1$snp) == nrow(D1$LD))
stopifnot(length(D2$snp) == nrow(D2$LD))

# ── Run SuSiE ────────────────────────────────────────────────────────────────
cat(sprintf("[%s] Running runsusie D1...\n", Sys.time()))
eqtl_s <- tryCatch(
  runsusie(D1, L=5, coverage=0.95, repeat_until_convergence=FALSE),
  error = function(e) { cat("  D1 L=5 error:", e$message, "\n"); NULL })

cat(sprintf("[%s] Running runsusie D2...\n", Sys.time()))
gwas_s <- tryCatch(
  runsusie(D2, L=5, coverage=0.95, repeat_until_convergence=FALSE),
  error = function(e) { cat("  D2 L=5 error:", e$message, "\n"); NULL })

cat(sprintf("  eqtl_s class: %s  cs: %s\n",
            class(eqtl_s), ifelse(inherits(eqtl_s,"susie"), length(eqtl_s$sets$cs), "N/A")))
cat(sprintf("  gwas_s class: %s  cs: %s\n",
            class(gwas_s), ifelse(inherits(gwas_s,"susie"), length(gwas_s$sets$cs), "N/A")))

if (inherits(eqtl_s,"susie") && inherits(gwas_s,"susie")) {
  cat(sprintf("[%s] Running coloc.susie...\n", Sys.time()))
  res <- tryCatch(coloc.susie(gwas_s, eqtl_s, p12=1e-5),
                  error = function(e) { cat("  coloc.susie error:", e$message, "\n"); NULL })
  if (!is.null(res)) print(res$summary)
}
cat(sprintf("[%s] Done.\n", Sys.time()))
