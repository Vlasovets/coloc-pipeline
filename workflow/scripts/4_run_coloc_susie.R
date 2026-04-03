#!/usr/bin/env Rscript

#' Stage 5: Coloc SuSiE Fine-mapping
#'
#' Runs coloc.susie for ambiguous colocalization cases (PP4 > 0.25, not already
#' significant from ABF). Re-extracts QTL data per-tissue from parquet files,
#' computes per-locus LD matrices with plink, then runs runsusie + coloc.susie.
#'
#' Adapted from src/4_run_coloc_susie.R for Snakemake modular pipeline.

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(dplyr)
  library(arrow)
})

source(file.path(snakemake@scriptdir, "coloc_helpers.R"))

###############################################################################
# Pipeline-compatible LD helpers (portable: accept ld_dir and plink_bin args)
###############################################################################

#' Compute LD matrix from plink bfile for a set of variants
#'
#' @param snp_df  data.table with columns SNP (rsid), A2 (ref allele), pos
#' @param chr     chromosome number
#' @param fname   base filename for temp plink outputs
#' @param bfile   plink bfile prefix (no extension)
#' @param ld_dir  output directory for plink LD files
#' @param plink_bin  path to plink binary
#' @return named LD matrix (SNP x SNP correlation) or NULL on failure
compute_ld_matrix <- function(snp_df, chr, fname, bfile, ld_dir, plink_bin) {
  dir.create(ld_dir, showWarnings = FALSE, recursive = TRUE)
  base   <- file.path(ld_dir, sub("\\.[^.]*$", "", fname))
  f_ref  <- paste0(base, ".SNP.REF.txt")
  f_ext  <- paste0(base, ".SNP.extract")
  f_ld   <- paste0(base, ".LD")

  write.table(snp_df[, c("SNP", "A2")], f_ref,
              col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(snp_df[, "SNP"],           f_ext,
              col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

  cmd <- paste(plink_bin,
    "--bfile", bfile,
    "--const-fid",
    "--extract", f_ext,
    "--threads 1",
    "--update-ref-allele", f_ref, "2 1",
    "--r square",
    "--out", f_ld,
    "--make-just-bim",
    "--silent"
  )
  ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (ret != 0 || !file.exists(paste0(f_ld, ".ld"))) return(NULL)

  bim <- tryCatch(
    read.table(paste0(f_ld, ".bim"), header = FALSE,
               col.names = c("chr", "rsid", "cM", "position", "a1", "a2")),
    error = function(e) NULL
  )
  if (is.null(bim)) return(NULL)

  ld <- tryCatch(
    as.matrix(read.table(paste0(f_ld, ".ld"), header = FALSE)),
    error = function(e) NULL
  )
  if (is.null(ld)) return(NULL)

  colnames(ld) <- bim$rsid
  rownames(ld) <- bim$rsid

  keep <- rownames(ld) %in% snp_df$SNP
  ld[keep, keep, drop = FALSE]
}

###############################################################################
# Snakemake parameters
###############################################################################
abf_results_file  <- snakemake@input$abf_results
gwas_data_file    <- snakemake@input$gwas_data
qtl_subset_file   <- snakemake@input$qtl_subset
output_file       <- snakemake@output$results
ld_output_dir     <- snakemake@output$ld_dir

GWAS_ID           <- snakemake@params$gwas_id
tissue            <- snakemake@params$tissue
QTL_n             <- snakemake@params$qtl_n
GWAS_n            <- snakemake@params$gwas_n      # c(cases, controls) or single N
qtl_nominal_dir   <- snakemake@params$qtl_nominal_dir
qtl_bfile         <- snakemake@params$qtl_bfile   # plink bfile prefix (tissue-specific)
gwas_bfile_prefix <- snakemake@params$gwas_bfile_prefix  # e.g. /path/to/CM/chr (chr added per locus)
pp4_threshold     <- snakemake@params$pp4_threshold    # default 0.25
pp4_pp3_ratio     <- snakemake@params$pp4_pp3_ratio    # default 2.0
cores             <- snakemake@threads

plink_bin         <- snakemake@params$plink_bin

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, append = TRUE)
sink(log_con, append = TRUE, type = "message")

cat(sprintf("[%s] Stage 5: Coloc SuSiE — %s / %s\n", Sys.time(), GWAS_ID, tissue))

###############################################################################
# 1. Load ABF results and filter to ambiguous candidates
###############################################################################
cat(sprintf("[%s] Loading ABF results from %s\n", Sys.time(), abf_results_file))
abf <- fread(abf_results_file)
if (nrow(abf) == 0) {
  cat(sprintf("[%s] No ABF results — skipping SuSiE\n", Sys.time()))
  fwrite(data.table(), output_file, sep = "\t")
  dir.create(ld_output_dir, showWarnings = FALSE, recursive = TRUE)
  quit(save = "no", status = 0)
}

type <- ifelse(length(GWAS_n) == 1, "quant", "cc")

# Already-significant: PP4 >= 0.8 OR (PP4 > 0.6 AND PP4/PP3 > pp4_pp3_ratio)
abf_sign <- abf[PP4 >= 0.8 | (PP4 > 0.6 & PP4 < 0.8 & PP4 / PP3 > pp4_pp3_ratio), ]

# Ambiguous candidates: PP4 > threshold AND not already significant for same gene-signal
# Note: per-tissue ABF files do not have a 'tissue' column — match on cpg + gwas_signal only
abf_susie <- abf %>%
  filter(
    !paste0(cpg, gwas_signal) %in% paste0(abf_sign$cpg, abf_sign$gwas_signal),
    PP4 > pp4_threshold
  ) %>%
  as.data.table()

cat(sprintf("[%s] ABF total: %d | significant: %d | SuSiE candidates: %d\n",
            Sys.time(), nrow(abf), nrow(abf_sign), nrow(abf_susie)))

if (nrow(abf_susie) == 0) {
  cat(sprintf("[%s] No ambiguous candidates — skipping SuSiE\n", Sys.time()))
  fwrite(data.table(), output_file, sep = "\t")
  dir.create(ld_output_dir, showWarnings = FALSE, recursive = TRUE)
  quit(save = "no", status = 0)
}

###############################################################################
# 2. Load GWAS data (from Stage 2 gwas_data.rda)
###############################################################################
cat(sprintf("[%s] Loading GWAS data from %s\n", Sys.time(), gwas_data_file))
load(gwas_data_file)  # loads GWAS_associations

###############################################################################
# 3. Re-extract QTL data from parquet for this tissue
###############################################################################
cat(sprintf("[%s] Loading QTL data for %s from %s\n", Sys.time(), tissue, qtl_nominal_dir))
qtl_genes_all <- fread(qtl_subset_file, header = TRUE)$gene_id
susie_genes  <- unique(abf_susie$cpg)
qtl_genes    <- intersect(qtl_genes_all, susie_genes)

mqtl_df <- data.table()
for (chr_n in seq(1, 22)) {
  qtl_file <- file.path(qtl_nominal_dir,
                        paste0("cis_qtl.cis_qtl_pairs.chr", chr_n, ".parquet"))
  if (!file.exists(qtl_file)) next
  qtl <- as.data.table(arrow::read_parquet(qtl_file))
  qtl <- qtl[phenotype_id %in% qtl_genes]
  if (nrow(qtl) == 0) next
  qtl[, c("chr", "pos", "REF", "ALT", "rsid") := tstrsplit(variant_id, "_")]
  qtl[, maf  := ifelse(af < 0.5, af, 1 - af)]
  qtl[, pos  := as.integer(pos)]
  qtl[, chr  := as.integer(sub("chr", "", chr))]
  qtl[, gene_id := phenotype_id]
  qtl <- qtl[rsid != "."]
  mqtl_df <- rbind(mqtl_df, qtl)
}
cat(sprintf("[%s] QTL data: %d variants across %d genes\n",
            Sys.time(), nrow(mqtl_df), length(unique(mqtl_df$gene_id))))

###############################################################################
# 4. Load variant annotation for allele harmonisation
###############################################################################
variant_ann <- fread(snakemake@params$variant_ann, data.table = FALSE)
variant_ann$chr <- as.integer(variant_ann$chr)
variant_ann <- variant_ann[!is.na(variant_ann$chr) &
                             nchar(variant_ann$ref) == 1 &
                             nchar(variant_ann$alt) == 1, ]
if (!"variant_id_chrpos" %in% colnames(variant_ann)) {
  variant_ann$variant_id_chrpos <- paste(variant_ann$chr, variant_ann$position,
                                         variant_ann$alt, sep = "_")
}

###############################################################################
# 5. Set up LD output directory
###############################################################################
dir.create(ld_output_dir, showWarnings = FALSE, recursive = TRUE)

MHC <- data.frame(chr = 6, start = 28477797, end = 33448354)

###############################################################################
# 6. Run SuSiE per gene-signal pair
###############################################################################
susie_pairs <- unique(abf_susie[, .(cpg, gwas_signal)])
cat(sprintf("[%s] Running SuSiE for %d gene-signal pairs\n", Sys.time(), nrow(susie_pairs)))

results_list <- lapply(seq_len(nrow(susie_pairs)), function(i) {
  gene    <- susie_pairs$cpg[i]
  signal  <- susie_pairs$gwas_signal[i]
  # gwas_signal is e.g. "1:2319680" — extract chromosome for parquet lookup
  signal_chr_from_locus <- as.integer(sub(":.*", "", signal))
  fname   <- sprintf("%s_%s_mQTL_%s_GWAS_%s.txt", gene, signal, tissue, GWAS_ID)

  cat(sprintf("[%s] [%d/%d] %s — %s\n", Sys.time(), i, nrow(susie_pairs), gene, signal))

  # ── QTL data for this gene ──────────────────────────────────────────────
  QTL <- mqtl_df[gene_id == gene]
  if (nrow(QTL) == 0) {
    cat(sprintf("  SKIP: no QTL variants for %s\n", gene))
    return(NULL)
  }
  QTL[, SNP := rsid]
  QTL[, ID  := variant_id]
  QTL[, A1  := ALT]
  QTL[, A2  := REF]
  QTL[, A1_freq := af]
  QTL[, beta := slope]      # parquet column is 'slope', not 'b'
  QTL[, se   := slope_se]   # parquet column is 'slope_se', not 'b_se'
  QTL[, p    := pval_nominal]
  QTL[, n    := QTL_n]

  # ── GWAS data for this signal ────────────────────────────────────────────
  signal_pos <- as.integer(sub(".*:", "", sub("_.*", "", signal)))

  GWAS_win <- GWAS_associations[
    chr == signal_chr_from_locus &
    position >= (signal_pos - 1e6) &
    position <= (signal_pos + 1e6)
  ]
  if (nrow(GWAS_win) == 0) {
    cat(sprintf("  SKIP: no GWAS variants in window for %s\n", signal))
    return(NULL)
  }
  GWAS_win[, SNP     := rsid]
  GWAS_win[, A1      := ea]
  GWAS_win[, A2      := nea]
  GWAS_win[, A1_freq := eaf]
  # beta, se, p already present in GWAS_associations
  if (type == "cc") {
    n_cases    <- GWAS_n[1]
    n_controls <- GWAS_n[2]
    GWAS_win[, n := n_cases + n_controls]
    GWAS_win[, s := n_cases / (n_cases + n_controls)]
  } else {
    GWAS_win[, n := GWAS_n[1]]
  }

  chr_n <- signal_chr_from_locus  # use signal chromosome, not QTL (avoids multi-chr issues)

  # ── Skip MHC ─────────────────────────────────────────────────────────────
  if (chr_n == MHC$chr &&
      any(QTL$pos > MHC$start) && any(QTL$pos < MHC$end)) {
    cat(sprintf("  SKIP: MHC region (chr6)\n"))
    return(NULL)
  }

  # ── Compute LD matrices with plink ───────────────────────────────────────
  QTL_for_ld  <- data.frame(SNP = QTL$SNP,  A2 = QTL$A2,  pos = QTL$pos,
                             stringsAsFactors = FALSE)
  GWAS_for_ld <- data.frame(SNP = GWAS_win$SNP, A2 = GWAS_win$A2,
                             pos = GWAS_win$position, stringsAsFactors = FALSE)
  QTL_for_ld  <- QTL_for_ld[!duplicated(QTL_for_ld$SNP), ]
  GWAS_for_ld <- GWAS_for_ld[!duplicated(GWAS_for_ld$SNP), ]

  ld_qtl  <- compute_ld_matrix(QTL_for_ld,  chr_n, fname,
                                bfile     = qtl_bfile,
                                ld_dir    = ld_output_dir,
                                plink_bin = plink_bin)

  gwas_bfile <- paste0(gwas_bfile_prefix, chr_n)
  ld_gwas <- compute_ld_matrix(GWAS_for_ld, chr_n,
                                paste0("GWAS_", fname),
                                bfile     = gwas_bfile,
                                ld_dir    = ld_output_dir,
                                plink_bin = plink_bin)

  if (is.null(ld_qtl) || is.null(ld_gwas)) {
    cat(sprintf("  SKIP: LD matrix is NULL (plink failed or no variants)\n"))
    return(NULL)
  }

  # ── Align QTL and GWAS to LD SNP sets ────────────────────────────────────
  qtl_snps  <- rownames(ld_qtl)
  gwas_snps <- rownames(ld_gwas)
  QTL_sub   <- QTL[SNP %in% qtl_snps]
  GWAS_sub  <- GWAS_win[SNP %in% gwas_snps]
  GWAS_sub  <- GWAS_sub[!is.na(A1_freq)]
  GWAS_sub  <- unique(GWAS_sub)

  # ── Build coloc dataset lists ─────────────────────────────────────────────
  D1 <- list(
    N       = QTL_n,
    MAF     = QTL_sub$A1_freq[match(qtl_snps, QTL_sub$SNP)],
    beta    = QTL_sub$beta[match(qtl_snps, QTL_sub$SNP)],
    varbeta = QTL_sub$se[match(qtl_snps, QTL_sub$SNP)]^2,
    type    = "quant",
    snp     = qtl_snps,
    pvalue  = QTL_sub$p[match(qtl_snps, QTL_sub$SNP)],
    position= QTL_sub$pos[match(qtl_snps, QTL_sub$SNP)],
    LD      = ld_qtl
  )
  D1$MAF <- ifelse(D1$MAF <= 0.5, D1$MAF, 1 - D1$MAF)

  g_snps <- rownames(ld_gwas)
  D2 <- list(
    N       = if (type == "cc") GWAS_sub$n[1] else GWAS_n[1],
    MAF     = GWAS_sub$A1_freq[match(g_snps, GWAS_sub$SNP)],
    beta    = GWAS_sub$beta[match(g_snps, GWAS_sub$SNP)],
    varbeta = GWAS_sub$se[match(g_snps, GWAS_sub$SNP)]^2,
    type    = type,
    snp     = g_snps,
    position= GWAS_sub$position[match(g_snps, GWAS_sub$SNP)],
    LD      = ld_gwas
  )
  D2$MAF <- ifelse(D2$MAF <= 0.5, D2$MAF, 1 - D2$MAF)
  if (type == "cc") D2$s <- GWAS_sub$s[1]

  # Remove NAs
  keep1 <- !is.na(D1$MAF) & !is.na(D1$beta) & !is.na(D1$varbeta)
  keep2 <- !is.na(D2$MAF) & !is.na(D2$beta) & !is.na(D2$varbeta)
  for (field in c("MAF", "beta", "varbeta", "snp", "pvalue", "position")) {
    if (!is.null(D1[[field]])) D1[[field]] <- D1[[field]][keep1]
  }
  D1$LD <- D1$LD[keep1, keep1]
  for (field in c("MAF", "beta", "varbeta", "snp", "position")) {
    if (!is.null(D2[[field]])) D2[[field]] <- D2[[field]][keep2]
  }
  D2$LD <- D2$LD[keep2, keep2]

  if (length(D1$snp) < 5 || length(D2$snp) < 5) {
    cat(sprintf("  SKIP: too few variants after filtering (QTL: %d, GWAS: %d)\n",
                length(D1$snp), length(D2$snp)))
    return(NULL)
  }

  # ── Run runsusie + coloc.susie ────────────────────────────────────────────
  eqtl_s <- tryCatch(
    runsusie(D1, L = 5, coverage = 0.95, repeat_until_convergence = FALSE),
    error = function(e) {
      cat(sprintf("  QTL L=5 failed (%s), trying L=1\n", e$message))
      tryCatch(runsusie(D1, L = 1, coverage = 0.95, repeat_until_convergence = FALSE),
               error = function(e2) { cat(sprintf("  QTL L=1 failed: %s\n", e2$message)); NULL })
    }
  )
  if (!is.null(eqtl_s) && inherits(eqtl_s, "susie") && is.null(eqtl_s$sets$cs)) {
    eqtl_s <- tryCatch(
      runsusie(D1, L = 1, coverage = 0.95, repeat_until_convergence = FALSE),
      error = function(e) NULL
    )
  }

  gwas_s <- tryCatch(
    runsusie(D2, L = 5, coverage = 0.95, repeat_until_convergence = FALSE),
    error = function(e) {
      cat(sprintf("  GWAS L=5 failed (%s), trying L=1\n", e$message))
      tryCatch(runsusie(D2, L = 1, coverage = 0.95, repeat_until_convergence = FALSE),
               error = function(e2) { cat(sprintf("  GWAS L=1 failed: %s\n", e2$message)); NULL })
    }
  )
  if (!is.null(gwas_s) && inherits(gwas_s, "susie") && is.null(gwas_s$sets$cs)) {
    gwas_s <- tryCatch(
      runsusie(D2, L = 1, coverage = 0.95, repeat_until_convergence = FALSE),
      error = function(e) NULL
    )
  }

  if (!inherits(eqtl_s, "susie") || !inherits(gwas_s, "susie")) {
    cat(sprintf("  SKIP: SuSiE did not converge for both datasets\n"))
    return(NULL)
  }
  if (is.null(eqtl_s$sets$cs) || is.null(gwas_s$sets$cs)) {
    cat(sprintf("  SKIP: no credible sets found\n"))
    return(NULL)
  }

  res <- tryCatch({
    sr <- coloc.susie(gwas_s, eqtl_s, p12 = 1e-05)
    sr$summary$gene_id    <- gene
    sr$summary$gwas_signal <- signal
    sr$summary$tissue     <- tissue
    sr$summary$GWAS_ID    <- GWAS_ID
    sr$summary$Dataset    <- fname
    sr$summary$H4_H3_ratio <- sr$summary$PP.H4.abf / sr$summary$PP.H3.abf
    sr$summary$Colocalise  <- sr$summary$PP.H4.abf > 0.8 |
      (sr$summary$PP.H4.abf > 0.6 & sr$summary$PP.H4.abf < 0.8 &
         sr$summary$H4_H3_ratio > 2)
    as.data.table(sr$summary)
  }, error = function(e) {
    cat(sprintf("  coloc.susie error: %s\n", e$message))
    NULL
  })

  return(res)
})

###############################################################################
# 7. Combine and write results
###############################################################################
results <- rbindlist(Filter(Negate(is.null), results_list), fill = TRUE)
cat(sprintf("[%s] SuSiE complete: %d results from %d/%d pairs\n",
            Sys.time(), nrow(results), sum(!sapply(results_list, is.null)),
            nrow(susie_pairs)))

fwrite(results, output_file, sep = "\t")
cat(sprintf("[%s] Results written to %s\n", Sys.time(), output_file))

sink()
sink(type = "message")
