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

#' Compute QTL LD matrix using variant_id (chr_pos_ref_alt_rsid) as plink ID
#'
#' @param snp_df  data.table with columns SNP (= variant_id), A2 (REF), pos
#' @param chr     chromosome number
#' @param fname   base filename for temp plink outputs
#' @param bfile   plink bfile prefix (no extension); bim uses chr_pos_ref_alt_rsid IDs
#' @param ld_dir  output directory for plink LD files
#' @param plink_bin  path to plink binary
#' @return named LD matrix (SNP x SNP correlation, names = variant_id) or NULL
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

  ld[is.nan(ld) | is.na(ld)] <- 0  # plink writes nan for monomorphic variants
  diag(ld) <- 1                     # ensure self-correlation is always 1

  colnames(ld) <- bim$rsid
  rownames(ld) <- bim$rsid

  keep <- rownames(ld) %in% snp_df$SNP
  ld[keep, keep, drop = FALSE]
}

#' Compute GWAS LD matrix using a two-step approach (avoids rsid requirement)
#'
#' Step 1: Extract region bim only (--make-just-bim, no LD) to get variant IDs+positions.
#' Step 2: Match bim to GWAS_win by position+allele to identify relevant variants (~hundreds).
#' Step 3: Extract only matched variants and compute --r square LD.
#'
#' This avoids computing a full LD matrix for all ~60K+ chip variants in a 2Mb window.
#'
#' @param gwas_dt  data.table with chr, position, ea, nea columns (GWAS variants)
#' @param chr      chromosome number (integer)
#' @param bp_start region start (bp)
#' @param bp_end   region end (bp)
#' @param fname    base filename for temp plink outputs
#' @param bfile    plink bfile prefix (bim uses rsids as variant IDs)
#' @param ld_dir   output directory
#' @param plink_bin path to plink binary
#' @return list(ld = LD matrix, gwas_sub = matched GWAS variants) or NULL
compute_gwas_ld_region <- function(gwas_dt, chr, bp_start, bp_end,
                                   fname, bfile, ld_dir, plink_bin) {
  dir.create(ld_dir, showWarnings = FALSE, recursive = TRUE)
  base      <- file.path(ld_dir, paste0("GWAS_", sub("\\.[^.]*$", "", fname)))
  f_bim_tmp <- paste0(base, ".region")
  f_ext     <- paste0(base, ".SNP.extract")
  f_ld      <- paste0(base, ".LD")

  # Step 1: Get bim for full region (no LD yet — cheap)
  cmd1 <- paste(plink_bin,
    "--bfile", bfile, "--const-fid",
    "--chr", chr, "--from-bp", bp_start, "--to-bp", bp_end,
    "--threads 1", "--make-just-bim",
    "--out", f_bim_tmp, "--silent"
  )
  system(cmd1, ignore.stdout = TRUE, ignore.stderr = TRUE)

  bim_full <- tryCatch(
    as.data.table(read.table(paste0(f_bim_tmp, ".bim"), header = FALSE,
               col.names = c("chr", "snp_id", "cM", "position", "a1", "a2"))),
    error = function(e) NULL
  )
  if (is.null(bim_full) || nrow(bim_full) == 0) return(NULL)

  # Step 2: Match GWAS variants to bim by position + allele
  gwas_dt <- copy(gwas_dt)
  gwas_dt[, match_key_fwd := paste(position, toupper(ea),  toupper(nea), sep = "_")]
  gwas_dt[, match_key_rev := paste(position, toupper(nea), toupper(ea),  sep = "_")]
  bim_full[, match_key := paste(position, toupper(a1), toupper(a2), sep = "_")]

  gwas_dt[, snp_id := bim_full$snp_id[match(match_key_fwd, bim_full$match_key)]]
  gwas_dt[is.na(snp_id),
          snp_id := bim_full$snp_id[match(match_key_rev[is.na(snp_id)], bim_full$match_key)]]
  gwas_dt <- gwas_dt[!is.na(snp_id)]

  if (nrow(gwas_dt) == 0) return(NULL)

  # Step 3: Extract only matched variants and compute LD
  write.table(gwas_dt$snp_id, f_ext,
              col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

  cmd2 <- paste(plink_bin,
    "--bfile", bfile, "--const-fid",
    "--extract", f_ext,
    "--threads 1", "--r square",
    "--out", f_ld, "--make-just-bim", "--silent"
  )
  ret <- system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (ret != 0 || !file.exists(paste0(f_ld, ".ld"))) return(NULL)

  bim <- tryCatch(
    as.data.table(read.table(paste0(f_ld, ".bim"), header = FALSE,
               col.names = c("chr", "snp_id", "cM", "position", "a1", "a2"))),
    error = function(e) NULL
  )
  if (is.null(bim)) return(NULL)

  ld <- tryCatch(
    as.matrix(read.table(paste0(f_ld, ".ld"), header = FALSE)),
    error = function(e) NULL
  )
  if (is.null(ld)) return(NULL)

  ld[is.nan(ld) | is.na(ld)] <- 0  # plink writes nan for monomorphic variants
  diag(ld) <- 1

  colnames(ld) <- bim$snp_id
  rownames(ld) <- bim$snp_id

  gwas_dt <- gwas_dt[snp_id %in% rownames(ld)]
  if (nrow(gwas_dt) == 0) return(NULL)

  keep <- rownames(ld) %in% gwas_dt$snp_id
  list(
    ld      = ld[keep, keep, drop = FALSE],
    gwas_sub = gwas_dt
  )
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
# 4. Save data checkpoint (for fast re-runs / interactive debugging)
###############################################################################
checkpoint_file <- paste0(output_file, ".checkpoint.rds")
saveRDS(list(GWAS_associations = GWAS_associations,
             mqtl_df           = mqtl_df,
             abf_susie         = abf_susie),
        checkpoint_file)
cat(sprintf("[%s] Checkpoint saved to %s\n", Sys.time(), checkpoint_file))

###############################################################################
# 5. Set up LD output directory
###############################################################################
dir.create(ld_output_dir, showWarnings = FALSE, recursive = TRUE)

MHC <- data.frame(chr = 6, start = 28477797, end = 33448354)

###############################################################################
# 6. Run SuSiE per gene-signal pair
###############################################################################
susie_pairs <- unique(abf_susie[, .(cpg, gwas_signal, coloc_region_chr, coloc_region_start, coloc_region_end)])
cat(sprintf("[%s] Running SuSiE for %d gene-signal pairs\n", Sys.time(), nrow(susie_pairs)))

results_list <- lapply(seq_len(nrow(susie_pairs)), function(i) {
  gene    <- susie_pairs$cpg[i]
  signal  <- susie_pairs$gwas_signal[i]
  signal_chr_from_locus <- susie_pairs$coloc_region_chr[i]
  region_start          <- susie_pairs$coloc_region_start[i]
  region_end            <- susie_pairs$coloc_region_end[i]
  fname   <- sprintf("%s_%s_mQTL_%s_GWAS_%s.txt", gene, signal, tissue, GWAS_ID)

  cat(sprintf("[%s] [%d/%d] %s — %s\n", Sys.time(), i, nrow(susie_pairs), gene, signal))

  # ── QTL data for this gene, restricted to the coloc region ────────────
  # Restrict to the same region used for GWAS LD to keep D1/D2 symmetric.
  # Using all ~4000 cis-window variants inflates the region and makes SuSiE
  # unable to form credible sets (prior too diffuse; signal strength diluted).
  QTL <- mqtl_df[gene_id == gene &
                 pos >= region_start &
                 pos <= region_end]
  if (nrow(QTL) == 0) {
    cat(sprintf("  SKIP: no QTL variants in region %s:%d-%d for %s\n",
                signal_chr_from_locus, region_start, region_end, gene))
    return(NULL)
  }
  # SNP = variant_id (full bim ID: chr_pos_ref_alt_rsid) — required for plink --extract
  # The QTL bim uses this format, not plain rsids
  QTL[, SNP := variant_id]
  QTL[, ID  := variant_id]
  QTL[, A1  := ALT]
  QTL[, A2  := REF]
  QTL[, A1_freq := af]
  QTL[, beta := slope]      # parquet column is 'slope', not 'b'
  QTL[, se   := slope_se]   # parquet column is 'slope_se', not 'b_se'
  QTL[, p    := pval_nominal]
  QTL[, n    := QTL_n]

  # ── GWAS data for this signal ────────────────────────────────────────────
  GWAS_win <- GWAS_associations[
    chr == signal_chr_from_locus &
    position >= region_start &
    position <= region_end
  ]
  if (nrow(GWAS_win) == 0) {
    cat(sprintf("  SKIP: no GWAS variants in window for %s\n", signal))
    return(NULL)
  }
  GWAS_win[, A1_freq := eaf]
  if (type == "cc") {
    n_cases    <- GWAS_n[1]
    n_controls <- GWAS_n[2]
    GWAS_win[, n := n_cases + n_controls]
    GWAS_win[, s := n_cases / (n_cases + n_controls)]
  } else {
    GWAS_win[, n := GWAS_n[1]]
  }

  chr_n <- signal_chr_from_locus

  # ── Skip MHC ─────────────────────────────────────────────────────────────
  if (chr_n == MHC$chr &&
      any(QTL$pos > MHC$start) && any(QTL$pos < MHC$end)) {
    cat(sprintf("  SKIP: MHC region (chr6)\n"))
    return(NULL)
  }

  # ── QTL LD: extract by variant_id (full bim ID: chr_pos_ref_alt_rsid) ────
  QTL_for_ld <- data.frame(SNP = QTL$SNP, A2 = QTL$A2, pos = QTL$pos,
                            stringsAsFactors = FALSE)
  QTL_for_ld <- QTL_for_ld[!duplicated(QTL_for_ld$SNP), ]

  ld_qtl <- compute_ld_matrix(QTL_for_ld, chr_n, fname,
                               bfile     = qtl_bfile,
                               ld_dir    = ld_output_dir,
                               plink_bin = plink_bin)
  if (is.null(ld_qtl)) {
    cat(sprintf("  SKIP: QTL LD matrix is NULL (plink failed or no variants)\n"))
    return(NULL)
  }

  # ── GWAS LD: region-based extraction (UKB bfile uses rsids; GWAS data lacks them)
  gwas_bfile <- paste0(gwas_bfile_prefix, chr_n)
  gwas_ld_res <- compute_gwas_ld_region(
    gwas_dt   = GWAS_win,
    chr       = chr_n,
    bp_start  = region_start,
    bp_end    = region_end,
    fname     = fname,
    bfile     = gwas_bfile,
    ld_dir    = ld_output_dir,
    plink_bin = plink_bin
  )
  if (is.null(gwas_ld_res)) {
    cat(sprintf("  SKIP: GWAS LD matrix is NULL (plink failed or no variants)\n"))
    return(NULL)
  }
  ld_gwas  <- gwas_ld_res$ld
  GWAS_sub <- gwas_ld_res$gwas_sub   # already matched to bim by position+allele

  # ── Align QTL to LD SNP set ──────────────────────────────────────────────
  qtl_snps <- rownames(ld_qtl)
  QTL_sub  <- QTL[SNP %in% qtl_snps]
  GWAS_sub <- GWAS_sub[!is.na(A1_freq)]
  GWAS_sub <- GWAS_sub[!duplicated(snp_id)]  # ensure 1:1 snp_id → data row
  # Keep only GWAS_sub rows whose snp_id is in the LD matrix
  GWAS_sub <- GWAS_sub[snp_id %in% rownames(ld_gwas)]
  # Trim LD to only matched GWAS variants
  ld_gwas  <- ld_gwas[GWAS_sub$snp_id, GWAS_sub$snp_id, drop = FALSE]

  cat(sprintf("  QTL: %d vars in LD, GWAS: %d vars in LD\n",
              length(qtl_snps), nrow(GWAS_sub)))
  # Diagnostic: signal strength
  qtl_z  <- abs(QTL_sub$beta / QTL_sub$se)
  gwas_z <- abs(GWAS_sub$beta / GWAS_sub$se)
  cat(sprintf("  QTL top |z|: %.2f  GWAS top |z|: %.2f\n",
              max(qtl_z, na.rm=TRUE), max(gwas_z, na.rm=TRUE)))

  # ── Build coloc dataset lists ─────────────────────────────────────────────
  # Order QTL_sub to match LD rownames exactly
  QTL_sub <- QTL_sub[match(qtl_snps, SNP)]
  D1 <- list(
    N       = QTL_n,
    MAF     = QTL_sub$A1_freq,
    beta    = QTL_sub$beta,
    varbeta = QTL_sub$se^2,
    type    = "quant",
    snp     = QTL_sub$SNP,
    pvalue  = QTL_sub$p,
    position= QTL_sub$pos,
    LD      = ld_qtl
  )
  D1$MAF <- ifelse(D1$MAF <= 0.5, D1$MAF, 1 - D1$MAF)

  # GWAS_sub rows are already in ld_gwas row order after the trim above
  D2 <- list(
    N       = if (type == "cc") GWAS_sub$n[1] else GWAS_n[1],
    MAF     = GWAS_sub$A1_freq,
    beta    = GWAS_sub$beta,
    varbeta = GWAS_sub$se^2,
    type    = type,
    snp     = GWAS_sub$snp_id,
    position= GWAS_sub$position,
    LD      = ld_gwas
  )
  D2$MAF <- ifelse(D2$MAF <= 0.5, D2$MAF, 1 - D2$MAF)
  if (type == "cc") D2$s <- GWAS_sub$s[1]

  # Remove NAs — must subset LD and all data vectors consistently
  keep1 <- !is.na(D1$MAF) & !is.na(D1$beta) & !is.na(D1$varbeta) &
           !is.na(D1$position)
  keep2 <- !is.na(D2$MAF) & !is.na(D2$beta) & !is.na(D2$varbeta) &
           !is.na(D2$position)
  for (field in c("MAF", "beta", "varbeta", "snp", "pvalue", "position")) {
    if (!is.null(D1[[field]])) D1[[field]] <- D1[[field]][keep1]
  }
  D1$LD <- D1$LD[keep1, keep1, drop = FALSE]
  for (field in c("MAF", "beta", "varbeta", "snp", "position")) {
    if (!is.null(D2[[field]])) D2[[field]] <- D2[[field]][keep2]
  }
  D2$LD <- D2$LD[keep2, keep2, drop = FALSE]

  # Verify consistency: runsusie requires length(snp) == nrow(LD)
  stopifnot(length(D1$snp) == nrow(D1$LD))
  stopifnot(length(D2$snp) == nrow(D2$LD))

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
  qtl_pip_max  <- if (inherits(eqtl_s,"susie")) max(eqtl_s$pip,  na.rm=TRUE) else NA
  gwas_pip_max <- if (inherits(gwas_s, "susie")) max(gwas_s$pip, na.rm=TRUE) else NA
  cat(sprintf("  QTL max PIP: %.4f  GWAS max PIP: %.4f  QTL CS: %s  GWAS CS: %s\n",
              qtl_pip_max, gwas_pip_max,
              ifelse(is.null(eqtl_s$sets$cs), "NULL", length(eqtl_s$sets$cs)),
              ifelse(is.null(gwas_s$sets$cs),  "NULL", length(gwas_s$sets$cs))))
  if (is.null(eqtl_s$sets$cs) || is.null(gwas_s$sets$cs)) {
    cat(sprintf("  SKIP: no credible sets found\n"))
    return(NULL)
  }

  # ── Harmonise lbf_variable column names before coloc.susie ───────────────
  # QTL snp names: chr_pos_REF_ALT_rsid  (from variant_id in QTL bim)
  # GWAS snp names: plain rsid           (from UKB chip bim)
  # coloc.susie uses intersect(colnames(bf1), colnames(bf2)) to find shared
  # SNPs. With 0 overlap it returns data.table(nsnps=NA) and the downstream
  # := call crashes. Fix: rename both to "pos_NNNNN" positional keys.
  qtl_pos_key <- paste0("pos_",
    sub("^[^_]+_([0-9]+)_.*", "\\1", colnames(eqtl_s$lbf_variable)))
  colnames(eqtl_s$lbf_variable) <- qtl_pos_key

  # For GWAS, look up position from GWAS_sub (snp_id → position)
  gwas_pos_key <- paste0("pos_",
    GWAS_sub$position[match(colnames(gwas_s$lbf_variable), GWAS_sub$snp_id)])
  colnames(gwas_s$lbf_variable) <- gwas_pos_key

  n_shared <- length(intersect(qtl_pos_key, gwas_pos_key))
  cat(sprintf("  Shared positional SNPs for coloc.susie: %d\n", n_shared))
  if (n_shared == 0) {
    cat(sprintf("  SKIP: no positional overlap between QTL and GWAS lbf matrices\n"))
    return(NULL)
  }

  sr <- try(coloc.susie(gwas_s, eqtl_s, p12 = 1e-05), silent = TRUE)
  if (inherits(sr, "try-error")) {
    cat(sprintf("  coloc.susie error: %s\n",
                conditionMessage(attr(sr, "condition"))))
    res <- NULL
  } else {
    cat(sprintf("  coloc.susie OK — summary class: %s  nrow: %d  cols: %s\n",
                paste(class(sr$summary), collapse="/"),
                ifelse(is.data.frame(sr$summary), nrow(sr$summary), NA),
                paste(names(sr$summary), collapse=",")))
    if (is.data.frame(sr$summary) && nrow(sr$summary) > 0)
      cat(sprintf("  PP.H4 values: %s\n",
                  paste(round(sr$summary$PP.H4.abf, 4), collapse=",")))
    smry <- as.data.frame(sr$summary)
    smry$gene_id     <- gene
    smry$gwas_signal <- signal
    smry$tissue      <- tissue
    smry$GWAS_ID     <- GWAS_ID
    smry$Dataset     <- fname
    smry$H4_H3_ratio <- smry$PP.H4.abf / smry$PP.H3.abf
    smry$Colocalise  <- smry$PP.H4.abf > 0.8 |
      (smry$PP.H4.abf > 0.6 & smry$PP.H4.abf < 0.8 &
         smry$H4_H3_ratio > 2)
    res <- as.data.table(smry)
  }

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
