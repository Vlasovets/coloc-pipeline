#!/usr/bin/env Rscript
#' Stage 7: Combine ABF and SuSiE colocalization results
#'
#' Logic:
#'   - ABF PP4 >= 0.8               → colocalizing (use ABF result)
#'   - ABF PP4 0.25–0.8 (ambiguous) → resolved by SuSiE where available;
#'                                     fallback to ABF if SuSiE had no credible sets
#'   - ABF PP4 < 0.25               → not colocalizing
#'
#' Outputs:
#'   {output_dir}/results/{trait}_combined_coloc_results.txt  — all tests
#'   {output_dir}/results/{trait}_combined_significant.txt    — Colocalise == TRUE
#'
#' Called by Snakemake or directly:
#'   Rscript 7_combine_results.R KNEE \
#'     /lustre/scratch/.../coloc-pipeline/results

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
if (exists("snakemake")) {
  trait      <- snakemake@wildcards[["trait"]]
  output_dir <- snakemake@params[["output_dir"]]
  out_all    <- snakemake@output[["all_results"]]
  out_sig    <- snakemake@output[["sig_results"]]
  log_file   <- snakemake@log[[1]]
  sink(log_file, append = FALSE, split = TRUE)
} else {
  args       <- commandArgs(trailingOnly = TRUE)
  trait      <- if (length(args) >= 1) args[1] else "KNEE"
  output_dir <- if (length(args) >= 2) args[2] else
    "/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"
  out_all    <- file.path(output_dir, "results",
                          sprintf("%s_combined_coloc_results.txt", trait))
  out_sig    <- file.path(output_dir, "results",
                          sprintf("%s_combined_significant.txt", trait))
}

dir.create(dirname(out_all), showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[%s] Combining ABF + SuSiE results for trait: %s\n", Sys.time(), trait))

TISSUES <- c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad")

# ---------------------------------------------------------------------------
# 1. Load ABF results
# ---------------------------------------------------------------------------
abf_files <- file.path(output_dir, "coloc_abf",
                       sprintf("%s.%s.colocABF_results.txt", trait, TISSUES))
names(abf_files) <- TISSUES

abf_list <- lapply(TISSUES, function(tis) {
  f <- abf_files[tis]
  if (!file.exists(f) || file.size(f) == 0) {
    cat(sprintf("  [ABF] %s: not found / empty — skipping\n", tis)); return(NULL)
  }
  dt <- fread(f)
  dt[, tissue  := tis]
  dt[, GWAS_ID := trait]
  cat(sprintf("  [ABF] %s: %d rows\n", tis, nrow(dt)))
  dt
})
abf_all <- rbindlist(abf_list, fill = TRUE)
cat(sprintf("[%s] ABF total: %d rows across %d tissue(s)\n",
            Sys.time(), nrow(abf_all), sum(!sapply(abf_list, is.null))))

if (nrow(abf_all) == 0) stop("No ABF results found — cannot proceed.")

# Standardise PP column name
if (!"PP4" %in% names(abf_all) && "PP.H4.abf" %in% names(abf_all))
  setnames(abf_all, "PP.H4.abf", "PP4")

# ---------------------------------------------------------------------------
# 2. Load SuSiE results
# ---------------------------------------------------------------------------
susie_files <- file.path(output_dir, "coloc_susie",
                         sprintf("%s.%s.colocSuSiE_results.txt", trait, TISSUES))
names(susie_files) <- TISSUES

susie_list <- lapply(TISSUES, function(tis) {
  f <- susie_files[tis]
  if (!file.exists(f) || file.size(f) == 0) {
    cat(sprintf("  [SuSiE] %s: not found / empty — skipping\n", tis)); return(NULL)
  }
  dt <- fread(f)
  cat(sprintf("  [SuSiE] %s: %d rows\n", tis, nrow(dt)))
  dt
})
susie_all <- rbindlist(susie_list, fill = TRUE)
cat(sprintf("[%s] SuSiE total: %d rows\n", Sys.time(), nrow(susie_all)))

# ---------------------------------------------------------------------------
# 3. Classify ABF results
# ---------------------------------------------------------------------------
abf_all[, abf_class := fcase(
  PP4 >= 0.8,              "strong",
  PP4 >= 0.25 & PP4 < 0.8, "ambiguous",
  default = "not_coloc"
)]

cat(sprintf("[%s] ABF classes: strong=%d  ambiguous=%d  not_coloc=%d\n",
            Sys.time(),
            sum(abf_all$abf_class == "strong"),
            sum(abf_all$abf_class == "ambiguous"),
            sum(abf_all$abf_class == "not_coloc")))

# ---------------------------------------------------------------------------
# 4. Resolve ambiguous pairs with SuSiE
# ---------------------------------------------------------------------------
# Build lookup key: gene × gwas_signal × tissue
abf_all[, pair_key := paste(cpg, gwas_signal, tissue, sep = "|")]

if (nrow(susie_all) > 0) {
  # Best SuSiE result per gene × signal × tissue = highest PP.H4.abf
  susie_best <- susie_all[, .(
    susie_PP4        = max(PP.H4.abf, na.rm = TRUE),
    susie_Colocalise = any(Colocalise, na.rm = TRUE),
    susie_nresults   = .N
  ), by = .(gene_id, gwas_signal, tissue)]
  susie_best[, pair_key := paste(gene_id, gwas_signal, tissue, sep = "|")]

  abf_all <- merge(abf_all, susie_best[, .(pair_key, susie_PP4,
                                            susie_Colocalise, susie_nresults)],
                   by = "pair_key", all.x = TRUE)
  cat(sprintf("[%s] Ambiguous pairs resolved by SuSiE: %d / %d\n",
              Sys.time(),
              sum(!is.na(abf_all$susie_PP4) & abf_all$abf_class == "ambiguous"),
              sum(abf_all$abf_class == "ambiguous")))
} else {
  abf_all[, c("susie_PP4", "susie_Colocalise", "susie_nresults") := .(NA_real_, NA, NA_integer_)]
}

# ---------------------------------------------------------------------------
# 5. Assign final PP4 and Colocalise flag
# ---------------------------------------------------------------------------
abf_all[, final_PP4 := fcase(
  abf_class == "strong",                      PP4,
  abf_class == "ambiguous" & !is.na(susie_PP4), susie_PP4,
  abf_class == "ambiguous" & is.na(susie_PP4),  PP4,   # SuSiE not run / no CS
  default = PP4
)]

abf_all[, final_method := fcase(
  abf_class == "strong",                       "ABF",
  abf_class == "ambiguous" & !is.na(susie_PP4), "SuSiE",
  abf_class == "ambiguous" & is.na(susie_PP4),  "ABF_fallback",
  default = "ABF"
)]

abf_all[, Colocalise := fcase(
  final_method == "ABF",         final_PP4 >= 0.8,
  final_method == "SuSiE",       susie_Colocalise == TRUE,
  final_method == "ABF_fallback", final_PP4 >= 0.8,
  default = FALSE
)]

# ---------------------------------------------------------------------------
# 6. Summary
# ---------------------------------------------------------------------------
cat(sprintf("\n[%s] Final summary:\n", Sys.time()))
cat(sprintf("  Total tests:             %d\n", nrow(abf_all)))
cat(sprintf("  Colocalise = TRUE:       %d\n", sum(abf_all$Colocalise, na.rm = TRUE)))
cat(sprintf("    via ABF (strong):      %d\n",
            sum(abf_all$Colocalise & abf_all$final_method == "ABF", na.rm = TRUE)))
cat(sprintf("    via SuSiE:             %d\n",
            sum(abf_all$Colocalise & abf_all$final_method == "SuSiE", na.rm = TRUE)))
cat(sprintf("    via ABF fallback:      %d\n",
            sum(abf_all$Colocalise & abf_all$final_method == "ABF_fallback", na.rm = TRUE)))
cat(sprintf("  Unique signals:          %d\n",
            length(unique(abf_all$gwas_signal[abf_all$Colocalise]))))
cat(sprintf("  Unique genes:            %d\n",
            length(unique(abf_all$cpg[abf_all$Colocalise]))))

by_tissue <- abf_all[Colocalise == TRUE, .N, by = .(GWAS_ID, tissue)]
setorder(by_tissue, GWAS_ID, tissue)
cat("\n  Significant by tissue:\n")
print(by_tissue)

# ---------------------------------------------------------------------------
# 7. Write outputs
# ---------------------------------------------------------------------------
cols_out <- c("cpg", "gwas_signal", "tissue", "GWAS_ID",
              "PP4", "abf_class",
              "susie_PP4", "susie_nresults",
              "final_PP4", "final_method", "Colocalise",
              # retain any other ABF columns that exist
              setdiff(names(abf_all),
                      c("cpg", "gwas_signal", "tissue", "GWAS_ID",
                        "PP4", "abf_class", "susie_PP4", "susie_nresults",
                        "final_PP4", "final_method", "Colocalise", "pair_key",
                        "susie_Colocalise")))
cols_out <- intersect(cols_out, names(abf_all))

fwrite(abf_all[, ..cols_out], out_all, sep = "\t")
fwrite(abf_all[Colocalise == TRUE, ..cols_out], out_sig, sep = "\t")

cat(sprintf("\n[%s] Written:\n  %s\n  %s\n", Sys.time(), out_all, out_sig))
cat(sprintf("[%s] Done.\n", Sys.time()))
