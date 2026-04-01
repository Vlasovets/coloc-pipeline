#!/usr/bin/env Rscript

#' Coloc ABF Analysis Module
#'
#' Performs colocalization analysis using coloc.abf method
#' Adapted from src/3_run_coloc_abf.R for Snakemake workflow

suppressPackageStartupMessages({
  library(data.table)
  library(VariantAnnotation)
  library(gwasvcf)
  library(magrittr)
  library(GenomicRanges)
  library(rtracklayer)
  library(coloc)
  library(dplyr)
  library(parallel)
})

# Source helper functions
source(file.path(snakemake@scriptdir, "Coloc_helper_functions.R"))

# Get parameters from Snakemake
GWAS_ID <- snakemake@params$gwas_id
tissue <- snakemake@params$tissue
QTL_n <- snakemake@params$qtl_n
GWAS_n <- snakemake@params$gwas_n
overlap_file <- snakemake@input$overlaps
gwas_vcf <- snakemake@input$gwas_vcf
output_file <- snakemake@output[[1]]
cores <- snakemake@threads

# Load variant annotation for strand flipping
variant_ann <- fread(snakemake@params$variant_ann, data.table=FALSE)
variant_ann$chr <- as.integer(variant_ann$chr)
variant_ann <- variant_ann[!is.na(variant_ann$chr) & nchar(variant_ann$ref)==1 & nchar(variant_ann$alt)==1,]
if(!"variant_id_chrpos" %in% colnames(variant_ann)) {
  variant_ann$variant_id_chrpos <- paste(variant_ann$chr, variant_ann$position, variant_ann$alt, sep="_")
}

# Define paths — use qtl_nominal_dir from config (correct per-tissue path)
qtl.nominal.dir <- snakemake@params$qtl_nominal_dir
window_length <- 1e6  # 1Mb
type <- ifelse(length(GWAS_n) == 1, "quant", "cc")

cat(sprintf("[%s] Starting coloc ABF analysis for %s - %s\n", 
            Sys.time(), GWAS_ID, tissue))

#############################################################################
# 1. Load overlap data
#############################################################################
cat(sprintf("[%s] Loading overlap data from %s\n", Sys.time(), overlap_file))
load(overlap_file)  # This loads 'overlaps', 'gwas_windows', 'qtl_windows'

if(!exists("overlaps") || nrow(overlaps) == 0) {
  cat(sprintf("[%s] No overlaps found, creating empty result file\n", Sys.time()))
  empty_result <- data.frame()
  fwrite(empty_result, output_file, sep="\t")
  quit(save="no", status=0)
}

# Create overlap_df structure expected by perform_coloc()
# overlaps columns: chr, mqtl_start_coord, mqtl_end_coord, pos, gene_id, qval, position, rsid, Loci
overlap_df_coloc <- data.frame(
  cpg.id = overlaps$gene_id,
  cpg.pos = overlaps$pos,
  signal_gwas.rsid = overlaps$rsid,
  signal_gwas.chr = overlaps$chr,
  signal_gwas.position = overlaps$position,
  stringsAsFactors = FALSE
)
cat(sprintf("[%s] Overlap data: %d gene-signal pairs\n", Sys.time(), nrow(overlap_df_coloc)))

#############################################################################
# 2. Load GWAS data (already prepared from Stage 2)
#############################################################################
cat(sprintf("[%s] Loading GWAS data\n", Sys.time()))
gwas_associations_file <- snakemake@input$gwas_data
load(gwas_associations_file)  # This loads 'GWAS_associations'

#############################################################################
# 3. Load QTL data for this tissue
#############################################################################
cat(sprintf("[%s] Loading QTL data for %s\n", Sys.time(), tissue))

# Read all QTL chromosomes and subset to genes with overlaps
# File pattern in qtl_nominal_dir: cis_qtl.cis_qtl_pairs.chr{N}.parquet
qtl_subset_file <- snakemake@input$qtl_subset

# If we have pre-saved QTL data from Stage 2, use it
if(file.exists(qtl_subset_file)) {
  # The Stage 2 output lists genes to extract (has header 'gene_id')
  qtl_genes <- fread(qtl_subset_file, header=TRUE)$gene_id

  cat(sprintf("[%s] Extracting QTL data for %d genes from %s\n", Sys.time(), length(qtl_genes), qtl.nominal.dir))

  tmp.dt <- data.table()
  for(c in seq(1,22)){
    qtl_file <- file.path(qtl.nominal.dir, paste0("cis_qtl.cis_qtl_pairs.chr", c, ".parquet"))
    if(!file.exists(qtl_file)) next
    
    qtl <- as.data.table(arrow::read_parquet(qtl_file))
    qtl <- qtl[phenotype_id %in% qtl_genes]
    if(nrow(qtl)==0) next
    
    qtl[, c("chr", "pos", "REF", "ALT", "rsid"):=tstrsplit(variant_id, "_")]
    qtl[, maf:=ifelse(af<0.5, af, 1-af)]
    qtl[, pos:=as.integer(pos)]
    qtl[, chr:=as.integer(sub("chr", "", chr))]
    qtl[, gene_id:=phenotype_id]
    qtl <- qtl[rsid!="."]
    tmp.dt <- rbind(tmp.dt, qtl)
  }
  mqtl_df <- tmp.dt
  rm(tmp.dt)
} else {
  stop("QTL subset file not found. Stage 2 must be run first.")
}

#############################################################################
# 4. Allele flipping for QTL data
#############################################################################
cat(sprintf("[%s] Performing allele flipping for QTL data\n", Sys.time()))

mqtl_df$variant_id_chrpos = paste(mqtl_df$chr, mqtl_df$pos, mqtl_df$ALT, sep="_")
mqtl_df$variant_id_rev = paste(mqtl_df$chr, mqtl_df$pos, mqtl_df$REF, sep="_")

tmp = ifelse(mqtl_df$variant_id_chrpos %in% variant_ann$variant_id_chrpos, "orig",
             ifelse(mqtl_df$variant_id_rev %in% variant_ann$variant_id_chrpos, "rev", NA))

mqtl_df$varid_for_coloc = ifelse(tmp == "orig" & !is.na(tmp), mqtl_df$variant_id_chrpos,
                                 ifelse(tmp == "rev" & !is.na(tmp), mqtl_df$variant_id_rev, NA))
mqtl_df$slope = ifelse(tmp == "orig" & !is.na(tmp), mqtl_df$slope,
                       ifelse(tmp == "rev" & !is.na(tmp), (-1) * (mqtl_df$slope), NA))
mqtl_df$ALT_coloc = ifelse(tmp == "orig" & !is.na(tmp), mqtl_df$ALT,
                           ifelse(tmp == "rev" & !is.na(tmp), mqtl_df$REF, NA))
mqtl_df$REF_coloc = ifelse(tmp == "orig" & !is.na(tmp), mqtl_df$REF,
                           ifelse(tmp == "rev" & !is.na(tmp), mqtl_df$ALT, NA))

mqtl_df$ALT <- mqtl_df$ALT_coloc
mqtl_df$REF <- mqtl_df$REF_coloc
mqtl_df = mqtl_df[!is.na(mqtl_df$varid_for_coloc),]

# Split QTL data by gene for perform_coloc function (expects a list)
in_m_qtl <- split(mqtl_df, mqtl_df$gene_id)

#############################################################################
# 5. Run colocalization using helper function
#############################################################################
cat(sprintf("[%s] Running colocalization analysis\n", Sys.time()))

# Create output directories
coloc_dir <- dirname(output_file)
out_path <- dirname(coloc_dir)
dir.create(coloc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_path, paste0(GWAS_ID, "_coloc_rda_files")), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_path, "coloc", "Coloc_sumstats", paste0("Coloc_mQTL_", GWAS_ID)), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_path, "coloc", "Coloc_sumstats", paste0("Coloc_", GWAS_ID)), showWarnings = FALSE, recursive = TRUE)

# Run coloc using the helper function
colocFAST_df_results <- perform_coloc(
  overlap_df = overlap_df_coloc,
  in_m_qtl = in_m_qtl,
  out_path = out_path,
  gwas_n = GWAS_n,
  in_gwas = GWAS_associations,
  GWAS_type = type,
  GWAS_ID = GWAS_ID,
  QTL_n = QTL_n,
  is.eQTL = TRUE,
  cores = cores
)

#############################################################################
# 6. Save results
#############################################################################
cat(sprintf("[%s] Saving results to %s\n", Sys.time(), output_file))
fwrite(colocFAST_df_results, output_file, sep="\t")

cat(sprintf("[%s] Coloc ABF analysis complete for %s - %s\n", 
            Sys.time(), GWAS_ID, tissue))
cat(sprintf("[%s] Total tests: %d\n", Sys.time(), nrow(colocFAST_df_results)))
