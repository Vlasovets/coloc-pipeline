.libPaths(c("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/R", .libPaths()))

project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc_clean/"
data.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/"
setwd(project.path)

library(data.table)
library(VariantAnnotation)
library(gwasvcf)
library(magrittr)
library(GenomicRanges)
library(rtracklayer)
library(coloc)
library(dplyr)
library(httr)
library(jsonlite)
library(parallel)

source("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/scripts/Coloc_helper_functions.R")

variant_ann <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/variant_ann_hg38.txt.gz", data.table=FALSE)
variant_ann$chr <- as.integer(variant_ann$chr)
variant_ann <- variant_ann[!is.na(variant_ann$chr) & nchar(variant_ann$ref)==1 & nchar(variant_ann$alt)==1,]

GWAS_list <- list("KNEE" = "GO2_b38_KNEE_ana_new.vcf.bgz",
                  "TKR" = "GO2_b38_TKR_ana_new.vcf.bgz",
                  "ALLOA" = "GO2_b38_ALLOA_ana_new.vcf.bgz")

gwas.b38.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/GO2_sumstats"
qtl.data.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/"

# Lift over independent signals
indep.locus <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/RESPONSE_GO2-Suppl-Tables_31Oct2024.csv", skip=2)
indep.locus <- indep.locus[`Osteoarthritis phenotype` %in% c("KNEE", "TKR", "ALLOA")]
indep.locus[, chromosome:=as.numeric(sub(":.*", "", Variant))]
indep.locus[, position:=as.numeric(sub("_.*", "", sub(".*:", "", Variant)))]
signals <- indep.locus[, .(Loci, rsid, chromosome, position, phenotype=`Osteoarthritis phenotype`)]
signals.b38 <- lifOverFunction_vcf(signals[, .(Loci, rsid, chrom=chromosome, pos=position, phenotype)], lift_down=FALSE)
signals.b38 <- signals.b38[, .(Loci, rsid, chromosome=as.numeric(sub("chr", "", seqnames)), position=as.numeric(start), phenotype)]
fwrite(signals.b38, "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/GO2_index_signals_b38.csv")

for(GWAS_ID in unique(signals.b38$phenotype)){ 
  top <- signals.b38[phenotype==GWAS_ID]
  top$chr <- top$chromosome
  ncases <- ifelse(GWAS_ID=="KNEE", 172256, ifelse(GWAS_ID=="TKR", 48161, 489952))
  ncontrols <- ifelse(GWAS_ID=="KNEE", 1144244, ifelse(GWAS_ID=="TKR", 958463, 1472094))
  GWAS_n <- c(ncases, ncontrols)
  
  print(paste0("Extracting full summary stats from the windows of interest for ", GWAS_ID))
  # Get vcf file path
  if(grepl("vcf", basename(GWAS_list[[GWAS_ID]]))){ # handle cases where I downloaded the vcf directly from the ieugwas database
    vcfname <- file.path(gwas.b38.path, basename(GWAS_list[[GWAS_ID]]))
  } else {
    vcfname <- file.path(gwas.b38.path, paste0(gsub(".gz",".vcf",basename(GWAS_list[[GWAS_ID]])), ".bgz"))
  }
  
  # Define the window length here
  window_length= 1e6 # 1Mb
  
  if(file.exists(file.path("GO2_GWAS_data", paste0(GWAS_ID,".rda")))) load(paste0("GO2_GWAS_data/", GWAS_ID, ".rda"))
  else GWAS_associations <- GWAS_sumstats_extract(vcfname, top[, .(rsid, chr=chromosome, position)], build="hg38", window_length, output_path="GO2_GWAS_data")
  
  #############################################################################
  ############### Now we create a window around the QTL lead variants and find the QTL windows that overlap the GWAS windows
  ############################################################################# 
  for(tissue in c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad")){  # 
    if(!file.exists(paste0("overlap_df_", tissue, "_", GWAS_ID,".rda"))){
      print(paste0("Loading QTL and creating window around lead variant for ", GWAS_ID))
      if(tissue=="high_grade_cartilage") ab <- "hg"
      else if(tissue=="low_grade_cartilage") ab <- "lg"
      else if(tissue=="synovium") ab <- "synovium"
      else if(tissue=="fat_pad") ab <- "fat_pad"
      PATH_IN_M <- paste0(data.path, tissue, "/tensorqtl/output_tensorqtl/nom/")
      mqtl <- fread(paste0(data.path, tissue, "/tensorqtl/output_tensorqtl/perm/egenes_df_", ab, ".txt")) 
      mqtl[, c("chr", "pos", "REF", "ALT", "rsid"):=data.table::tstrsplit(variant_id, "_")]
      mqtl[, pos:=as.integer(pos)]
      mqtl[, chr:=as.integer(sub("chr", "", chr))]
      mqtl[, gene_id:=phenotype_id]
      
      # Filter for qval < 0.05
      mqtl <- mqtl[mqtl$qval < 0.05,]
      
      mqtl_g <- mqtl
      # Find mQTLs that overlap the GWAS regions  
      # Create QTL windows around the lead variant
      mqtl_g <- mqtl_g %>%
        mutate(mqtl_start_coord = ifelse(pos - window_length <= 0, 1, pos - window_length),
               mqtl_end_coord = pos + window_length)
      
      
      # Find overlap between the QTL and GWAS windows
      top_tmp <- top
      top_tmp$position2 <- top_tmp$position
      top_tmp$chr <- as.numeric(as.character(top_tmp$chr))
      
      setDT(mqtl_g)  
      setDT(top_tmp)  
      setkey(mqtl_g, "chr", "mqtl_start_coord", "mqtl_end_coord")
      setkey(top_tmp, "chr", "position", "position2")
      
      
      mqtl_g_subset<-foverlaps(mqtl_g, top_tmp, type="any", nomatch=NULL, mult="all")
      
      top <- as.data.frame(top_tmp)
      
      
      fwrite(data.frame(CpG=unique(mqtl_g_subset$gene_id)), paste0(tissue, "_qtl_GWAS_uniqueCpG_overlap_",GWAS_ID,".txt"), sep="\t")
      
      # Extract the full summary stats for this set of traits. The full data would be too big
      print(paste0("Extracting full QTL sumstats in overlapping window for ",
                   length(unique(mqtl_g_subset$gene_id)),
                   " traits for GWAS ID ", GWAS_ID))
      
      tmp.dt <- data.table()
      for(c in seq(1,22)){
        print(c)
        qtl <- as.data.table(arrow::read_parquet(paste0(PATH_IN_M, ab, ".cis_qtl_pairs.chr", c, ".parquet")))
        qtl <- qtl[phenotype_id %in% unique(mqtl_g_subset$gene_id)]
        if(nrow(qtl)==0) next
        qtl[, c("chr", "pos", "REF", "ALT", "rsid"):=tstrsplit(variant_id, "_")]
        qtl[, maf:=ifelse(af<0.5, af, 1-af)]
        qtl[, pos:=as.integer(pos)]
        qtl[, chr:=as.integer(sub("chr", "", chr))]
        qtl[, gene_id:=phenotype_id]
        qtl <- qtl[rsid!="."]
        tmp.dt <- rbind(tmp.dt, qtl)
      }
      rm(qtl)
      fwrite(tmp.dt, paste0(tissue, "_qtl_GWAS_overlappingSNPs_",GWAS_ID,".txt"))
      
      
      # Clean up and clear some memory
      rm(mqtl)
      rm(tmp.dt)
      rm(mqtl_g_subset)
      gc()
      
      # prepare GWAS genomic ranges
      print(paste0("Preparing GWAS genomic ranges for ", GWAS_ID))
      GWAS_signals_coloc = prepare_gwas_sig(in_gwas_signals = top, 
                                            in_mb_dist = window_length)
      
      
      # Load QTL full summary stats  
      print(paste0("Preparing QTL genomic ranges for ", GWAS_ID))
      mqtl_df <- fread(paste0(tissue, "_qtl_GWAS_overlappingSNPs_",GWAS_ID,".txt"), data.table=FALSE)
      
      ###################################
      # Allele flip
      #####################################
      # Check that the QTL is on the same strand ad the GWAS. If not, flip the allele and the sign of beta
      mqtl_df$variant_id_chrpos = paste(mqtl_df$chr, mqtl_df$pos, mqtl_df$ALT, sep="_")
      mqtl_df$variant_id_rev = paste(mqtl_df$chr, mqtl_df$pos, mqtl_df$REF, sep="_")
      nrow(mqtl_df)
      sum(mqtl_df$variant_id_chrpos %in% variant_ann$variant_id_chrpos)
      sum(mqtl_df$variant_id_rev %in% variant_ann$variant_id_chrpos)
      
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
      mqtl_df = mqtl_df[!is.na(mqtl_df$varid_for_coloc),]  # 158023
      
      
      # GWAS
      GWAS_associations$variant_id_chrpos = paste(GWAS_associations$chr, GWAS_associations$position, GWAS_associations$ea, sep="_")
      GWAS_associations$variant_id_rev = paste(GWAS_associations$chr, GWAS_associations$position, GWAS_associations$nea, sep="_")
      nrow(GWAS_associations)
      sum(GWAS_associations$variant_id_chrpos %in% variant_ann$variant_id_chrpos)
      sum(GWAS_associations$variant_id_rev %in% variant_ann$variant_id_chrpos)
      
      tmp = ifelse(GWAS_associations$variant_id_chrpos %in% variant_ann$variant_id_chrpos, "orig",
                   ifelse(GWAS_associations$variant_id_rev %in% variant_ann$variant_id_chrpos, "rev", NA))
      
      GWAS_associations$varid_for_coloc = ifelse(tmp == "orig" & !is.na(tmp), GWAS_associations$variant_id_chrpos,
                                                 ifelse(tmp == "rev" & !is.na(tmp), GWAS_associations$variant_id_rev, NA))
      GWAS_associations$beta = ifelse(tmp == "orig" & !is.na(tmp), GWAS_associations$beta,
                                      ifelse(tmp == "rev" & !is.na(tmp), (-1) * (GWAS_associations$beta), NA))
      GWAS_associations$ea_coloc = ifelse(tmp == "orig" & !is.na(tmp), GWAS_associations$ea,
                                          ifelse(tmp == "rev" & !is.na(tmp), GWAS_associations$nea, NA))
      GWAS_associations$nea_coloc = ifelse(tmp == "orig" & !is.na(tmp), GWAS_associations$nea,
                                           ifelse(tmp == "rev" & !is.na(tmp), GWAS_associations$ea, NA))
      
      
      GWAS_associations$ea <- GWAS_associations$ea_coloc
      GWAS_associations$nea <- GWAS_associations$nea_coloc
      GWAS_associations = GWAS_associations[!is.na(GWAS_associations$varid_for_coloc),]
      
      save(GWAS_associations, file = file.path("GO2_GWAS_data", paste0(GWAS_ID,"_", tissue, "_GWAS_associations.rda")))
      
      ########################
      # generate genomic ranges for QTL
      cpg_df <- mqtl_g[mqtl_g$gene_id %in% unique(mqtl_df$gene_id), c("gene_id", "chr", "pos")]
      cpg_df$id = cpg_df$gene_id
      cpg_df$start = cpg_df$pos
      cpg_df$end = cpg_df$start
      cpg_df$chr <- gsub("chr", "", cpg_df$chr)
      cpg_df$pos <- cpg_df$start
      
      cpg_df_gr = makeGRangesFromDataFrame(cpg_df,
                                           keep.extra.columns=FALSE,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field=c("chr"),
                                           start.field="start",
                                           end.field=c("end"),
                                           strand.field="strand",
                                           starts.in.df.are.0based=FALSE)                                           
      
      
      
      save(cpg_df_gr, file = paste0("mqtl_df_", tissue, "_", GWAS_ID,".rda"))
      save(mqtl_df, file = paste0("cpg_df_gr_", tissue, "_", GWAS_ID,".rda"))
      
      
      print(paste0("Making QTL-GWAS overlap file for ", GWAS_ID))
      overlap_df = get_overlap_df_small(in_signals_gr = GWAS_signals_coloc$signals_gr, 
                                        in_cpg_df_gr = cpg_df_gr, 
                                        in_cpg_df = cpg_df, 
                                        in_peer_dat = mqtl_df, 
                                        in_gwas_sig = GWAS_signals_coloc$gwas_sigs)
      
      dir.create(paste0("PWCOCO/PWCOCO_", GWAS_ID))
      dir.create(paste0("PWCOCO/PWCOCO_mQTL_", tissue, "_", GWAS_ID))
      dir.create(file.path("PWCOCO/", paste0(GWAS_ID, "_", tissue, "_coloc_rda_files")))
      dir.create(file.path(paste0(GWAS_ID, "_coloc_rda_files")))
      
      save(overlap_df, file= paste0("overlap_df_", tissue, "_", GWAS_ID,".rda"))
    }
    
    else load(file=paste0("overlap_df_", tissue, "_", GWAS_ID,".rda"))
    # for a case control dataset,the proportion of samples in dataset XXX that are cases
    type = ifelse(length(GWAS_n) == 1, "quant", "cc")
    QTL_n=ifelse(tissue=="high_grade_cartilage", 216, ifelse(tissue=="low_grade_cartilage", 263, ifelse(tissue=="synovium", 278, 94)))
    print(paste0("Performing coloc for ", GWAS_ID, " and ", tissue))
    colocFAST_df_results = perform_coloc(overlap_df = overlap_df$overlap_df, in_m_qtl = overlap_df$m_qtl, 
                                         out_path = project.path, tissue = tissue, gwas_n = GWAS_n, in_gwas = GWAS_associations,
                                         gwas_cc_ratio = NULL, GTEX_APP = FALSE, GWAS_type=type, GWAS_ID=GWAS_ID, QTL_n=QTL_n, is.eQTL=TRUE, cores=5)
    
    print(paste0("Writing coloc results for ", GWAS_ID, " in ", tissue))
    fwrite(colocFAST_df_results, file= file.path(paste0(GWAS_ID, "_coloc_rda_files"), paste0("colocABF_df_", tissue, "_results.txt")), sep="\t")
    
    
    # Clean up
    rm(overlap_df)
    rm(mqtl_df)
    rm(cpg_df_gr)
    gc()
  }
}