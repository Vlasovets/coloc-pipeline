.libPaths(c("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/R", .libPaths()))

library(data.table)
library(magrittr)
library("GenomicRanges")
library(rtracklayer)
library(coloc)
library(dplyr)
library(Hmisc)

project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/"
data.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/"
setwd(project.path)

source("scripts/Coloc_helper_functions.R")

GWAS_list_nameDescription <- list("KNEE" = "KNEE",
                                  "TKR" = "TKR",
                                  "ALLOA" = "ALLOA")

coloc_res <- c()
for(tissue in c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad")){
  for(GWAS_ID in names(GWAS_list_nameDescription)){
    tmp.dt <- fread(file.path(paste0(GWAS_ID, "_coloc_rda_files"), paste0("colocABF_df_", tissue, "_results.txt")),data.table=FALSE)
    tmp.dt$tissue=tissue
    tmp.dt$GWAS_ID=GWAS_ID
    coloc_res <- rbind(coloc_res, tmp.dt)
  }
}

coloc_res_sign <- coloc_res[coloc_res$PP4 >= 0.8 | (coloc_res$PP4 > 0.6 & coloc_res$PP4 < 0.8 & coloc_res$PP4/coloc_res$PP3 > 2),]

# exclude the results where we already have a colocalizing signal, then extract the other signals to colocalize with susie
coloc_res_susie <- coloc_res %>%
  dplyr::filter(!paste0(coloc_res$cpg, coloc_res$gwas_lead_snp_rs, coloc_res$GWAS_ID, coloc_res$tissue) %in% paste0(coloc_res_sign$cpg, coloc_res_sign$gwas_lead_snp_rs, coloc_res_sign$GWAS_ID, coloc_res_sign$tissue),
                PP4 > 0.25
  )

dim(coloc_res_susie)

coloc_res_susie.summary <- coloc_res_susie %>%
  dplyr::select(cpg, GWAS_ID, tissue) %>%
  unique() %>%
  dplyr::count(GWAS_ID, tissue) %>%
  arrange(desc(n))

dir.create(file.path("Coloc_susie", "LD"))
for(tissue in c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad")){ 
  if(tissue=="high_grade_cartilage") ab <- "cartilage_hg"
  else if(tissue=="low_grade_cartilage") ab <- "cartilage_lg"
  else if(tissue=="synovium") ab <- "synovium"
  else if(tissue=="fat_pad") ab <- "fat_pad"
  
  for(x in seq(1:length(GWAS_list_nameDescription))) {
    
    GWAS_ID <- names(GWAS_list_nameDescription)[x]
    
    coloc_res_susie_tmp <- coloc_res_susie[coloc_res_susie$GWAS_ID == GWAS_ID & coloc_res_susie$tissue==tissue,]
    
    QTL_dir <- paste0("PWCOCO_mQTL_", tissue, "_", GWAS_ID)
    GWAS_dir <- paste0("PWCOCO_", GWAS_ID)
    
    fnames.qtl <- list.files(file.path("PWCOCO", QTL_dir))
    fnames.qtl <- fnames.qtl[grepl(paste(coloc_res_susie_tmp$cpg, collapse="|"),fnames.qtl)]
    
    fnames.gwas <- list.files(file.path("PWCOCO", GWAS_dir))
    fnames.gwas <- fnames.gwas[grepl(tissue ,fnames.gwas)]
    fnames.gwas <- fnames.gwas[grepl(paste(coloc_res_susie_tmp$cpg, collapse="|"),fnames.gwas)]

    susie.res.full <- c()
    count=0
    
    tmp.res <- lapply(1:length(fnames.gwas), FUN = function(i){
      
      print(paste0("Processing CpG ",i, " of ", length(fnames.gwas), " for GWAS ", GWAS_list_nameDescription[[GWAS_ID]]))
      
      QTL <- fread(file.path("PWCOCO", QTL_dir, fnames.qtl[i]), data.table=TRUE)
      GWAS <- fread(file.path("PWCOCO", GWAS_dir, fnames.gwas[i]), data.table=TRUE)

      print("Annotating QTL with chrom and pos")
      QTL[, ID:=SNP]
      QTL[, c("chr", "pos", "SNP"):=tstrsplit(ID, "_")[c(1,2,5)]]
      QTL[, chr:=as.integer(sub("chr", "", chr))]
      QTL[, pos:=as.integer(pos)]
      QTL[, cptid:=paste(paste(chr, pos, sep=":"), A1, A2, sep="_")]
      QTL[, cptid_rev:=paste(paste(chr, pos, sep=":"), A2, A1, sep="_")]
      
      print("Annotating GWAS with rsID from QTL")
      GWAS[, cptid:=SNP]
      GWAS[, SNP:=NULL]
      GWAS <- rbind(merge(GWAS, QTL[, .(SNP, cptid)], by="cptid"), merge(GWAS, QTL[, .(SNP, cptid=cptid_rev)], by="cptid"))
      setkey(QTL, SNP)
      setkey(GWAS, SNP)
      
      #######################
      positions = QTL$pos
      chr = unique(QTL$chr)
      susie.res <- data.frame()
      
      MHCregion=data.frame(chr="6", start=28477797, end=33448354)
      # Only run if not MHC region
      if(!(chr == MHCregion$chr & (sum(positions > MHCregion$start) > 0 & sum(positions < MHCregion$end) > 0))){
        
        print("Loading LD matrix")
        # bfile is the path and prefix of genotype data in plink format    
        if(tissue=="low_grade_cartilage") ld.file <- paste0("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/", tissue, "/CleanGenotypes/P21083A_P21083B.indel_snp.recalibrated.autosomal.hardfiltered.samqc.", ab, "_samples.MAF005CR999.biall_rename")
        else ld.file <- paste0("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/", tissue, "/CleanGenotypes/P21083A_P21083B.indel_snp.recalibrated.autosomal.hardfiltered.samqc.", ab, "_samples.MAF005CR99.biall_rename")
        ld_subset.QTL <- ld_matrix_local_MT(QTL, chr, fnames.gwas[i], bfile=ld.file)
        ld_subset.GWAS <- ld_matrix_local_GO(GWAS, chr, fnames.gwas[i], bfile=paste0("/lustre/groups/itg/shared/referenceData/ukbiobank/chip/bgen2plink/CM/chr", chr))

        # Rename row and column ld_subset.QTL with rsid
        colnames(ld_subset.QTL) <- QTL$SNP[match(colnames(ld_subset.QTL), QTL$ID)]
        rownames(ld_subset.QTL) <- QTL$SNP[match(rownames(ld_subset.QTL), QTL$ID)]
        
        # QTL
        D1 = list(N = QTL$n[1], 
                  MAF = QTL$A1_freq[QTL$SNP %in% rownames(ld_subset.QTL)], 
                  beta = QTL$beta[QTL$SNP %in% rownames(ld_subset.QTL)], 
                  varbeta = (QTL$se[QTL$SNP %in% rownames(ld_subset.QTL)])^2, 
                  type = "quant", 
                  snp = QTL$SNP[QTL$SNP %in% rownames(ld_subset.QTL)], 
                  pvalue =  QTL$p[QTL$SNP %in% rownames(ld_subset.QTL)],
                  stringsAsFactors = F,
                  position= QTL$pos[QTL$SNP %in% rownames(ld_subset.QTL)],
                  LD = ld_subset.QTL
        )
        D1$MAF = ifelse(D1$MAF <= 0.5, D1$MAF, 1 - D1$MAF)
        
        # GWAS
        GWAS <- GWAS[!is.na(GWAS$A1_freq),]
        GWAS <- unique(GWAS)
        D2 = list(N = GWAS$n[1], 
                  MAF = GWAS$A1_freq[GWAS$SNP %in% rownames(ld_subset.GWAS)], 
                  beta = GWAS$beta[GWAS$SNP %in% rownames(ld_subset.GWAS)], 
                  varbeta = (GWAS$se[GWAS$SNP %in% rownames(ld_subset.GWAS)])^2, 
                  type = "quant", 
                  snp = GWAS$SNP[GWAS$SNP %in% rownames(ld_subset.GWAS)], 
                  stringsAsFactors = F,
                  position= GWAS$p[GWAS$SNP %in% rownames(ld_subset.GWAS)],
                  LD = ld_subset.GWAS
        )
        D2$MAF = ifelse(D2$MAF <= 0.5, D2$MAF, 1 - D2$MAF)
        # Change type and add ratio of cases if cc
        if("s" %in% colnames(GWAS)){
          D2$type = "cc"
          D2$s <- GWAS$s[1]
          
        }
        
        susie.res <- data.frame()
        # First, try to finemap with up to 5 credible sets.
        # If this fails, finemap with 1 credible set, which is equivalent to what coloc does
        print(paste0("Running coloc.susie"))
        eqtl.s = tryCatch(runsusie(D1, L=5, coverage=0.95, repeat_until_convergence = FALSE), error = function(e){print("QTL did not converge")})
        gwas.s = tryCatch(runsusie(D2, L=5, coverage=0.95, repeat_until_convergence = FALSE), error = function(e){print("GWAS did not converge")})
        
        if(!class(eqtl.s) == "susie"){
          eqtl.s = tryCatch(runsusie(D1, L=1, coverage=0.95, repeat_until_convergence = FALSE), error = function(e){print("QTL did not converge")})
        }
        else if(is.null(eqtl.s$sets$cs)){
          eqtl.s = tryCatch(runsusie(D1, L=1, coverage=0.95, repeat_until_convergence = FALSE), error = function(e){print("QTL did not converge")})
        }
        
        if(!class(gwas.s) == "susie"){
          gwas.s = tryCatch(runsusie(D2, L=1, coverage=0.95, repeat_until_convergence = FALSE), error = function(e){print("GWAS did not converge")})
        }
        else if(is.null(gwas.s$sets$cs)){
          gwas.s = tryCatch(runsusie(D2, L=1, coverage=0.95, repeat_until_convergence = FALSE), error = function(e){print("GWAS did not converge")})
        }
        
        if(class(eqtl.s) == "susie" & class(gwas.s) == "susie"){
          if(!is.null(eqtl.s$sets$cs) & !is.null(gwas.s$sets$cs)){
            print("Running coloc.susie")
            susie.res=coloc.susie(gwas.s, eqtl.s, p12=1e-05)
            susie.res$summary$Dataset <- fnames.gwas[i]
            susie.res$summary$H4_H3_ratio <- susie.res$summary[,8]/susie.res$summary[,7]
            susie.res$summary$Colocalise <- ifelse(susie.res$summary$PP.H4.abf > 0.8 | (susie.res$summary$PP.H4.abf > 0.6 & susie.res$summary$PP.H4.abf < 0.8 & susie.res$summary$H4_H3_ratio > 2), TRUE, FALSE)
            print(susie.res$summary)
            susie.res <- as.data.frame(susie.res$summary)
          } else{
            susie.res=coloc.susie(gwas.s, eqtl.s, p12=1e-05)
            susie.res$Colocalise <- NA
          }
        } 
        
        return(data.table::as.data.table(susie.res))
        count=count+1
      }
    }
    )
    susie.res.full <- rbindlist(tmp.res, fill=TRUE)
    fwrite(susie.res.full, file= file.path("Coloc_susie", paste0("ColocSusie_results_",GWAS_ID, "_", tissue, ".txt")), sep="\t")
  }
}

