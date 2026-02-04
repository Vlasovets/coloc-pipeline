project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/"
data.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/"
setwd(project.path)

library(data.table)
library(dplyr)
library(rtracklayer)

GWAS_list <- list("KNEE" = "GO2_b38_KNEE_ana_new.vcf.bgz",
                  "TKR" = "GO2_b38_TKR_ana_new.vcf.bgz",
                  "ALLOA" = "GO2_b38_ALLOA_ana_new.vcf.bgz")

# Need: "Gene", "Chrom", "Start", "End", "Coding"
source(paste0(project.path, "scripts/plot_functions.R"))
GRCh38_Genes <- as.data.table(rtracklayer::import("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/gencode.v46.basic.annotation.gtf.gz"))
GRCh38_Genes <- GRCh38_Genes[type=="gene" & !(seqnames %in% c("chrX", "chrY", "chrM")), .(Chrom=as.integer(sub("chr", "", seqnames)), Start=start, End=end, Gene=gene_name, Coding="Coding")]
GRCh38_Genes <- GRCh38_Genes[!grepl("ENSG", Gene) & End-Start>5000]
GRCh38_Genes <- as.data.frame(GRCh38_Genes)

plot.coloc <- function(gwas.data, qtl.data, traits, res.coloc, genes.data, plot.path="plots/", ld.path="LD/", out.file.suffix="", range=1e+6){
  sel.snp <- res.coloc$coloc_lead_snp_rs
  sel.chr <- res.coloc$coloc_lead_snp_chr
  sel.pos <- res.coloc$coloc_lead_snp_pos
  sel.PP4 <- res.coloc$PP4
  
  # Get data
  trait1.region <- unique(gwas.data[CHR==sel.chr & POS<sel.pos+range & POS>sel.pos-range, .(CHR, SNP, P, BP=POS, logP=-log10(as.numeric(P)))], by="SNP")
  trait2.region <- unique(qtl.data[chr==sel.chr & pos<sel.pos+range & pos>sel.pos-range, .(CHR=chr, SNP=rsid, P=pval_nominal, BP=pos, logP=-log10(as.numeric(pval_nominal)))], by="SNP")
  data.lst <- list(trait1.region, trait2.region)
  names(data.lst) <-  c(traits[1], traits[2])

  # Get LD matrix
  if(!file.exists(paste0(ld.path, sel.snp, "_LD.ld"))){
    print("Calculating LD with plink")
    tryCatch({
      system(paste0("plink --bim /lustre/groups/itg/teams/data-informatics/projects/ukbb/genotypes/imputed/EGAD00010001474/GO2_var_ids_fake_sample_IDs/ukbb.imputed.v3.chr", sel.chr, ".bim.orig --bed /lustre/groups/itg/teams/data-informatics/projects/ukbb/genotypes/imputed/EGAD00010001474/GO2_var_ids_fake_sample_IDs/ukbb.imputed.v3.chr", sel.chr, ".bed --fam /lustre/groups/itg/teams/data-informatics/projects/ukbb/genotypes/imputed/EGAD00010001474/GO2_var_ids_fake_sample_IDs/ukbb.imputed.v3.chr", sel.chr, ".fam --from ",  trait1.region[BP==min(trait1.region$BP), SNP], " --to ", trait1.region[BP==max(trait1.region$BP), SNP], " --threads 10 --r2 inter-chr --ld-snp ", sel.snp, "  --ld-window-r2 0 --out ", ld.path, sel.snp, "_LD"))
    })
  }
  if(!file.exists(paste0(ld.path, sel.snp, "_LD.ld"))) ld.file <- data.table::fread(paste0(ld.path, "rs10407349_LD.ld"))[SNP_A==sel.snp | SNP_B==sel.snp]
  else ld.file <- data.table::fread(paste0(ld.path, sel.snp, "_LD.ld"))[SNP_A==sel.snp | SNP_B==sel.snp]
  ld.file[SNP_B==sel.snp, c("SNP_B", "SNP_A") := .(SNP_A, SNP_B)]

  # Regional association plot
  locus.zoom(data = data.lst,
             offset_bp = range,
             genes.data = genes.data,
             file.name = paste0(plot.path, paste(traits[1], traits[2], sel.snp, sep="_"), out.file.suffix, ".png"),
             snp=sel.snp,
             ignore.lead=TRUE,
             ld.file=ld.file,
             pp="PP4",
             pp.value=round(unname(sel.PP4), digits=3),
             nplots=TRUE)
}

plot.coloc.susie <- function(gwas.data, qtl.data, traits, res.coloc, genes.data, plot.path="plots_susie/", ld.path="LD/", out.file.suffix="", range=1e+6){
  sel.snp <- unique(res.coloc$hit2)
  sel.chr <- unique(res.coloc$chr)
  sel.pos <- unique(res.coloc$pos)
  sel.PP4 <- unique(res.coloc$PP.H4.abf)
  
  # Get data
  trait1.region <- unique(gwas.data[chr==sel.chr & position<sel.pos+range & position>sel.pos-range, .(CHR=chr, SNP=rsid, P=p, BP=position, logP=-log10(as.numeric(p)))], by="SNP")
  trait2.region <- unique(qtl.data[chr==sel.chr & pos<sel.pos+range & pos>sel.pos-range, .(CHR=chr, SNP=rsid, P=pval_nominal, BP=pos, logP=-log10(as.numeric(pval_nominal)))], by="SNP")
  data.lst <- list(trait1.region, trait2.region)
  names(data.lst) <-  c(traits[1], traits[2])
  
  # Get LD matrix
  if(!file.exists(paste0(ld.path, sel.snp, "_LD.ld"))){
    print("Calculating LD with plink")
    tryCatch({
      system(paste0("~/plink --bim /lustre/groups/itg/teams/data-informatics/projects/ukbb/genotypes/imputed/EGAD00010001474/GO2_var_ids_fake_sample_IDs/ukbb.imputed.v3.chr", sel.chr, ".bim.orig --bed /lustre/groups/itg/teams/data-informatics/projects/ukbb/genotypes/imputed/EGAD00010001474/GO2_var_ids_fake_sample_IDs/ukbb.imputed.v3.chr", sel.chr, ".bed --fam /lustre/groups/itg/teams/data-informatics/projects/ukbb/genotypes/imputed/EGAD00010001474/GO2_var_ids_fake_sample_IDs/ukbb.imputed.v3.chr", sel.chr, ".fam --from ",  trait1.region[BP==min(trait1.region$BP), SNP], " --to ", trait1.region[BP==max(trait1.region$BP), SNP], " --threads 10 --r2 inter-chr --ld-snp ", sel.snp, "  --ld-window-r2 0 --out ", ld.path, sel.snp, "_LD"))
    })
  }
  if(!file.exists(paste0(ld.path, sel.snp, "_LD.ld"))) ld.file <- data.table::fread(paste0(ld.path, "rs2076328_LD.ld"))[SNP_A==sel.snp | SNP_B==sel.snp]
  else ld.file <- data.table::fread(paste0(ld.path, sel.snp, "_LD.ld"))[SNP_A==sel.snp | SNP_B==sel.snp]
  ld.file[SNP_B==sel.snp, c("SNP_B", "SNP_A") := .(SNP_A, SNP_B)]
  
  # Regional association plot
  locus.zoom(data = data.lst,
             offset_bp = range,
             genes.data = genes.data,
             file.name = paste0(plot.path, paste(traits[1], traits[2], sel.snp, sep="_"), out.file.suffix, "_susie.png"),
             snp = sel.snp,
             ignore.lead = TRUE,
             ld.file = ld.file,
             pp = "PP4",
             pp.value = round(unname(sel.PP4), digits=3),
             nplots = TRUE)
}

signals.b38 <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/GO2_index_signals_b38.csv")

#########################################################################
############# Check coloc.abf results
#########################################################################
coloc_res <- c()
for(tissue in c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad")){
  for(GWAS_ID in names(GWAS_list)){
    tmp.dt <- fread(file.path(paste0(GWAS_ID, "_coloc_rda_files"), paste0("colocABF_df_", tissue, "_results.txt")))
    tmp.dt$tissue=tissue
    tmp.dt$GWAS_ID=GWAS_ID
    coloc_res <- rbind(coloc_res, tmp.dt)
  }
}

coloc_res <- merge(coloc_res, unique(signals.b38[, .(Loci, rsid)]), by.x="gwas_signal", by.y="rsid", all.x=T)
length(unique(coloc_res$gwas_signal))  # 357
length(unique(coloc_res$Loci))  # 211
length(unique(coloc_res$cpg))  # 2483
fwrite(coloc_res, "all_coloc_abf_results.csv")

coloc_res_sign <- coloc_res[PP4 >= 0.8 | (PP4 > 0.6 & PP4 < 0.8 & PP4/PP3 > 2)]
fwrite(coloc_res_sign, "GO2_significant_coloc_results_abf.txt")
fwrite(coloc_res_sign, "GO2_significant_coloc_results_abf.csv")
length(unique(coloc_res_sign$gwas_signal))  # 147 --> 109
length(unique(coloc_res_sign$Loci))  # 
length(unique(coloc_res_sign$cpg))  # 164 --> 110

#########################################################################
############# Plot coloc.abf() results
#########################################################################
for(g in unique(coloc_res_sign$GWAS_ID)){
  load(paste0("GO2_GWAS_data/", g, ".rda"))  # GWAS_associations
  gwas <- as.data.table(GWAS_associations)
  gwas[, chr:=as.numeric(levels(chr))[chr]]
  gwas$cptid <- gwas$rsid
  gwas$rsid <- NULL
  
  for(t in unique(coloc_res_sign[GWAS_ID==g]$tissue)){
    tmp.dt <- data.table::as.data.table(coloc_res_sign[GWAS_ID==g & tissue==t])
    if(t=="high_grade_cartilage") ab <- "hg"
    else if(t=="low_grade_cartilage") ab <- "lg"
    else ab <- t
    
    load(paste0("overlap_df_", t, "_", g, ".rda"))  # overlap_df
    for(gene in unique(tmp.dt$cpg)){
      qtl <- as.data.table(overlap_df$m_qtl[[gene]])
      merged1 <- merge(gwas, unique(qtl[, .(id=paste(paste(chr, pos, sep=":"), ALT, REF, sep="_"), rsid)]), by.x="cptid", by.y="id")
      merged2 <- merge(gwas, unique(qtl[, .(id=paste(paste(chr, pos, sep=":"), REF, ALT, sep="_"), rsid)]), by.x="cptid", by.y="id")
      
      cur.gwas <- rbind(merged1, merged2)
      
      if(nrow(tmp.dt[cpg==gene]) > 1) plot.coloc(cur.gwas, qtl, traits=c(g, paste(t, gene, sep="_")), res.coloc=tmp.dt[cpg==gene]%>% .[1], genes.data=GRCh38_Genes)
      else plot.coloc(cur.gwas, qtl, traits=c(g, paste(t, gene, sep="_")), res.coloc=tmp.dt[cpg==gene], genes.data=GRCh38_Genes)
    }
  }
}

#########################################################################
############# Check coloc.susie results
#########################################################################
coloc_susie_res <- c()
for(tissue in c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad")){
  for(GWAS_ID in names(GWAS_list)){
    tmp.dt <- fread(file.path("Coloc_susie", paste0("ColocSusie_results_", GWAS_ID, "_", tissue, ".txt")))
    tmp.dt$tissue=tissue
    tmp.dt$GWAS_ID=GWAS_ID
    coloc_susie_res <- rbind(coloc_susie_res, tmp.dt, fill=T)
  }
}
coloc_susie_res[, `:=` (cpg=sub("_.*", "", Dataset), gwas_signal=gsub(".*(rs[0-9]+).*", "\\1", Dataset))]
coloc_susie_res <- merge(coloc_susie_res, unique(signals.b38[, .(Loci, rsid)]), by.x="gwas_signal", by.y="rsid", all.x=T)
fwrite(coloc_susie_res, "all_coloc_susie_results.csv")

coloc_susie_res_sign <- coloc_susie_res[Colocalise==TRUE]
nrow(coloc_susie_res_sign)  # 15
length(unique(coloc_susie_res_sign$cpg))  # 11

# Any new genes through coloc.susie()? YES! 7
unique(coloc_susie_res_sign$cpg) %in% unique(coloc_res_sign$cpg)

# Annotate coloc_susie_res_sign with chr and pos using variant_ann
variant_ann <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/variant_ann_hg38.txt.gz")
variant_ann$chr <- as.integer(variant_ann$chr)
variant_ann <- variant_ann[!is.na(chr) & nchar(ref)==1 & nchar(alt)==1,]
coloc_susie_res_sign <- merge(coloc_susie_res_sign, variant_ann[, .(rsid, chr, pos)], by.x="hit2", by.y="rsid")
fwrite(coloc_susie_res_sign, "GO2_significant_coloc_results_susie.txt")
fwrite(coloc_susie_res_sign, "GO2_significant_coloc_results_susie.csv")

#########################################################################
############# Plot coloc.susie results
#########################################################################
for(g in unique(coloc_susie_res_sign$GWAS_ID)){
  load(paste0("GO2_GWAS_data/", g, ".rda"))  # GWAS_associations
  gwas <- as.data.table(GWAS_associations)
  gwas[, chr:=as.numeric(levels(chr))[chr]]
  gwas$cptid <- gwas$rsid
  gwas$rsid <- NULL
  
  for(t in unique(coloc_susie_res_sign[GWAS_ID==g]$tissue)){
    tmp.dt <- data.table::as.data.table(coloc_susie_res_sign[GWAS_ID==g & tissue==t])
    if(t=="high_grade_cartilage") ab <- "hg"
    else if(t=="low_grade_cartilage") ab <- "lg"
    else ab <- t
    
    load(paste0("overlap_df_", t, "_", g, ".rda"))  # overlap_df
    for(gene in unique(tmp.dt$cpg)){
      qtl <- as.data.table(overlap_df$m_qtl[[gene]])
      merged1 <- merge(gwas, unique(qtl[, .(id=paste(paste(chr, pos, sep=":"), ALT, REF, sep="_"), rsid)]), by.x="cptid", by.y="id")
      merged2 <- merge(gwas, unique(qtl[, .(id=paste(paste(chr, pos, sep=":"), REF, ALT, sep="_"), rsid)]), by.x="cptid", by.y="id")
      
      cur.gwas <- rbind(merged1, merged2)
      
      for(i in 1:nrow(tmp.dt[cpg==gene])){
        plot.coloc.susie(cur.gwas, qtl, traits=c(g, paste(t, gene, sep="_")), res.coloc=tmp.dt[cpg==gene]%>% .[i], genes.data=GRCh38_Genes)
      }
    }
  }
}

#########################################################################
############# Output colocalizing genes
#########################################################################
coloc.abf <- data.table::fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/GO2_significant_coloc_results_abf.txt")
coloc.susie <- data.table::fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/GO2_significant_coloc_results_susie.txt")
coloc.loci <-  unique(c(coloc.abf$Loci, coloc.susie$Loci))  # 74/211
coloc.signals <-  unique(c(coloc.abf$gwas_signal, coloc.susie$gwas_signal)) # 114/357
coloc.genes <- unique(c(coloc.abf$cpg, coloc.susie$cpg))  # 117

# Add gene name
mart <- biomaRt::useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 105)
coloc.genes.dt <- biomaRt::getBM(filters = "ensembl_gene_id", 
                                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                                 values = coloc.genes,
                                 mart = mart)
coloc.genes.dt <- data.table::as.data.table(coloc.genes.dt)
fwrite(coloc.genes.dt, "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/all_colocalizing_genes_with_names.csv")

output.dt <- rbind(coloc.abf[, .(gene=cpg, tissue, GWAS_ID, method="abf")], coloc.susie[,.(gene=cpg, tissue, GWAS_ID, method="susie")])
fwrite(output.dt, "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/coloc_per_gene_tissue_trait.csv")
