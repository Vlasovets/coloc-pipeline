# Change allele2 to be also a SNP (%in% c("A", "T", "C", "G"))
# Use CPTID instead of SNP

.libPaths(c("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/R", .libPaths()))
project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc"
data.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/"
setwd(project.path)

source("scripts/Coloc_helper_functions.R")

library(data.table)
library(VariantAnnotation)
library(gwasvcf)
library(magrittr)
library("GenomicRanges")
library(rtracklayer)
library(coloc)
library(dplyr)
library(parallel)

variant_ann <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/variant_ann_hg38.txt.gz")
variant_ann$chr <- as.integer(variant_ann$chr)
variant_ann <- variant_ann[!is.na(variant_ann$chr) & nchar(variant_ann$ref)==1 & nchar(variant_ann$alt)==1,]

for(trait in c("KNEE", "TKR", "ALLOA")){
  ncases <- ifelse(trait=="KNEE", 172256, ifelse(trait=="TKR", 48161, 489952))
  ncontrols <- ifelse(trait=="KNEE", 1144244, ifelse(trait=="TKR", 958463, 1471094))
  GWAS_n <- c(ncases, ncontrols)
  GWASfile <- paste0("/lustre/groups/itg/teams/zeggini/projects/GO2/GO2SummStats/ALL.MAIN/GO.FILTER.GW.final.meta.results.ALL.", trait, ".FULL.MAFless0.01.Nupdated.txt.gz")
  make_vcf(GWASfile=GWASfile,
           chrom="CHR",
           pos="POS",
           nea="NEA",
           ea="EA",
           snp="CPTID",
           ea_af="EAF",
           effect="BETA",
           se="SE",
           pval="P",
           hg="hg19", # hg build of summary stats !! For GO2 it's hg19 so change here and then change the want to lift over here (not sure where the chain file has to be)
           lift_down=FALSE, # if lift_down=TRUE then hg38 is converted to hg19. If FALSE, hg19 is converted to hg38. Ignored if WantToLiftOver is set ot FALSE
           WantToLiftOver=TRUE, # Whether or not you want to lift over the coordinates
           GWAS_n=GWAS_n, # a vector of one or two elements. If quantitative trait, total sample size, if case control, number of cases and controls
           variant_ann=variant_ann, # reference file to map the missing values
           # output=paste0("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/GO2_sumstats/GO2_b38_", trait, "_ody_cptid.vcf"))
           output=paste0("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/GO2_sumstats/GO2_b38_", trait, "_ana_new.vcf"))
}
