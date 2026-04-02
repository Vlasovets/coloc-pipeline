"""
Colocalization Pipeline
Performs colocalization analysis between GWAS and molecular QTL data

Modular Snakemake pipeline with stage-specific rule files:
- Stage 1: GWAS VCF conversion (stage1_gwas_vcf.smk)
- Stage 2: QTL-GWAS overlap detection (stage2_overlaps.smk)
- Stage 3: Coloc ABF analysis (stage3_coloc_abf.smk)
- Stage 4: Coloc SuSiE fine-mapping (stage4_coloc_susie.smk)
- Stage 5: Aggregation and plotting (stage5_postprocess.smk)
"""

import os

# Load configuration
configfile: "config.yaml"

# Define tissues and traits
TISSUES = config.get("tissues", ["high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad"])
TRAITS = config["traits"]

# Import modular rules
include: "workflow/rules/stage1_gwas_vcf.smk"
include: "workflow/rules/stage2_overlaps.smk"
include: "workflow/rules/stage3_coloc_abf.smk"
include: "workflow/rules/stage4_coloc_susie.smk"
include: "workflow/rules/stage5_postprocess.smk"

# Define output files
rule all:
    input:
        # Stage 1: GWAS VCF conversion
        expand(
            "{output_dir}/gwas_vcf/{trait}_hg38.vcf.bgz",
            output_dir=config["output_dir"],
            trait=TRAITS
        ),
        # Stage 2: QTL-GWAS overlap detection
        expand(
            "{output_dir}/overlaps/{trait}.{tissue}.overlaps.rda",
            output_dir=config["output_dir"],
            trait=TRAITS,
            tissue=TISSUES
        ),
        # Stage 3: Coloc ABF results
        expand(
            "{output_dir}/coloc_abf/{trait}.{tissue}.colocABF_results.txt",
            output_dir=config["output_dir"],
            trait=TRAITS,
            tissue=TISSUES
        ),
        # Stage 4: Aggregated results
        expand(
            "{output_dir}/results/{trait}_all_coloc_results.txt",
            output_dir=config["output_dir"],
            trait=TRAITS
        )
