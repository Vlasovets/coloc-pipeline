"""
Stage 1: Convert GWAS summary statistics to VCF format with liftover to hg38
"""

import os

rule convert_gwas_to_vcf:
    input:
        gwas=lambda wildcards: os.path.join(
            config["gwas_sumstats_dir"],
            config["gwas_file_pattern"].format(trait=wildcards.trait)
        ),
        variant_ann=config["variant_annotation_file"]
    output:
        vcf="{output_dir}/gwas_vcf/{trait}_hg38.vcf.bgz",
        index="{output_dir}/gwas_vcf/{trait}_hg38.vcf.bgz.tbi"
    params:
        n_cases=lambda wildcards: config["sample_sizes"][wildcards.trait]["cases"],
        n_controls=lambda wildcards: config["sample_sizes"][wildcards.trait]["controls"],
        helper_script=config["helper_script"],
        genome_build=config["genome_build"],
        liftover=config["liftover_to_hg38"]
    log:
        "{output_dir}/logs/convert_gwas_{trait}.log"
    conda:
        "../../envs/r_coloc.yml"
    threads: 1
    resources:
        mem_mb=32000,
        time="04:00:00",
        tmpdir=lambda wildcards: os.path.join(config["output_dir"], "tmp")
    script:
        "../scripts/1_convert_gwas_to_vcf.R"
