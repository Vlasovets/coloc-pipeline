"""
Colocalization Pipeline
Performs colocalization analysis between GWAS and molecular QTL data
"""

import os

# Load configuration
configfile: "config.yaml"

# Define tissues and traits
TISSUES = config.get("tissues", ["high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad"])
TRAITS = config["traits"]

# Define output files
rule all:
    input:
        # Stage 1: GWAS VCF conversion
        expand(
            "{output_dir}/gwas_vcf/{trait}_hg38.vcf.gz",
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
            "{output_dir}/results/{trait}_coloc_aggregated.txt",
            output_dir=config["output_dir"],
            trait=TRAITS
        )

# Rule 1: Convert GWAS summary statistics to VCF format with liftover
rule convert_gwas_to_vcf:
    input:
        gwas=lambda wildcards: os.path.join(
            config["gwas_sumstats_dir"],
            config["gwas_file_pattern"].format(trait=wildcards.trait)
        ),
        variant_ann=config["variant_annotation_file"]
    output:
        vcf="{output_dir}/gwas_vcf/{trait}_hg38.vcf.gz",
        index="{output_dir}/gwas_vcf/{trait}_hg38.vcf.gz.tbi"
    params:
        n_cases=lambda wildcards: config["sample_sizes"][wildcards.trait]["cases"],
        n_controls=lambda wildcards: config["sample_sizes"][wildcards.trait]["controls"],
        helper_script=config["helper_script"],
        genome_build=config["genome_build"],
        liftover=config["liftover_to_hg38"]
    log:
        "{output_dir}/logs/convert_gwas_{trait}.log"
    conda:
        "envs/r_coloc.yml"
    threads: 1
    resources:
        mem_mb=32000,
        time="04:00:00",
        tmpdir=lambda wildcards: os.path.join(config["output_dir"], "tmp")
    script:
        "workflow/scripts/1_convert_gwas_to_vcf.R"

# Rule 2: Find QTL-GWAS overlaps
rule find_overlaps:
    input:
        gwas_vcf="{output_dir}/gwas_vcf/{trait}_hg38.vcf.gz",
        qtl_perm=lambda wildcards: config["qtl_permutation_files"][wildcards.tissue],
        signals=config["gwas_signals_file"]
    output:
        overlaps="{output_dir}/overlaps/{trait}.{tissue}.overlaps.rda",
        qtl_subset="{output_dir}/overlaps/{trait}.{tissue}.qtl_subset.txt"
    params:
        window_size=config.get("window_size", 1000000),
        qval_threshold=config.get("qval_threshold", 0.05)
    log:
        "{output_dir}/logs/overlaps_{trait}.{tissue}.log"
    conda:
        "envs/r_coloc.yml"
    threads: 1
    resources:
        mem_mb=16000,
        time="02:00:00"
    script:
        "workflow/scripts/2_find_qtl_gwas_overlaps.R"

# Rule 3: Run Coloc ABF
rule run_coloc_abf:
    input:
        overlaps="{output_dir}/overlaps/{trait}.{tissue}.overlaps.rda",
        qtl_data="{output_dir}/overlaps/{trait}.{tissue}.qtl_subset.txt",
        gwas_vcf="{output_dir}/gwas_vcf/{trait}_hg38.vcf.gz"
    output:
        results="{output_dir}/coloc_abf/{trait}.{tissue}.colocABF_results.txt",
        sumstats_dir=directory("{output_dir}/coloc_abf/{trait}.{tissue}.sumstats")
    params:
        qtl_sample_size=lambda wildcards: config["qtl_sample_sizes"][wildcards.tissue],
        gwas_sample_size=lambda wildcards: [
            config["sample_sizes"][wildcards.trait]["cases"],
            config["sample_sizes"][wildcards.trait]["controls"]
        ],
        cores=config.get("cores", 1)
    log:
        "{output_dir}/logs/coloc_abf_{trait}.{tissue}.log"
    conda:
        "envs/r_coloc.yml"
    threads: config.get("cores", 4)
    resources:
        mem_mb=64000,
        time="08:00:00"
    script:
        "workflow/scripts/3_run_coloc_abf.R"

# Rule 4: Run Coloc SuSiE (for ambiguous cases)
rule run_coloc_susie:
    input:
        abf_results="{output_dir}/coloc_abf/{trait}.{tissue}.colocABF_results.txt",
        sumstats_dir="{output_dir}/coloc_abf/{trait}.{tissue}.sumstats"
    output:
        results="{output_dir}/coloc_susie/{trait}.{tissue}.colocSuSiE_results.txt"
    params:
        pp4_threshold=0.25,
        ld_reference=lambda wildcards: config["ld_reference_files"][wildcards.tissue]
    log:
        "{output_dir}/logs/coloc_susie_{trait}.{tissue}.log"
    conda:
        "envs/r_coloc.yml"
    threads: 2
    resources:
        mem_mb=32000,
        time="04:00:00"
    script:
        "workflow/scripts/4_run_coloc_susie.R"

# Rule 5: Aggregate results
rule aggregate_results:
    input:
        abf_results=expand(
            "{{output_dir}}/coloc_abf/{{trait}}.{tissue}.colocABF_results.txt",
            tissue=TISSUES
        )
    output:
        aggregated="{output_dir}/results/{trait}_coloc_aggregated.txt",
        summary="{output_dir}/results/{trait}_coloc_summary.txt"
    params:
        pp4_threshold=0.8
    log:
        "{output_dir}/logs/aggregate_{trait}.log"
    conda:
        "envs/r_coloc.yml"
    threads: 1
    resources:
        mem_mb=8000,
        time="01:00:00"
    script:
        "workflow/scripts/5_postprocess_coloc.R"

# Rule 6: Generate plots (optional)
rule generate_plots:
    input:
        results="{output_dir}/results/{trait}_coloc_aggregated.txt",
        gwas_vcf="{output_dir}/gwas_vcf/{trait}_hg38.vcf.gz"
    output:
        plots=directory("{output_dir}/plots/{trait}")
    params:
        genes_annotation=config.get("genes_annotation_file")
    log:
        "{output_dir}/logs/plots_{trait}.log"
    conda:
        "envs/r_coloc.yml"
    threads: 1
    resources:
        mem_mb=8000,
        time="01:00:00"
    script:
        "workflow/scripts/6_generate_plots.R"
