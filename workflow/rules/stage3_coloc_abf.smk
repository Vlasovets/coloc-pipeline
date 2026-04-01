"""
Stage 3: Run Coloc ABF (Approximate Bayes Factor) analysis
Performs colocalization testing using Bayesian approach
"""

rule run_coloc_abf:
    input:
        overlaps="{output_dir}/overlaps/{trait}.{tissue}.overlaps.rda",
        qtl_subset="{output_dir}/overlaps/{trait}.{tissue}.qtl_subset.txt",
        gwas_vcf="{output_dir}/gwas_vcf/{trait}_hg38.vcf.bgz",
        gwas_data="{output_dir}/overlaps/{trait}.{tissue}.gwas_data.rda"
    output:
        "{output_dir}/coloc_abf/{trait}.{tissue}.colocABF_results.txt"
    params:
        gwas_id=lambda wildcards: wildcards.trait,
        tissue=lambda wildcards: wildcards.tissue,
        qtl_n=lambda wildcards: config["qtl_sample_sizes"][wildcards.tissue],
        gwas_n=lambda wildcards: [
            config["sample_sizes"][wildcards.trait]["cases"],
            config["sample_sizes"][wildcards.trait]["controls"]
        ],
        variant_ann=config["variant_annotation_file"],
        qtl_nominal_dir=lambda wildcards: config["qtl_nominal_dirs"][wildcards.tissue]
    log:
        "{output_dir}/logs/coloc_abf_{trait}.{tissue}.log"
    conda:
        "../../envs/r_coloc.yml"
    threads: config.get("cores", 4)
    resources:
        mem_mb=64000,
        time="08:00:00"
    script:
        "../scripts/3_run_coloc_abf.R"
