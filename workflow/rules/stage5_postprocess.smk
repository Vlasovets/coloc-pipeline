"""
Stage 5: Aggregate and postprocess colocalization results
Combines results across tissues and generates summary statistics
"""

def _available_coloc_results(wildcards):
    """Return only coloc result files that already exist on disk.
    Allows Stage 4 to run even when not all tissues have Stage 3 output."""
    candidates = expand(
        "{output_dir}/coloc_abf/{trait}.{tissue}.colocABF_results.txt",
        output_dir=wildcards.output_dir,
        trait=wildcards.trait,
        tissue=TISSUES,
    )
    return [f for f in candidates if os.path.exists(f)]


rule aggregate_results:
    input:
        _available_coloc_results
    output:
        all_results="{output_dir}/results/{trait}_all_coloc_results.txt",
        sig_results="{output_dir}/results/{trait}_significant_coloc_results.txt"
    params:
        pp4_threshold=config.get("pp4_threshold", 0.8),
        pp4_pp3_ratio=config.get("pp4_pp3_ratio", 2.0)
    log:
        "{output_dir}/logs/aggregate_{trait}.log"
    conda:
        "../../envs/r_coloc.yml"
    threads: 1
    resources:
        mem_mb=8000,
        time="01:00:00"
    script:
        "../scripts/5_postprocess_coloc.R"

rule generate_plots:
    input:
        results="{output_dir}/results/{trait}_significant_coloc_results.txt",
        gwas_vcf="{output_dir}/gwas_vcf/{trait}_hg38.vcf.bgz"
    output:
        plots=directory("{output_dir}/plots/{trait}")
    params:
        genes_annotation=config.get("genes_annotation_file")
    log:
        "{output_dir}/logs/plots_{trait}.log"
    conda:
        "../../envs/r_coloc.yml"
    threads: 1
    resources:
        mem_mb=8000,
        time="01:00:00"
    script:
        "../scripts/6_generate_plots.R"
