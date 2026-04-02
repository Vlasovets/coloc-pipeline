"""
Stage 2: Find QTL-GWAS overlaps
Identifies genomic regions where GWAS signals overlap with eQTLs
"""

rule find_overlaps:
    input:
        gwas_vcf="{output_dir}/gwas_vcf/{trait}_hg38.vcf.bgz",
        qtl_perm=lambda wildcards: config["qtl_permutation_files"][wildcards.tissue],
        signals=config["gwas_signals_file"]
    output:
        overlaps="{output_dir}/overlaps/{trait}.{tissue}.overlaps.rda",
        qtl_subset="{output_dir}/overlaps/{trait}.{tissue}.qtl_subset.txt",
        gwas_data="{output_dir}/overlaps/{trait}.{tissue}.gwas_data.rda"
    params:
        window_size=config.get("window_size", 1000000),
        qval_threshold=config.get("qval_threshold", 0.05),
        variant_ann=config["variant_annotation_file"]
    log:
        "{output_dir}/logs/overlaps_{trait}.{tissue}.log"
    conda:
        "../../envs/r_coloc.yml"
    threads: 1
    resources:
        mem_mb=16000,
        time="02:00:00"
    script:
        "../scripts/2_find_qtl_gwas_overlaps.R"
