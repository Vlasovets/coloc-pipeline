"""
Stage 4: Run Coloc SuSiE (Sum of Single Effects) analysis
Performs fine-mapping for ambiguous colocalization cases
"""

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
        "../../envs/r_coloc.yml"
    threads: 2
    resources:
        mem_mb=32000,
        time="04:00:00"
    script:
        "../scripts/4_run_coloc_susie.R"
