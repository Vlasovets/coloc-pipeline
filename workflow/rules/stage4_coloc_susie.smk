"""
Stage 5: Run Coloc SuSiE (Sum of Single Effects) analysis
Performs fine-mapping for ambiguous colocalization cases (PP4 > 0.25, not already significant).
Requires LD reference: tissue-specific plink bfiles + UKB chromosome bfiles for GWAS.
"""

rule run_coloc_susie:
    input:
        abf_results="{output_dir}/coloc_abf/{trait}.{tissue}.colocABF_results.txt",
        gwas_data="{output_dir}/overlaps/{trait}.{tissue}.gwas_data.rda",
        qtl_subset="{output_dir}/overlaps/{trait}.{tissue}.qtl_subset.txt"
    output:
        results="{output_dir}/coloc_susie/{trait}.{tissue}.colocSuSiE_results.txt",
        ld_dir=directory("{output_dir}/coloc_susie/LD/{trait}.{tissue}")
    params:
        gwas_id=lambda wildcards: wildcards.trait,
        tissue=lambda wildcards: wildcards.tissue,
        qtl_n=lambda wildcards: config["qtl_sample_sizes"][wildcards.tissue],
        gwas_n=lambda wildcards: [
            config["sample_sizes"][wildcards.trait]["cases"],
            config["sample_sizes"][wildcards.trait]["controls"]
        ],
        variant_ann=config["variant_annotation_file"],
        qtl_nominal_dir=lambda wildcards: config["qtl_nominal_dirs"][wildcards.tissue],
        qtl_bfile=lambda wildcards: config["qtl_genotype_bfiles"][wildcards.tissue],
        gwas_bfile_prefix=config["gwas_ld_bfile_prefix"],
        pp4_threshold=config.get("susie_pp4_threshold", 0.25),
        pp4_pp3_ratio=config.get("pp4_pp3_ratio", 2.0),
        plink_bin=config.get("plink_bin", "/lustre/groups/itg/shared/software/bin/plink")
    log:
        "{output_dir}/logs/coloc_susie_{trait}.{tissue}.log"
    conda:
        "../../envs/r_coloc.yml"
    threads: 2
    resources:
        mem_mb=32000,
        time="08:00:00"
    script:
        "../scripts/4_run_coloc_susie.R"
