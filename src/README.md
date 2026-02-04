# Source Code Archive

This directory contains the original working version of the colocalization pipeline scripts from the MTutino/Coloc repository. These files are kept as reference but are **NOT used** in the production Snakemake pipeline.

## Reference Implementation

The scripts in this directory represent the original research code that formed the basis for this pipeline. They demonstrate the full workflow but contain hardcoded paths and are not designed for reusability.

## Production Pipeline

The production pipeline with best practices is located in:
- `workflow/scripts/` - Clean, reusable scripts
- `Snakefile` - Main workflow orchestration
- `config.yaml` - Configuration without hardcoded paths

## Files

These files should be copied from: https://github.com/MTutino/Coloc/tree/main/coloc

- `Coloc_helper_functions.R` - Core helper functions
- `1_make_GO2_b38_vcf.R` - GWAS VCF conversion
- `3_run_coloc_abf.R` - Coloc ABF analysis
- `4_run_coloc_susie.R` - Coloc SuSiE analysis
- `5_postprocess_coloc.R` - Post-processing and plotting
- `plot_functions.R` - Plotting utilities
- Jupyter notebooks (`.ipynb`) - Interactive analysis examples

## Usage

**Do not run these scripts directly**. They contain absolute paths specific to the original computing environment. Use them only as reference for understanding the analysis logic.

For running analyses, use the Snakemake pipeline in the project root directory.
