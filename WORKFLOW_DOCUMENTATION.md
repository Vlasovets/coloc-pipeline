# Coloc Pipeline Workflow Documentation

## Overview

This document describes the refactored colocalization pipeline with modular, production-ready scripts.

## File Structure

```
coloc-pipeline/
├── Snakefile                    # Main workflow orchestration
├── config.yaml                  # Configuration file
├── submit_pipeline.sh           # SLURM submission script
├── environment.yml              # Conda environment
│
├── workflow/
│   ├── scripts/                 # Modular analysis scripts
│   │   ├── data_loader.R       # Data loading functions
│   │   ├── qtl_processor.R     # QTL-GWAS overlap detection
│   │   ├── coloc_helpers.R     # Helper functions (existing)
│   │   ├── convert_gwas_to_vcf.R  # GWAS VCF conversion (existing)
│   │   ├── run_coloc_abf.R     # Coloc ABF analysis
│   │   ├── run_coloc_susie.R   # Coloc SuSiE analysis
│   │   ├── postprocess_coloc.R # Result aggregation
│   │   └── generate_plots.R    # Visualization (to be created)
│   │
│   └── rules/                   # Additional Snakemake rules (optional)
│
├── src/                         # Original reference code (NOT USED)
│   ├── README.md
│   ├── Coloc_helper_functions.R
│   ├── 1_make_GO2_b38_vcf.R
│   ├── 3_run_coloc_abf.R
│   ├── 4_run_coloc_susie.R
│   └── 5_postprocess_coloc.R
│
├── results/                     # Output directory
│   ├── gwas_vcf/               # Converted GWAS VCFs
│   ├── overlaps/               # QTL-GWAS overlaps
│   ├── coloc_abf/              # Coloc ABF results
│   ├── coloc_susie/            # Coloc SuSiE results
│   ├── results/                # Aggregated results
│   ├── plots/                  # Visualization outputs
│   └── logs/                   # Execution logs
│
└── logs/                        # Snakemake logs

```

## Pipeline Stages

### Stage 1: GWAS VCF Conversion
**Rule:** `convert_gwas_to_vcf`

Converts GWAS summary statistics to standardized VCF format with genome build liftover.

**Inputs:**
- GWAS summary statistics (text/gz format)
- Variant annotation file

**Outputs:**
- `results/gwas_vcf/{trait}_hg38.vcf.gz`
- `results/gwas_vcf/{trait}_hg38.vcf.gz.tbi`

**Script:** `workflow/scripts/convert_gwas_to_vcf.R`

### Stage 2: QTL-GWAS Overlap Detection
**Rule:** `find_overlaps`

Identifies QTL genes that overlap with GWAS signal regions.

**Inputs:**
- GWAS VCF file
- QTL permutation results
- GWAS independent signals

**Outputs:**
- `results/overlaps/{trait}_{tissue}_overlaps.rda`
- `results/overlaps/{trait}_{tissue}_qtl_subset.txt`

**Script:** `workflow/scripts/find_qtl_gwas_overlaps.R` (to be created using qtl_processor.R)

### Stage 3: Colocalization Analysis (ABF)
**Rule:** `run_coloc_abf`

Performs Bayesian colocalization using coloc.abf method.

**Inputs:**
- QTL-GWAS overlaps
- QTL nominal pass data
- GWAS VCF

**Outputs:**
- `results/coloc_abf/{trait}_{tissue}_colocABF_results.txt`
- `results/coloc_abf/{trait}_{tissue}_sumstats/` (input files for plotting)

**Script:** `workflow/scripts/run_coloc_abf.R`

### Stage 4: Colocalization with SuSiE (Optional)
**Rule:** `run_coloc_susie`

Tests for colocalization with multiple causal variants using SuSiE.

**Inputs:**
- Coloc ABF results (filters for PP4 > 0.25)
- Summary statistics
- LD matrices

**Outputs:**
- `results/coloc_susie/{trait}_{tissue}_colocSuSiE_results.txt`

**Script:** `workflow/scripts/run_coloc_susie.R`

### Stage 5: Result Aggregation
**Rule:** `aggregate_results`

Combines results across tissues and generates summary statistics.

**Inputs:**
- All coloc ABF results for a trait

**Outputs:**
- `results/results/{trait}_coloc_aggregated.txt`
- `results/results/{trait}_coloc_summary.txt`

**Script:** `workflow/scripts/postprocess_coloc.R`

### Stage 6: Visualization (Optional)
**Rule:** `generate_plots`

Creates locus zoom plots for significant colocalizations.

**Inputs:**
- Aggregated results
- GWAS VCF

**Outputs:**
- `results/plots/{trait}/` (directory with plot files)

**Script:** `workflow/scripts/generate_plots.R` (to be created)

## Modular Scripts Overview

### 1. data_loader.R
**Purpose:** Centralized data loading with validation

**Key Functions:**
- `load_variant_annotation()` - Load and filter variant annotations
- `load_gwas_vcf()` - Extract GWAS data from VCF
- `extract_gwas_windows()` - Get GWAS data for multiple windows
- `load_qtl_permutations()` - Load significant QTL results
- `load_qtl_nominal()` - Load QTL nominal pass data
- `load_gwas_signals()` - Load independent GWAS signals

**Improvements from original:**
- Error handling with informative messages
- Data validation
- Consistent logging
- Memory-efficient loading

### 2. qtl_processor.R
**Purpose:** QTL-GWAS data harmonization

**Key Functions:**
- `create_gwas_windows()` - Create genomic windows around signals
- `create_qtl_windows()` - Create windows around QTL variants
- `find_qtl_gwas_overlaps()` - Detect overlapping regions
- `harmonize_qtl_gwas()` - Allele harmonization and filtering
- `prepare_coloc_qtl()` - Format data for coloc
- `prepare_coloc_gwas()` - Format GWAS for coloc

**Improvements from original:**
- Cleaner separation of concerns
- Better memory management
- Reproducible window creation
- Extensive validation

### 3. run_coloc_abf.R
**Purpose:** Execute coloc.abf analysis

**Key Functions:**
- `run_coloc_abf()` - Single coloc test with error handling
- `extract_credible_set()` - Get 95% credible set
- CLI interface for Snakemake integration

**Improvements from original:**
- Modular function design
- Progress tracking
- Parallel processing support
- Error recovery

### 4. run_coloc_susie.R
**Purpose:** Execute coloc.susie for complex loci

**Key Functions:**
- `run_coloc_susie()` - SuSiE colocalization
- LD matrix preparation
- CLI interface

**Improvements from original:**
- Cleaner LD handling
- Better error messages
- MHC region detection

### 5. postprocess_coloc.R
**Purpose:** Aggregate and summarize results

**Key Functions:**
- `aggregate_coloc_results()` - Combine all results
- `summarize_results()` - Generate summary statistics
- Filtering by PP4 threshold

**Improvements from original:**
- Flexible aggregation
- Summary statistics
- Export formats

## Execution Instructions

### 1. Configuration

Edit `config.yaml` to specify:
```yaml
# Data paths
gwas_sumstats_dir: "/path/to/gwas"
variant_annotation_file: "/path/to/variants.txt.gz"
gwas_signals_file: "/path/to/signals.csv"

# Tissues to analyze
tissues:
  - high_grade_cartilage
  - low_grade_cartilage
  - synovium
  - fat_pad

# Traits
traits:
  - KNEE
  - TKR
  - ALLOA

# Sample sizes
sample_sizes:
  KNEE:
    cases: 172256
    controls: 1144244

# QTL data
qtl_permutation_files:
  high_grade_cartilage: "/path/to/hg_perm.txt"

qtl_sample_sizes:
  high_grade_cartilage: 115

# Analysis parameters
window_size: 1000000
qval_threshold: 0.05
cores: 4
```

### 2. Dry Run

Test the workflow without execution:
```bash
snakemake --dry-run --cores 1
```

### 3. Local Execution

Run on local machine:
```bash
snakemake --cores 4 --use-conda
```

### 4. SLURM Cluster Execution

Submit to cluster:
```bash
sbatch submit_pipeline.sh
```

Or run interactively with cluster profile:
```bash
snakemake --profile profiles/slurm --cores 100
```

### 5. Specific Rules

Run only GWAS VCF conversion:
```bash
snakemake --cores 1 results/gwas_vcf/KNEE_hg38.vcf.gz
```

Run coloc for specific tissue:
```bash
snakemake --cores 4 results/coloc_abf/KNEE_high_grade_cartilage_colocABF_results.txt
```

### 6. Resume Failed Jobs

Snakemake automatically detects incomplete outputs:
```bash
snakemake --cores 4 --use-conda --rerun-incomplete
```

## Key Improvements from Original Code

### Code Organization
- **Before:** Monolithic scripts with hardcoded paths
- **After:** Modular functions with clear responsibilities

### Configuration
- **Before:** Paths hardcoded in scripts
- **After:** Centralized config.yaml

### Error Handling
- **Before:** Minimal error checking
- **After:** Try-catch blocks, validation, informative errors

### Logging
- **Before:** print() statements
- **After:** Timestamped logging with levels

### Reproducibility
- **Before:** Manual execution order
- **After:** Snakemake dependency management

### Scalability
- **Before:** Single-threaded execution
- **After:** Parallel processing, cluster support

### Documentation
- **Before:** Comments in code
- **After:** Docstrings, type hints, README

## Troubleshooting

### Issue: Missing input files
**Solution:** Check config.yaml paths and run `snakemake --list-input-files`

### Issue: Memory errors
**Solution:** Reduce cores, increase memory in submit_pipeline.sh

### Issue: Coloc errors
**Solution:** Check logs in `results/logs/`, validate input data formats

### Issue: Missing packages
**Solution:** Ensure conda environment is activated: `conda activate coloc-pipeline`

## Contact & Support

For questions about the refactored pipeline:
- Check logs in `results/logs/`
- Review original code in `src/` for algorithm details

---

**Last Updated:** 2026-02-04
**Pipeline Version:** 2.0 (Modular Refactor)
