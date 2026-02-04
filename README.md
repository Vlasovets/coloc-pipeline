# Colocalization Pipeline

Automated pipeline for colocalization analysis between GWAS and molecular QTL summary statistics.

## Overview

This pipeline performs two-sample Mendelian randomization and colocalization analysis between GWAS summary statistics and molecular QTL data generated from fastQTL/tensorQTL.

## Requirements

- GWAS summary statistics (hg19 or hg38)
- Table of GWAS significant independent signals
- Full molecular QTL summary statistics (fastQTL/tensorQTL format)
- Permutation results from fastQTL/tensorQTL
- Reference genotype file for LD calculations

## Pipeline Stages

1. **GWAS VCF Conversion**: Convert GWAS summary statistics to standardized VCF format with liftover (hg19→hg38)
2. **QTL-GWAS Overlap Detection**: Identify QTL-GWAS overlapping regions (±1Mb windows)
3. **Colocalization (ABF)**: Run coloc.abf for overlapping regions
4. **Colocalization (SuSiE)**: Run coloc.SuSiE for ambiguous cases (multiple causal variants)
5. **Result Aggregation**: Combine and summarize results across tissues
6. **Visualization**: Generate locus zoom plots

## Usage

### Quick Start with Automatic Checkpoints

The pipeline uses **automatic checkpoints** - it will skip completed stages and resume from the last successful step.

#### Submit to SLURM Cluster
```bash
sbatch submit_pipeline.sh
```

#### Run Locally (Interactive Mode)
```bash
bash run_pipeline_local.sh
```

### Check Pipeline Progress
```bash
bash check_pipeline_status.sh
```

This shows:
- Completed stages (✓)
- Pending stages (✗)
- Number of genes/tests per tissue
- Recent activity

### Resume After Failure

Simply rerun the same command - Snakemake will automatically:
- Skip completed steps (based on output file timestamps)
- Resume from last successful stage
- Rerun only failed/incomplete jobs

```bash
# Pipeline stopped? Just rerun:
sbatch submit_pipeline.sh
```

### Dry Run (Preview Without Execution)
```bash
snakemake --dry-run -n
```

### Force Rerun Specific Stages
```bash
# Rerun Stage 3 for all tissues
snakemake --forcerun run_coloc_abf --cores 8

# Rerun specific trait-tissue combination
snakemake results/coloc_abf/KNEE.high_grade_cartilage.colocABF_results.txt --forcerun
```

## Configuration

Edit `config.yaml` to specify:
- Input file paths (GWAS, QTL, annotations)
- Traits and tissues to analyze
- Sample sizes
- Analysis parameters (window size, thresholds)

## Output Structure

```
results/
├── gwas_vcf/              # Stage 1: Converted GWAS VCFs
├── overlaps/              # Stage 2: QTL-GWAS overlaps
├── coloc_abf/             # Stage 3: Coloc ABF results
├── coloc_susie/           # Stage 4: Coloc SuSiE results
├── results/               # Stage 5: Aggregated results
├── plots/                 # Stage 6: Visualizations
└── logs/                  # Execution logs
```

## Checkpoint System

The pipeline implements automatic checkpoints through Snakemake's DAG (Directed Acyclic Graph):

- **Output-based**: Each stage produces output files that serve as checkpoints
- **Timestamp tracking**: Snakemake compares input/output file timestamps
- **Automatic skip**: Completed stages with up-to-date outputs are skipped
- **Dependency tracking**: Stages run only when dependencies are satisfied
- **Parallelization**: Independent jobs (different tissues) run in parallel

### How It Works

1. **First run**: Executes all stages from scratch
2. **Interrupted run**: Resume picks up where it left off
3. **File modified**: Only affected downstream stages rerun
4. **Forced rerun**: Use `--forcerun` to override checkpoints

## Troubleshooting

### View logs
```bash
# Latest Snakemake log
ls -t logs/snakemake_*.out | head -1 | xargs tail

# Specific stage log
tail results/logs/coloc_abf_KNEE.high_grade_cartilage.log
```

### Clear incomplete runs
```bash
snakemake --cleanup-metadata
```

### Reset specific stages
```bash
rm results/overlaps/*.rda  # Remove Stage 2 outputs
sbatch submit_pipeline.sh  # Rerun from Stage 2
```

(C) ITG 2026
