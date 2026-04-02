# Pipeline Scripts Documentation

## Overview

The `scripts/` directory contains modular SLURM batch scripts for running the colocalization pipeline. Each stage can be run independently or as part of the complete pipeline.

## Directory Structure

```
scripts/
├── submit_pipeline.sh      # Main pipeline orchestrator
├── stage1_gwas_vcf.sh     # Stage 1: GWAS VCF conversion
├── stage2_overlaps.sh     # Stage 2: QTL-GWAS overlap detection  
├── stage3_coloc_abf.sh    # Stage 3: Coloc ABF analysis
├── stage4_aggregate.sh    # Stage 4: Aggregate results
├── check_status.sh        # Pipeline status checker
├── clean.sh               # Cleanup utility
└── QUICKSTART.sh          # Quick reference guide
```

## Usage

### Complete Pipeline

Run all stages sequentially:
```bash
# Sequential mode (stages run one after another)
sbatch scripts/submit_pipeline.sh KNEE sequential

# Snakemake mode (parallelized across tissues)
sbatch scripts/submit_pipeline.sh KNEE snakemake
```

### Individual Stages

Run specific stages independently:

**Stage 1: GWAS VCF Conversion**
```bash
sbatch scripts/stage1_gwas_vcf.sh KNEE
```

**Stage 2: QTL-GWAS Overlaps** 
```bash
sbatch scripts/stage2_overlaps.sh KNEE "high_grade_cartilage low_grade_cartilage"
```

**Stage 3: Coloc ABF Analysis**
```bash
sbatch scripts/stage3_coloc_abf.sh KNEE "high_grade_cartilage low_grade_cartilage"
```

**Stage 4: Aggregate Results**
```bash
sbatch scripts/stage4_aggregate.sh KNEE
```

### Utilities

**Check Pipeline Status**
```bash
bash scripts/check_status.sh KNEE
```

**Clean Intermediate Files**
```bash
bash scripts/clean.sh KNEE clean
```

**Purge All Outputs**
```bash
bash scripts/clean.sh KNEE purge
```

## Script Details

### submit_pipeline.sh

Main orchestrator that runs all pipeline stages.

**Arguments:**
- `TRAIT`: Trait to analyze (default: KNEE)
- `MODE`: Execution mode - `sequential` or `snakemake` (default: sequential)

**Features:**
- ✓ Dependency checking between stages
- ✓ Error handling and logging
- ✓ Two execution modes (sequential or parallelized)
- ✓ Automatic directory creation
- ✓ Summary reporting

**SLURM Resources:**
- Time: 24 hours
- Memory: 8G
- CPUs: 1

### stage1_gwas_vcf.sh

Converts GWAS summary statistics to VCF format with liftover to hg38.

**Arguments:**
- `TRAIT`: Trait to process (default: KNEE)
- `DRY_RUN`: Set to "true" for dry-run (default: false)

**SLURM Resources:**
- Time: 4 hours
- Memory: 32G
- CPUs: 1

**Outputs:**
- `${OUTPUT_DIR}/gwas_vcf/${TRAIT}_hg38.vcf.bgz`
- `${OUTPUT_DIR}/gwas_vcf/${TRAIT}_hg38.vcf.bgz.tbi`

### stage2_overlaps.sh

Identifies genomic regions where GWAS signals overlap with eQTLs.

**Arguments:**
- `TRAIT`: Trait to process (default: KNEE)
- `TISSUES`: Space-separated tissue list (default: all 4 tissues)
- `DRY_RUN`: Set to "true" for dry-run (default: false)

**SLURM Resources:**
- Time: 2 hours
- Memory: 16G
- CPUs: 4

**Outputs:**
- `${OUTPUT_DIR}/overlaps/${TRAIT}.${TISSUE}.overlaps.rda`
- `${OUTPUT_DIR}/overlaps/${TRAIT}.${TISSUE}.qtl_subset.txt`

### stage3_coloc_abf.sh

Performs colocalization testing using Approximate Bayes Factor.

**Arguments:**
- `TRAIT`: Trait to process (default: KNEE)
- `TISSUES`: Space-separated tissue list (default: all 4 tissues)
- `DRY_RUN`: Set to "true" for dry-run (default: false)

**SLURM Resources:**
- Time: 8 hours
- Memory: 64G
- CPUs: 4

**Outputs:**
- `${OUTPUT_DIR}/coloc_abf/${TRAIT}.${TISSUE}.colocABF_results.txt`
- `${OUTPUT_DIR}/coloc_abf/${TRAIT}.${TISSUE}.sumstats/`

### stage4_aggregate.sh

Combines colocalization results across tissues and generates summary.

**Arguments:**
- `TRAIT`: Trait to process (default: KNEE)
- `DRY_RUN`: Set to "true" for dry-run (default: false)

**SLURM Resources:**
- Time: 1 hour
- Memory: 8G
- CPUs: 1

**Outputs:**
- `${OUTPUT_DIR}/results/${TRAIT}_coloc_aggregated.txt`
- `${OUTPUT_DIR}/results/${TRAIT}_coloc_summary.txt`

### check_status.sh

Checks completion status of all pipeline stages.

**Usage:**
```bash
bash scripts/check_status.sh [TRAIT]
```

**Features:**
- ✓ Color-coded output
- ✓ File size and modification time
- ✓ Gene/result counts
- ✓ Overall completion status

### clean.sh

Removes intermediate files or purges all outputs.

**Usage:**
```bash
bash scripts/clean.sh [TRAIT] [MODE]
```

**Modes:**
- `clean`: Removes intermediate files (.snakemake, tmp, old logs)
- `purge`: Removes ALL outputs including final results (prompts for confirmation)

**Safety:**
- Prompts before purging
- Can target specific trait or all traits

## Examples

### Complete Workflow

```bash
# 1. Check current status
bash scripts/check_status.sh KNEE

# 2. Run complete pipeline
sbatch scripts/submit_pipeline.sh KNEE sequential

# 3. Monitor progress
squeue -u $USER

# 4. Check status again
bash scripts/check_status.sh KNEE
```

### Targeted Re-run

```bash
# Re-run only Stage 3 for specific tissues
sbatch scripts/stage3_coloc_abf.sh KNEE "high_grade_cartilage synovium"
```

### Testing/Debugging

```bash
# Dry-run Stage 2
bash scripts/stage2_overlaps.sh KNEE "high_grade_cartilage" true

# Clean and restart
bash scripts/clean.sh KNEE purge
sbatch scripts/submit_pipeline.sh KNEE sequential
```

## Error Handling

All scripts include:
- **Exit on error** (`set -euo pipefail`)
- **Prerequisite checking** (verifies previous stages completed)
- **Detailed logging** (timestamps, stage names, file paths)
- **Exit codes** (0 = success, non-zero = failure)

## Logging

Logs are written to:
- **SLURM logs**: `logs/stageN_${JOBID}.{out,err}`
- **Pipeline logs**: `${OUTPUT_DIR}/logs/${STAGE}_${TRAIT}.${TISSUE}.log`

## Best Practices

1. **Always check status first**: `bash scripts/check_status.sh TRAIT`
2. **Use dry-run for testing**: Add `true` as last argument
3. **Monitor SLURM jobs**: `squeue -u $USER`
4. **Check logs on failure**: `tail logs/stage*_*.err`
5. **Clean before full re-run**: `bash scripts/clean.sh TRAIT purge`

## Troubleshooting

**Pipeline fails at Stage 1:**
- Check input GWAS file exists
- Verify variant annotation file path
- Check conda environment is activated

**Pipeline fails at Stage 2:**
- Verify Stage 1 completed successfully
- Check QTL permutation files exist
- Verify GWAS signals file path

**Pipeline fails at Stage 3:**
- Ensure Stage 2 completed for all tissues
- Check sufficient memory allocation
- Verify QTL nominal directories exist

**Out of memory errors:**
- Increase `--mem` in SLURM directives
- Reduce number of parallel cores
- Process fewer tissues at once

## Advanced Usage

### Custom Tissue Subsets

```bash
# Only high and low grade cartilage
export TISSUES="high_grade_cartilage low_grade_cartilage"
sbatch scripts/stage2_overlaps.sh KNEE "$TISSUES"
```

### Job Dependencies

```bash
# Chain jobs with dependencies
JOB1=$(sbatch --parsable scripts/stage1_gwas_vcf.sh KNEE)
JOB2=$(sbatch --dependency=afterok:$JOB1 --parsable scripts/stage2_overlaps.sh KNEE)
JOB3=$(sbatch --dependency=afterok:$JOB2 scripts/stage3_coloc_abf.sh KNEE)
```

### Resource Adjustment

Edit SLURM directives in individual scripts:
```bash
#SBATCH --time=12:00:00  # Increase time limit
#SBATCH --mem=128G       # Increase memory
#SBATCH --cpus-per-task=8 # More CPUs
```
