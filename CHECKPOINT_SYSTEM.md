# Pipeline Checkpoint System

## Overview

The colocalization pipeline now has **automatic checkpoint functionality** that allows it to:
- Resume from the last successful stage after interruption
- Skip completed stages when rerunning
- Parallelize independent jobs across tissues and traits
- Track dependencies between stages

## How to Run the Complete Pipeline

### Option 1: Submit to SLURM Cluster (Recommended for Production)

```bash
cd /home/itg/oleg.vlasovets/projects/coloc-pipeline
sbatch submit_pipeline.sh
```

**Features:**
- Automatically submits jobs to SLURM cluster
- Parallelizes up to 20 jobs simultaneously
- Each rule gets appropriate resources (memory, time, CPUs)
- Automatic checkpointing built-in

**Configuration in `submit_pipeline.sh`:**
- Main Snakemake job: 32GB RAM, 24 hours, 8 CPUs
- Individual rules: Configured via Snakefile `resources` directives
- Cluster submission via `--cluster` flag
- Job monitoring via `slurm_status.py`

### Option 2: Run Locally (For Testing or Interactive Sessions)

```bash
cd /home/itg/oleg.vlasovets/projects/coloc-pipeline
bash run_pipeline_local.sh
```

**Features:**
- Runs on current machine (interactive node)
- Uses 8 cores by default
- Simpler execution without cluster overhead
- Still has automatic checkpointing

## Checkpoint Mechanism

### How Snakemake Checkpoints Work

1. **DAG Construction**: Snakemake builds a Directed Acyclic Graph of all tasks
2. **Output File Checking**: Before running each rule, checks if output files exist
3. **Timestamp Comparison**: Compares timestamps of input vs output files
4. **Smart Execution**: 
   - ✓ Output exists + newer than inputs → **SKIP** (checkpoint)
   - ✗ Output missing or older than inputs → **RUN**
   - ⚠ Upstream dependency changed → **RERUN** affected stages

### Pipeline Stages with Checkpoints

```
Stage 1: GWAS VCF Conversion
├─ Checkpoint: results/gwas_vcf/{trait}_hg38.vcf.gz
├─ Dependencies: GWAS summary stats, variant annotation
└─ Triggers Stage 2 for all tissues

Stage 2: QTL-GWAS Overlap Detection  
├─ Checkpoint: results/overlaps/{trait}.{tissue}.overlaps.rda
├─ Dependencies: Stage 1 VCF, QTL permutations
└─ Triggers Stage 3 for same tissue

Stage 3: Coloc ABF Analysis
├─ Checkpoint: results/coloc_abf/{trait}.{tissue}.colocABF_results.txt
├─ Dependencies: Stage 2 overlaps, QTL nominal data
└─ Triggers Stage 5 aggregation

Stage 4: Coloc SuSiE (Optional, conditional)
├─ Checkpoint: results/coloc_susie/{trait}.{tissue}.colocSuSiE_results.txt
├─ Dependencies: Stage 3 ABF results (PP.H4 < 0.25)
└─ Used for ambiguous colocalization cases

Stage 5: Result Aggregation
├─ Checkpoint: results/results/{trait}_coloc_aggregated.txt
├─ Dependencies: All Stage 3 results for all tissues
└─ Triggers Stage 6 plots

Stage 6: Visualization
├─ Checkpoint: results/plots/{trait}/
├─ Dependencies: Stage 5 aggregated results
└─ Final outputs
```

## Key Features

### 1. Automatic Resume After Failure

If pipeline fails at any stage, simply rerun the same command:

```bash
# Pipeline stopped at Stage 3? Just rerun:
sbatch submit_pipeline.sh

# Snakemake automatically:
# - Skips Stage 1 (already completed)
# - Skips Stage 2 (already completed)
# - Resumes Stage 3 from where it failed
```

### 2. Parallel Execution

Independent jobs run in parallel:

```
Stage 2: Find Overlaps (12 jobs run in parallel)
├─ KNEE + high_grade_cartilage
├─ KNEE + low_grade_cartilage
├─ KNEE + synovium
├─ KNEE + fat_pad
├─ TKR + high_grade_cartilage
├─ ... (all trait-tissue combinations)
└─ ALLOA + fat_pad
```

### 3. Smart Dependency Tracking

If you modify an input file, only affected stages rerun:

```bash
# Modify GWAS signals file
vim GO2_index_signals_b38.csv

# Rerun pipeline
sbatch submit_pipeline.sh

# Result: Only Stage 2-6 rerun, Stage 1 is skipped (unchanged)
```

### 4. Rerun Control

Force rerun specific stages:

```bash
# Rerun only Stage 3 for all tissues
snakemake --forcerun run_coloc_abf --cores 8

# Rerun everything for KNEE trait
snakemake results/results/KNEE_coloc_aggregated.txt --forcerun

# Rerun specific tissue only
snakemake results/coloc_abf/KNEE.high_grade_cartilage.colocABF_results.txt --forcerun
```

## Monitoring Progress

### Check Pipeline Status

```bash
bash check_pipeline_status.sh
```

**Output shows:**
- ✓ Completed stages with file sizes
- ✗ Pending stages
- Number of genes/tests processed
- Recent log activity

### View Execution Plan (Dry Run)

```bash
snakemake --dry-run -n
```

Shows what will run without actually executing.

### Monitor SLURM Jobs

```bash
# Check running jobs
squeue -u $USER

# Check completed jobs
sacct -j <JOBID> --format=JobID,JobName,State,ExitCode,Elapsed
```

### View Logs

```bash
# Main Snakemake log
tail -f logs/snakemake_<JOBID>.out

# Individual stage logs
tail -f results/logs/coloc_abf_KNEE.high_grade_cartilage.log

# SLURM rule logs
ls logs/slurm_*_run_coloc_abf_*.out
```

## Resource Configuration

Each rule has specific resource requirements set in `Snakefile`:

| Stage | Memory | Time | CPUs |
|-------|--------|------|------|
| GWAS VCF Conversion | 32 GB | 4h | 1 |
| Find Overlaps | 16 GB | 2h | 1 |
| Coloc ABF | 64 GB | 8h | 4 |
| Coloc SuSiE | 32 GB | 4h | 2 |
| Aggregation | 8 GB | 1h | 1 |
| Plots | 8 GB | 1h | 1 |

These can be adjusted in the Snakefile `resources` section of each rule.

## Troubleshooting

### Pipeline Won't Start

```bash
# Check Snakemake can see all inputs
snakemake --summary

# Check for syntax errors
snakemake --lint
```

### Jobs Keep Rerunning (Not Checkpointing)

Possible causes:
1. Output files deleted manually
2. Input files have newer timestamps than outputs
3. Rule definition changed in Snakefile

```bash
# Force rebuild DAG
snakemake --forceall --dry-run

# Check file timestamps
ls -lt results/overlaps/
```

### Incomplete Job Cleanup

If jobs were killed ungracefully:

```bash
# Mark incomplete jobs for rerun
snakemake --rerun-incomplete

# Clean up incomplete outputs
snakemake --cleanup-metadata
```

### Reset Specific Stages

```bash
# Remove Stage 2 outputs to force rerun
rm -rf results/overlaps/*.rda results/overlaps/*.txt

# Remove all outputs to start fresh
rm -rf results/gwas_vcf results/overlaps results/coloc_abf results/results

# Then rerun
sbatch submit_pipeline.sh
```

## Advanced Usage

### Run Specific Trait Only

```bash
# Edit config.yaml to specify single trait
traits: ["KNEE"]

# Run pipeline
sbatch submit_pipeline.sh
```

### Run Specific Stage Range

```bash
# Run only Stage 1 and 2
snakemake results/overlaps/*.overlaps.rda --cores 8

# Run only Stage 3 onwards (assuming 1-2 completed)
snakemake results/results/*_coloc_aggregated.txt --cores 8
```

### Custom Checkpoint File

Create manual checkpoint markers:

```bash
# After Stage 2 completes
touch .checkpoint_stage2_complete

# Check in scripts
if [ -f .checkpoint_stage2_complete ]; then
    echo "Stage 2 already done, skipping..."
fi
```

## Files Overview

```
coloc-pipeline/
├── submit_pipeline.sh              # SLURM cluster execution
├── run_pipeline_local.sh           # Local execution
├── check_pipeline_status.sh        # Progress monitoring
├── Snakefile                       # Pipeline definition with checkpoints
├── config.yaml                     # Configuration
├── workflow/scripts/
│   ├── slurm_status.py            # SLURM job status checker
│   ├── 1_convert_gwas_to_vcf.R    # Stage 1
│   ├── 2_find_qtl_gwas_overlaps.R # Stage 2
│   ├── 3_run_coloc_abf.R          # Stage 3
│   ├── 4_run_coloc_susie.R        # Stage 4
│   ├── 5_postprocess_coloc.R      # Stage 5
│   └── 6_generate_plots.R         # Stage 6
└── results/                        # Checkpoint outputs
    ├── gwas_vcf/                   # Stage 1 checkpoints
    ├── overlaps/                   # Stage 2 checkpoints
    ├── coloc_abf/                  # Stage 3 checkpoints
    ├── coloc_susie/                # Stage 4 checkpoints
    └── results/                    # Stage 5 checkpoints
```

## Summary

✅ **Automatic checkpointing** - No manual tracking needed
✅ **Resume anywhere** - Restart from last successful stage  
✅ **Parallel execution** - Multiple tissues processed simultaneously
✅ **Smart dependencies** - Only rerun what changed
✅ **Resource aware** - Optimal memory/CPU allocation per stage
✅ **Production ready** - Cluster integration with SLURM

Simply run: `sbatch submit_pipeline.sh` and let Snakemake handle the rest!
