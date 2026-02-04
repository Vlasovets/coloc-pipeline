#!/bin/bash
#SBATCH --job-name=coloc_pipeline
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err

# Change to project directory
cd /home/itg/oleg.vlasovets/projects/coloc-pipeline

# Create logs directory if it doesn't exist
mkdir -p logs results

# Load conda
source /home/itg/oleg.vlasovets/miniconda3/etc/profile.d/conda.sh

# Activate environment
conda activate coloc-pipeline

echo "============================================"
echo "Starting Colocalization Pipeline"
echo "Date: $(date)"
echo "============================================"

# Run complete pipeline with automatic checkpoints
# Snakemake will:
# - Skip completed steps (based on output file timestamps)
# - Resume from the last successful stage
# - Parallelize independent jobs across tissues
# - Use cluster submission for compute-intensive rules

snakemake \
    --use-conda \
    --cores 8 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    --reason \
    --latency-wait 60 \
    --cluster "sbatch --partition=cpu_p --qos=cpu_normal --mem={resources.mem_mb}M --time={resources.time} --cpus-per-task={threads} --output=logs/slurm_%j_{rule}_{wildcards}.out --error=logs/slurm_%j_{rule}_{wildcards}.err" \
    --cluster-status "python workflow/scripts/slurm_status.py" \
    --jobs 20 \
    --max-jobs-per-second 10

EXIT_CODE=$?

echo "============================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "Pipeline completed successfully at $(date)"
else
    echo "Pipeline failed with exit code $EXIT_CODE at $(date)"
    echo "Check logs in logs/ directory for details"
fi
echo "============================================"

exit $EXIT_CODE
