#!/bin/bash
# Run complete colocalization pipeline locally (without cluster submission)
# Use this for testing or when running in an interactive session with sufficient resources

set -e  # Exit on error

# Change to project directory
cd /home/itg/oleg.vlasovets/projects/coloc-pipeline

# Create output directories
mkdir -p logs results

# Load conda
source /home/itg/oleg.vlasovets/miniconda3/etc/profile.d/conda.sh

# Activate environment
conda activate coloc-pipeline

echo "============================================"
echo "Starting Colocalization Pipeline (Local Mode)"
echo "Date: $(date)"
echo "Traits: $(grep 'traits:' config.yaml | cut -d'[' -f2 | cut -d']' -f1)"
echo "============================================"
echo ""

# Run Snakemake locally with checkpoints
# Features:
# - Automatic checkpoints (skips completed stages based on file timestamps)
# - Rerun incomplete jobs
# - Keep going on soft failures
# - Print commands and reasons for execution

snakemake \
    --use-conda \
    --cores 8 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    --reason \
    --latency-wait 10 \
    all

EXIT_CODE=$?

echo ""
echo "============================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Pipeline completed successfully at $(date)"
    echo ""
    echo "Results available in:"
    echo "  - GWAS VCFs: results/gwas_vcf/"
    echo "  - Overlaps: results/overlaps/"
    echo "  - Coloc ABF: results/coloc_abf/"
    echo "  - Aggregated: results/results/"
else
    echo "✗ Pipeline failed with exit code $EXIT_CODE at $(date)"
    echo ""
    echo "To resume from last successful checkpoint, simply rerun this script."
    echo "Check logs in results/logs/ for error details."
fi
echo "============================================"

exit $EXIT_CODE
