#!/bin/bash
#SBATCH --job-name=coloc_stage4
#SBATCH --output=logs/stage4_%j.out
#SBATCH --error=logs/stage4_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal

################################################################################
# Stage 4: Aggregate Results
# Combines colocalization results across tissues and generates summary
################################################################################

# Use local scratch for temporary files (faster I/O)
export TMPDIR=/localscratch/${USER}
mkdir -p $TMPDIR

# Clean up local scratch on exit
trap "rm -rf $TMPDIR/*" EXIT

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Configuration
PIPELINE_DIR="/home/itg/oleg.vlasovets/projects/coloc-pipeline"
OUTPUT_DIR="/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"
CONDA_ENV="coloc-pipeline"

# Parse command line arguments
TRAIT="${1:-KNEE}"
DRY_RUN="${2:-false}"

# Functions
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

check_prerequisites() {
    log_info "Checking prerequisites..."
    
    # Check if Stage 3 outputs exist
    local tissue_count=0
    for tissue in high_grade_cartilage low_grade_cartilage synovium fat_pad; do
        if [[ -f "${OUTPUT_DIR}/coloc_abf/${TRAIT}.${tissue}.colocABF_results.txt" ]]; then
            tissue_count=$((tissue_count + 1))
        fi
    done
    
    if [[ ${tissue_count} -eq 0 ]]; then
        log_error "No Stage 3 outputs found"
        log_error "Please run Stage 3 first"
        exit 1
    fi
    
    log_info "Found ${tissue_count} tissue results"
}

# Main execution
main() {
    log_info "=========================================="
    log_info "Stage 4: Aggregate Results"
    log_info "=========================================="
    log_info "Trait: ${TRAIT}"
    log_info "Dry-run: ${DRY_RUN}"
    log_info ""
    
    # Change to pipeline directory
    cd "${PIPELINE_DIR}" || exit 1
    
    # Check prerequisites
    check_prerequisites
    
    # Load conda
    log_info "Activating conda environment: ${CONDA_ENV}"
    source /home/itg/oleg.vlasovets/miniconda3/etc/profile.d/conda.sh
    conda activate "${CONDA_ENV}"
    
    # Unlock directory in case of previous failed runs
    log_info "Unlocking Snakemake directory..."
    snakemake --unlock || log_info "Directory already unlocked or no lock present"
    
    # Build snakemake command
    # Target: _all_coloc_results.txt (aggregate_results rule in stage5_postprocess.smk)
    SNAKEMAKE_CMD="snakemake \
        ${OUTPUT_DIR}/results/${TRAIT}_all_coloc_results.txt \
        --use-conda \
        --cores ${SLURM_CPUS_PER_TASK:-1} \
        --rerun-incomplete \
        --rerun-triggers mtime \
        --printshellcmds \
        --latency-wait 60"
    
    if [[ "${DRY_RUN}" == "true" ]]; then
        SNAKEMAKE_CMD="${SNAKEMAKE_CMD} --dry-run"
        log_info "Running in dry-run mode"
    fi
    
    # Run snakemake
    log_info "Executing snakemake..."
    if eval "${SNAKEMAKE_CMD}"; then
        log_info "=========================================="
        log_info "Stage 4 completed successfully!"
        log_info "=========================================="
        log_info "Output files:"
        ls -lh "${OUTPUT_DIR}/results/${TRAIT}_all_coloc_results.txt" \
                "${OUTPUT_DIR}/results/${TRAIT}_significant_coloc_results.txt" 2>/dev/null || true

        log_info ""
        log_info "Significant colocalizations:"
        if [[ -f "${OUTPUT_DIR}/results/${TRAIT}_significant_coloc_results.txt" ]]; then
            NSIG=$(( $(wc -l < "${OUTPUT_DIR}/results/${TRAIT}_significant_coloc_results.txt") - 1 ))
            log_info "  ${NSIG} significant colocalization(s)"
        fi
        exit 0
    else
        EXIT_CODE=$?
        log_error "=========================================="
        log_error "Stage 4 failed with exit code ${EXIT_CODE}"
        log_error "=========================================="
        log_error "Check log: ${OUTPUT_DIR}/logs/aggregate_${TRAIT}.log"
        exit "${EXIT_CODE}"
    fi
}

# Run main function
main "$@"
