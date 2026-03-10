#!/bin/bash
#SBATCH --job-name=coloc_stage3
#SBATCH --output=logs/stage3_%j.out
#SBATCH --error=logs/stage3_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal

################################################################################
# Stage 3: Coloc ABF Analysis
# Performs colocalization testing using Approximate Bayes Factor approach
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
TISSUES="${2:-high_grade_cartilage low_grade_cartilage synovium fat_pad}"
DRY_RUN="${3:-false}"

# Functions
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

check_prerequisites() {
    log_info "Checking prerequisites..."
    
    # Check if Stage 2 outputs exist
    local missing_files=0
    for tissue in ${TISSUES}; do
        if [[ ! -f "${OUTPUT_DIR}/overlaps/${TRAIT}.${tissue}.overlaps.rda" ]]; then
            log_error "Stage 2 output not found for ${tissue}"
            missing_files=$((missing_files + 1))
        fi
    done
    
    if [[ ${missing_files} -gt 0 ]]; then
        log_error "Please run Stage 2 first for all tissues"
        exit 1
    fi
    
    log_info "Prerequisites OK"
}

# Main execution
main() {
    log_info "=========================================="
    log_info "Stage 3: Coloc ABF Analysis"
    log_info "=========================================="
    log_info "Trait: ${TRAIT}"
    log_info "Tissues: ${TISSUES}"
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
    
    # Build target list
    TARGETS=""
    for tissue in ${TISSUES}; do
        TARGETS="${TARGETS} ${OUTPUT_DIR}/coloc_abf/${TRAIT}.${tissue}.colocABF_results.txt"
    done
    
    # Build snakemake command
    SNAKEMAKE_CMD="snakemake \
        ${TARGETS} \
        --use-conda \
        --cores ${SLURM_CPUS_PER_TASK:-4} \
        --rerun-incomplete \
        --printshellcmds \
        --keep-going \
        --latency-wait 60"
    
    if [[ "${DRY_RUN}" == "true" ]]; then
        SNAKEMAKE_CMD="${SNAKEMAKE_CMD} --dry-run"
        log_info "Running in dry-run mode"
    fi
    
    # Run snakemake
    log_info "Executing snakemake..."
    log_info "This may take several hours..."
    
    if eval "${SNAKEMAKE_CMD}"; then
        log_info "=========================================="
        log_info "Stage 3 completed successfully!"
        log_info "=========================================="
        log_info "Output files:"
        ls -lh "${OUTPUT_DIR}/coloc_abf/${TRAIT}."*.colocABF_results.txt 2>/dev/null || true
        log_info ""
        log_info "Summary:"
        for tissue in ${TISSUES}; do
            if [[ -f "${OUTPUT_DIR}/coloc_abf/${TRAIT}.${tissue}.colocABF_results.txt" ]]; then
                RESULT_COUNT=$(($(wc -l < "${OUTPUT_DIR}/coloc_abf/${TRAIT}.${tissue}.colocABF_results.txt") - 1))
                log_info "  ${tissue}: ${RESULT_COUNT} colocalization tests"
            fi
        done
        exit 0
    else
        EXIT_CODE=$?
        log_error "=========================================="
        log_error "Stage 3 failed with exit code ${EXIT_CODE}"
        log_error "=========================================="
        log_error "Check logs in: ${OUTPUT_DIR}/logs/coloc_abf_${TRAIT}.*.log"
        exit "${EXIT_CODE}"
    fi
}

# Run main function
main "$@"
