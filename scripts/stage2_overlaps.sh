#!/bin/bash
#SBATCH --job-name=coloc_stage2
#SBATCH --output=logs/stage2_%j.out
#SBATCH --error=logs/stage2_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal

################################################################################
# Stage 2: QTL-GWAS Overlap Detection
# Identifies genomic regions where GWAS signals overlap with eQTLs
################################################################################
# Use local scratch for temporary files (faster I/O)
export TMPDIR=/localscratch/${USER}
mkdir -p $TMPDIR

# Set R to use local scratch for cache and temporary files
export R_USER_CACHE_DIR=${TMPDIR}/R_cache
export R_USER_DATA_DIR=${TMPDIR}/R_data
mkdir -p ${R_USER_CACHE_DIR} ${R_USER_DATA_DIR}

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
TEST_MODE="${4:-false}"  # Set to "true" to process only first 10 signals for debugging

# Export test mode for R script
export COLOC_TEST_MODE="${TEST_MODE}"

# Functions
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

check_prerequisites() {
    log_info "Checking prerequisites..."
    
    # Check if Stage 1 output exists
    if [[ ! -f "${OUTPUT_DIR}/gwas_vcf/${TRAIT}_hg38.vcf.bgz" ]]; then
        log_error "Stage 1 output not found: ${OUTPUT_DIR}/gwas_vcf/${TRAIT}_hg38.vcf.bgz"
        log_error "Please run Stage 1 first"
        exit 1
    fi
    
    log_info "Prerequisites OK"
}

# Main execution
main() {
    log_info "=========================================="
    log_info "Stage 2: QTL-GWAS Overlap Detection"
    log_info "=========================================="
    log_info "Trait: ${TRAIT}"
    log_info "Tissues: ${TISSUES}"
    log_info "Dry-run: ${DRY_RUN}"
    log_info "Test mode: ${TEST_MODE}"
    if [[ "${TEST_MODE}" == "true" ]]; then
        log_info "*** TESTING: Will process only first 10 GWAS signals ***"
    fi
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
    
    # Build target list (include gwas_data.rda to force rerun if missing)
    TARGETS=""
    for tissue in ${TISSUES}; do
        TARGETS="${TARGETS} ${OUTPUT_DIR}/overlaps/${TRAIT}.${tissue}.overlaps.rda"
        TARGETS="${TARGETS} ${OUTPUT_DIR}/overlaps/${TRAIT}.${tissue}.gwas_data.rda"
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
    if eval "${SNAKEMAKE_CMD}"; then
        log_info "=========================================="
        log_info "Stage 2 completed successfully!"
        log_info "=========================================="
        log_info "Output files:"
        ls -lh "${OUTPUT_DIR}/overlaps/${TRAIT}."*.overlaps.rda 2>/dev/null || true
        log_info ""
        log_info "Summary:"
        for tissue in ${TISSUES}; do
            if [[ -f "${OUTPUT_DIR}/overlaps/${TRAIT}.${tissue}.qtl_subset.txt" ]]; then
                GENE_COUNT=$(($(wc -l < "${OUTPUT_DIR}/overlaps/${TRAIT}.${tissue}.qtl_subset.txt") - 1))
                log_info "  ${tissue}: ${GENE_COUNT} genes with overlaps"
            fi
        done
        exit 0
    else
        EXIT_CODE=$?
        log_error "=========================================="
        log_error "Stage 2 failed with exit code ${EXIT_CODE}"
        log_error "=========================================="
        log_error "Check logs in: ${OUTPUT_DIR}/logs/overlaps_${TRAIT}.*.log"
        exit "${EXIT_CODE}"
    fi
}

# Run main function
main "$@"
