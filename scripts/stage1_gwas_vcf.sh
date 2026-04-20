#!/bin/bash
#SBATCH --job-name=coloc_stage1
#SBATCH --output=logs/stage1_%j.out
#SBATCH --error=logs/stage1_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal

################################################################################
# Stage 1: GWAS VCF Conversion
# Converts GWAS summary statistics to VCF format with liftover to hg38
################################################################################

# Use local scratch for temporary files (faster I/O)
if [[ -d "/localscratch/${USER}" ]] || mkdir -p "/localscratch/${USER}" 2>/dev/null; then
    export TMPDIR="/localscratch/${USER}"
else
    export TMPDIR="/tmp/${USER}"
fi
mkdir -p "${TMPDIR}"

# Clean up local scratch on exit
trap 'rm -rf "${TMPDIR}"/*' EXIT

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Configuration
PIPELINE_DIR="/home/itg/oleg.vlasovets/projects/coloc-pipeline"
OUTPUT_DIR="/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"
CONDA_ENV="coloc-pipeline"

# Parse command line arguments
TRAIT="${1:-KNEE}"  # Default to KNEE if not provided
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
    
    if [[ ! -f "${PIPELINE_DIR}/config.yaml" ]]; then
        log_error "Config file not found: ${PIPELINE_DIR}/config.yaml"
        exit 1
    fi
    
    log_info "Prerequisites OK"
}

# Main execution
main() {
    log_info "=========================================="
    log_info "Stage 1: GWAS VCF Conversion"
    log_info "=========================================="
    log_info "Trait: ${TRAIT}"
    log_info "Output: ${OUTPUT_DIR}/gwas_vcf/${TRAIT}_hg38.vcf.bgz"
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
    
    # Build snakemake command
    SNAKEMAKE_CMD="snakemake \
        ${OUTPUT_DIR}/gwas_vcf/${TRAIT}_hg38.vcf.bgz \
        --use-conda \
        --cores ${SLURM_CPUS_PER_TASK:-1} \
        --rerun-incomplete \
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
        log_info "Stage 1 completed successfully!"
        log_info "=========================================="
        log_info "Output files:"
        ls -lh "${OUTPUT_DIR}/gwas_vcf/${TRAIT}_hg38.vcf.bgz"* 2>/dev/null || true
        exit 0
    else
        EXIT_CODE=$?
        log_error "=========================================="
        log_error "Stage 1 failed with exit code ${EXIT_CODE}"
        log_error "=========================================="
        log_error "Check log: ${OUTPUT_DIR}/logs/convert_gwas_${TRAIT}.log"
        exit "${EXIT_CODE}"
    fi
}

# Run main function
main "$@"
