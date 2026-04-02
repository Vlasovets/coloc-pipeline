#!/bin/bash
#SBATCH --job-name=coloc_pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal

################################################################################
# Complete Colocalization Pipeline Orchestrator
# Executes all stages sequentially with dependency checking and error handling
#
# Usage:
#   sbatch scripts/submit_pipeline.sh [TRAIT] [MODE]
#
# Arguments:
#   TRAIT: Trait to analyze (default: KNEE)
#   MODE:  execution mode - sequential or snakemake (default: sequential)
#
# Examples:
#   sbatch scripts/submit_pipeline.sh KNEE sequential
#   sbatch scripts/submit_pipeline.sh TKR snakemake
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
SCRIPT_DIR="${PIPELINE_DIR}/scripts"

# Parse arguments
TRAIT="${1:-KNEE}"
MODE="${2:-sequential}"  # sequential or snakemake
TISSUES="high_grade_cartilage low_grade_cartilage synovium fat_pad"

# Functions
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

log_stage() {
    echo ""
    echo "=========================================="
    echo "$*"
    echo "=========================================="
}

run_stage() {
    local stage_name=$1
    local stage_script=$2
    shift 2
    local stage_args=("$@")
    
    log_stage "Running ${stage_name}"
    
    if bash "${stage_script}" "${stage_args[@]}"; then
        log_info "${stage_name} completed successfully"
        return 0
    else
        local exit_code=$?
        log_error "${stage_name} failed with exit code ${exit_code}"
        return ${exit_code}
    fi
}

run_sequential_mode() {
    log_info "Running in SEQUENTIAL mode with explicit stage execution"
    
    # Stage 1: GWAS VCF Conversion
    if ! run_stage "Stage 1: GWAS VCF Conversion" \
        "${SCRIPT_DIR}/stage1_gwas_vcf.sh" "${TRAIT}" "false"; then
        log_error "Pipeline stopped at Stage 1"
        return 1
    fi
    
    # Stage 2: QTL-GWAS Overlap Detection
    if ! run_stage "Stage 2: QTL-GWAS Overlap Detection" \
        "${SCRIPT_DIR}/stage2_overlaps.sh" "${TRAIT}" "${TISSUES}" "false"; then
        log_error "Pipeline stopped at Stage 2"
        return 1
    fi
    
    # Stage 3: Coloc ABF Analysis
    if ! run_stage "Stage 3: Coloc ABF Analysis" \
        "${SCRIPT_DIR}/stage3_coloc_abf.sh" "${TRAIT}" "${TISSUES}" "false"; then
        log_error "Pipeline stopped at Stage 3"
        return 1
    fi
    
    # Stage 4: Aggregate Results
    if ! run_stage "Stage 4: Aggregate Results" \
        "${SCRIPT_DIR}/stage4_aggregate.sh" "${TRAIT}" "false"; then
        log_error "Pipeline stopped at Stage 4"
        return 1
    fi
    
    return 0
}

run_snakemake_mode() {
    log_info "Running in SNAKEMAKE mode (parallelized)"
    
    cd "${PIPELINE_DIR}" || exit 1
    
    source /home/itg/oleg.vlasovets/miniconda3/etc/profile.d/conda.sh
    conda activate coloc-pipeline
    
    # Run all stages for the trait using Snakemake's DAG
    snakemake \
        "${OUTPUT_DIR}/results/${TRAIT}_coloc_aggregated.txt" \
        --use-conda \
        --cores 8 \
        --keep-going \
        --rerun-incomplete \
        --printshellcmds \
        --latency-wait 60
    
    return $?
}

# Main execution
main() {
    log_info "=============================================="
    log_info "Colocalization Pipeline"
    log_info "=============================================="
    log_info "Start time: $(date)"
    log_info "Trait: ${TRAIT}"
    log_info "Mode: ${MODE}"
    log_info "Tissues: ${TISSUES}"
    log_info "Pipeline directory: ${PIPELINE_DIR}"
    log_info "Output directory: ${OUTPUT_DIR}"
    log_info ""
    
    # Create necessary directories
    mkdir -p "${OUTPUT_DIR}"/{gwas_vcf,overlaps,coloc_abf,coloc_susie,results,logs,tmp}
    mkdir -p "${PIPELINE_DIR}/logs"
    
    # Execute based on mode
    if [[ "${MODE}" == "sequential" ]]; then
        if run_sequential_mode; then
            EXIT_CODE=0
        else
            EXIT_CODE=$?
        fi
    elif [[ "${MODE}" == "snakemake" ]]; then
        if run_snakemake_mode; then
            EXIT_CODE=0
        else
            EXIT_CODE=$?
        fi
    else
        log_error "Invalid mode: ${MODE}. Use 'sequential' or 'snakemake'"
        EXIT_CODE=1
    fi
    
    # Final summary
    log_info ""
    log_info "=============================================="
    if [[ ${EXIT_CODE} -eq 0 ]]; then
        log_info "Pipeline completed successfully!"
        log_info "End time: $(date)"
        log_info "=============================================="
        log_info ""
        log_info "Final outputs:"
        ls -lh "${OUTPUT_DIR}/results/${TRAIT}_coloc"*.txt 2>/dev/null || true
    else
        log_error "Pipeline failed with exit code ${EXIT_CODE}"
        log_error "End time: $(date)"
        log_error "=============================================="
        log_error "Check logs in: ${PIPELINE_DIR}/logs/"
    fi
    
    exit ${EXIT_CODE}
}

# Run main function
main "$@"
