#!/bin/bash

################################################################################
# Clean Pipeline Outputs
# Removes intermediate files and optionally all results for a fresh start
#
# Usage: 
#   bash scripts/clean.sh [TRAIT] [MODE]
#
# Arguments:
#   TRAIT: Trait to clean (default: KNEE, use 'all' for all traits)
#   MODE:  clean or purge (default: clean)
#          - clean: removes only intermediate files (.snakemake, tmp)
#          - purge: removes all outputs including final results
################################################################################

set -euo pipefail

# Configuration
PIPELINE_DIR="/home/itg/oleg.vlasovets/projects/coloc-pipeline"
OUTPUT_DIR="/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"

TRAIT="${1:-KNEE}"
MODE="${2:-clean}"

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

log_info() {
    echo -e "${GREEN}[INFO]${NC} $*"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $*"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $*"
}

clean_intermediate() {
    log_info "Cleaning intermediate files..."
    
    # Remove Snakemake metadata
    if [[ -d "${PIPELINE_DIR}/.snakemake" ]]; then
        log_info "Removing .snakemake directory"
        rm -rf "${PIPELINE_DIR}/.snakemake"
    fi
    
    # Remove temp files
    if [[ -d "${OUTPUT_DIR}/tmp" ]]; then
        log_info "Removing temporary files"
        rm -rf "${OUTPUT_DIR}/tmp"/*
    fi
    
    # Remove old job logs (keep last 10)
    if [[ -d "${PIPELINE_DIR}/logs" ]]; then
        log_info "Cleaning old job logs"
        cd "${PIPELINE_DIR}/logs"
        ls -t | tail -n +11 | xargs -r rm --
    fi
    
    log_info "Intermediate files cleaned"
}

purge_outputs() {
    local trait=$1
    
    log_warn "PURGING all outputs for trait: ${trait}"
    read -p "Are you sure? This will delete all results! (yes/no): " -r
    
    if [[ $REPLY != "yes" ]]; then
        log_info "Purge cancelled"
        return 0
    fi
    
    if [[ "${trait}" == "all" ]]; then
        log_warn "Purging ALL traits..."
        rm -rf "${OUTPUT_DIR}"/{gwas_vcf,overlaps,coloc_abf,coloc_susie,results}/*
    else
        log_warn "Purging ${trait}..."
        rm -f "${OUTPUT_DIR}/gwas_vcf/${trait}"*
        rm -f "${OUTPUT_DIR}/overlaps/${trait}."*
        rm -f "${OUTPUT_DIR}/coloc_abf/${trait}."*
        rm -f "${OUTPUT_DIR}/coloc_susie/${trait}."*
        rm -f "${OUTPUT_DIR}/results/${trait}_"*
        rm -f "${OUTPUT_DIR}/logs/*${trait}*"
    fi
    
    log_info "Purge complete"
}

# Main
echo "================================================"
echo "Pipeline Cleanup Utility"
echo "================================================"
echo "Trait: ${TRAIT}"
echo "Mode:  ${MODE}"
echo ""

case "${MODE}" in
    clean)
        clean_intermediate
        ;;
    purge)
        clean_intermediate
        purge_outputs "${TRAIT}"
        ;;
    *)
        log_error "Invalid mode: ${MODE}"
        log_error "Use 'clean' or 'purge'"
        exit 1
        ;;
esac

echo ""
log_info "Cleanup complete!"
