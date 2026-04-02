#!/bin/bash
#SBATCH --job-name=coloc_susie
#SBATCH --output=logs/stage5_susie_%j.out
#SBATCH --error=logs/stage5_susie_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

################################################################################
# Stage 5: Coloc SuSiE Fine-mapping
# Runs coloc.susie for ambiguous ABF cases (PP4 > 0.25, not already significant)
# Requires: Stage 3 output (ABF results) to exist
# Usage: sbatch scripts/stage5_susie.sh [TRAIT] [TISSUE]
################################################################################

export TMPDIR=/localscratch/${USER}
mkdir -p $TMPDIR
trap "rm -rf $TMPDIR/*" EXIT
set -euo pipefail

PIPELINE_DIR="/home/itg/oleg.vlasovets/projects/coloc-pipeline"
OUTPUT_DIR="/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"
CONDA_ENV="coloc-pipeline"
TRAIT="${1:-KNEE}"
TISSUE="${2:-high_grade_cartilage}"

log_info()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; }

log_info "=========================================="
log_info "Stage 5: Coloc SuSiE Fine-mapping"
log_info "=========================================="
log_info "Trait:  ${TRAIT}"
log_info "Tissue: ${TISSUE}"

# Verify Stage 3 input exists
ABF="${OUTPUT_DIR}/coloc_abf/${TRAIT}.${TISSUE}.colocABF_results.txt"
GWAS_DATA="${OUTPUT_DIR}/overlaps/${TRAIT}.${TISSUE}.gwas_data.rda"
QTL_SUBSET="${OUTPUT_DIR}/overlaps/${TRAIT}.${TISSUE}.qtl_subset.txt"

for f in "$ABF" "$GWAS_DATA" "$QTL_SUBSET"; do
    if [[ ! -f "$f" ]]; then
        log_error "Required input missing: $f"
        log_error "Run Stage 3 first for ${TRAIT} / ${TISSUE}"
        exit 1
    fi
done

# Count SuSiE candidates
NCAND=$(awk -F'\t' 'NR>1 && $22 > 0.25 && $22 < 0.8' "$ABF" | wc -l)
log_info "Candidate gene-signal pairs for SuSiE: ${NCAND}"

if [[ $NCAND -eq 0 ]]; then
    log_info "No candidates — Stage 5 not needed for ${TRAIT} / ${TISSUE}"
    exit 0
fi

TARGET="${OUTPUT_DIR}/coloc_susie/${TRAIT}.${TISSUE}.colocSuSiE_results.txt"
log_info "Target: ${TARGET}"

cd "${PIPELINE_DIR}"
source /home/itg/oleg.vlasovets/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}"

log_info "Unlocking Snakemake..."
snakemake --unlock 2>/dev/null || true

log_info "Running Stage 5 (SuSiE)..."
if snakemake \
    "${TARGET}" \
    --use-conda \
    --rerun-triggers mtime \
    --cores "${SLURM_CPUS_PER_TASK:-2}" \
    --rerun-incomplete \
    --printshellcmds \
    --latency-wait 60; then
    log_info "=========================================="
    log_info "Stage 5 completed!"
    log_info "=========================================="
    ls -lh "${TARGET}"
    if [[ -f "${TARGET}" ]]; then
        NRESULTS=$(( $(wc -l < "${TARGET}") - 1 ))
        NSIG=$(awk -F'\t' 'NR>1 && $5 > 0.8' "${TARGET}" | wc -l)
        log_info "SuSiE tests:            ${NRESULTS}"
        log_info "Significant (PP4>0.8):  ${NSIG}"
    fi
    exit 0
else
    log_error "Stage 5 failed. Check: ${OUTPUT_DIR}/logs/coloc_susie_${TRAIT}.${TISSUE}.log"
    exit 1
fi
