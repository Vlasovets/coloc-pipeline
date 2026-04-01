#!/bin/bash
#SBATCH --job-name=coloc_stage3_only
#SBATCH --output=logs/stage3_only_%j.out
#SBATCH --error=logs/stage3_only_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal

################################################################################
# Stage 3: Coloc ABF Analysis — skips Stage 2 re-run via --rerun-triggers mtime
# Use this script after Stage 2 outputs already exist but Stage 3 hasn't run.
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

log_info() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; }

TARGET="${OUTPUT_DIR}/coloc_abf/${TRAIT}.${TISSUE}.colocABF_results.txt"

log_info "Stage 3 (coloc ABF only — no Stage 2 re-run)"
log_info "Trait: ${TRAIT}  Tissue: ${TISSUE}"
log_info "Target: ${TARGET}"

# Verify Stage 2 inputs exist
for f in \
    "${OUTPUT_DIR}/overlaps/${TRAIT}.${TISSUE}.overlaps.rda" \
    "${OUTPUT_DIR}/overlaps/${TRAIT}.${TISSUE}.qtl_subset.txt" \
    "${OUTPUT_DIR}/overlaps/${TRAIT}.${TISSUE}.gwas_data.rda" \
    "${OUTPUT_DIR}/gwas_vcf/${TRAIT}_hg38.vcf.bgz"; do
    if [[ ! -f "$f" ]]; then
        log_error "Required input missing: $f"
        log_error "Run Stage 2 first (scripts/stage2_overlaps.sh ${TRAIT} ${TISSUE})"
        exit 1
    fi
done
log_info "All Stage 2 inputs present"

cd "${PIPELINE_DIR}"
source /home/itg/oleg.vlasovets/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}"

log_info "Unlocking Snakemake directory..."
snakemake --unlock 2>/dev/null || true

log_info "Running Stage 3 (mtime triggers only — skips Stage 2 re-run)..."
if snakemake \
    "${TARGET}" \
    --use-conda \
    --rerun-triggers mtime \
    --cores "${SLURM_CPUS_PER_TASK:-4}" \
    --rerun-incomplete \
    --printshellcmds \
    --latency-wait 60; then
    log_info "Stage 3 completed!"
    ls -lh "${TARGET}"
    NRESULTS=$(( $(wc -l < "${TARGET}") - 1 ))
    log_info "Colocalization tests: ${NRESULTS}"
    exit 0
else
    log_error "Stage 3 failed. Check: ${OUTPUT_DIR}/logs/coloc_abf_${TRAIT}.${TISSUE}.log"
    exit 1
fi
