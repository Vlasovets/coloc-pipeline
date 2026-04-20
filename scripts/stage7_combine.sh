#!/bin/bash
#SBATCH --job-name=coloc_combine
#SBATCH --output=logs/stage7_combine_%j.out
#SBATCH --error=logs/stage7_combine_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000

################################################################################
# Stage 7: Combine ABF + SuSiE results into a single final output
# Usage: sbatch scripts/stage7_combine.sh [TRAIT]
# Example: sbatch scripts/stage7_combine.sh KNEE
################################################################################

set -euo pipefail

PIPELINE_DIR="/home/itg/oleg.vlasovets/projects/coloc-pipeline"
OUTPUT_DIR="/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"
CONDA_ENV="coloc-pipeline"
TRAIT="${1:-KNEE}"

log_info() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"; }

log_info "Stage 7: Combine ABF + SuSiE — ${TRAIT}"

cd "${PIPELINE_DIR}"

# Use the Rscript from the Snakemake-managed r_coloc conda env (has data.table, coloc)
RSCRIPT=$(find .snakemake/conda -name "Rscript" -path "*/bin/Rscript" | head -1)
if [[ -z "${RSCRIPT}" ]]; then
    log_info "ERROR: Snakemake conda Rscript not found — run a Snakemake stage first"
    exit 1
fi
log_info "Using Rscript: ${RSCRIPT}"

"${RSCRIPT}" workflow/scripts/7_combine_results.R "${TRAIT}" "${OUTPUT_DIR}"

TARGET="${OUTPUT_DIR}/results/${TRAIT}_combined_significant.txt"
if [[ -f "${TARGET}" ]]; then
    NSIG=$(( $(wc -l < "${TARGET}") - 1 ))
    log_info "Done — ${NSIG} significant colocalizations written to ${TARGET}"
fi
