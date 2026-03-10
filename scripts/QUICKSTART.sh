#!/bin/bash
# Quick reference script for running different pipeline stages

PIPELINE_DIR="/home/itg/oleg.vlasovets/projects/coloc-pipeline"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Colocalization Pipeline - Quick Start${NC}"
echo "========================================"
echo ""

# Check if conda environment is activated
if [[ "$CONDA_DEFAULT_ENV" != "coloc-pipeline" ]]; then
    echo "Activating conda environment..."
    source /home/itg/oleg.vlasovets/miniconda3/etc/profile.d/conda.sh
    conda activate coloc-pipeline
fi

echo -e "${GREEN}Available Commands:${NC}"
echo ""
echo "1. Dry-run entire pipeline:"
echo "   snakemake --dry-run all"
echo ""
echo "2. Run Stage 1 (GWAS VCF conversion):"
echo "   snakemake results/gwas_vcf/KNEE_hg38.vcf.bgz --use-conda --cores 4"
echo ""
echo "3. Run Stage 2 (QTL-GWAS overlaps) for KNEE:"
echo "   snakemake results/overlaps/KNEE.high_grade_cartilage.overlaps.rda --use-conda --cores 4"
echo ""
echo "4. Run Stage 2 for all tissues (KNEE):"
echo "   snakemake results/overlaps/KNEE.{high_grade_cartilage,low_grade_cartilage,synovium,fat_pad}.overlaps.rda --use-conda --cores 4"
echo ""
echo "5. Submit complete pipeline as SLURM job:"
echo "   sbatch scripts/submit_pipeline.sh"
echo ""
echo "6. Check pipeline status:"
echo "   bash scripts/check_pipeline_status.sh"
echo ""
echo -e "${GREEN}Pipeline Structure:${NC}"
echo "  workflow/rules/     - Modular Snakemake rules (by stage)"
echo "  workflow/scripts/   - R analysis scripts"
echo "  scripts/            - Submission and utility scripts"
echo "  envs/               - Conda environment definitions"
echo "  results/            - Pipeline outputs"
echo "  logs/               - Job logs"
echo ""
