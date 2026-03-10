#!/bin/bash

################################################################################
# Pipeline Status Checker
# Checks the completion status of all pipeline stages for a given trait
#
# Usage: bash scripts/check_status.sh [TRAIT]
################################################################################

set -euo pipefail

# Configuration
OUTPUT_DIR="/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"
TRAIT="${1:-KNEE}"
TISSUES=(high_grade_cartilage low_grade_cartilage synovium fat_pad)

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

check_file() {
    if [[ -f "$1" ]]; then
        echo -e "${GREEN}✓${NC}"
        return 0
    else
        echo -e "${RED}✗${NC}"
        return 1
    fi
}

get_file_info() {
    if [[ -f "$1" ]]; then
        ls -lh "$1" | awk '{print $5, $6, $7, $8}'
    else
        echo "Not found"
    fi
}

echo -e "${BLUE}================================================${NC}"
echo -e "${BLUE}Pipeline Status Check: ${TRAIT}${NC}"
echo -e "${BLUE}================================================${NC}"
echo ""

# Stage 1: GWAS VCF
echo -e "${YELLOW}Stage 1: GWAS VCF Conversion${NC}"
VCF_FILE="${OUTPUT_DIR}/gwas_vcf/${TRAIT}_hg38.vcf.bgz"
echo -n "  VCF file: "
if check_file "${VCF_FILE}"; then
    echo "    $(get_file_info ${VCF_FILE})"
fi
echo -n "  Index:    "
check_file "${VCF_FILE}.tbi" && echo "    $(get_file_info ${VCF_FILE}.tbi)"
echo ""

# Stage 2: Overlaps
echo -e "${YELLOW}Stage 2: QTL-GWAS Overlaps${NC}"
stage2_complete=0
for tissue in "${TISSUES[@]}"; do
    OVERLAP_FILE="${OUTPUT_DIR}/overlaps/${TRAIT}.${tissue}.overlaps.rda"
    QTL_FILE="${OUTPUT_DIR}/overlaps/${TRAIT}.${tissue}.qtl_subset.txt"
    
    echo -n "  ${tissue}: "
    if check_file "${OVERLAP_FILE}"; then
        if [[ -f "${QTL_FILE}" ]]; then
            GENE_COUNT=$(($(wc -l < "${QTL_FILE}") - 1))
            echo "    ${GENE_COUNT} genes with overlaps"
            stage2_complete=$((stage2_complete + 1))
        fi
    fi
done
echo ""

# Stage 3: Coloc ABF
echo -e "${YELLOW}Stage 3: Coloc ABF Analysis${NC}"
stage3_complete=0
for tissue in "${TISSUES[@]}"; do
    ABF_FILE="${OUTPUT_DIR}/coloc_abf/${TRAIT}.${tissue}.colocABF_results.txt"
    echo -n "  ${tissue}: "
    if check_file "${ABF_FILE}"; then
        RESULT_COUNT=$(($(wc -l < "${ABF_FILE}") - 1))
        echo "    ${RESULT_COUNT} colocalization tests"
        stage3_complete=$((stage3_complete + 1))
    fi
done
echo ""

# Stage 4: Aggregated results
echo -e "${YELLOW}Stage 4: Aggregated Results${NC}"
AGG_FILE="${OUTPUT_DIR}/results/${TRAIT}_coloc_aggregated.txt"
SUM_FILE="${OUTPUT_DIR}/results/${TRAIT}_coloc_summary.txt"

echo -n "  Aggregated: "
check_file "${AGG_FILE}" && echo "    $(get_file_info ${AGG_FILE})"
echo -n "  Summary:    "
check_file "${SUM_FILE}" && echo "    $(get_file_info ${SUM_FILE})"
echo ""

# Summary
echo -e "${BLUE}================================================${NC}"
echo -e "${BLUE}Summary${NC}"
echo -e "${BLUE}================================================${NC}"
echo -e "Stage 1 (GWAS VCF):       $(check_file ${VCF_FILE} && echo -e ${GREEN}Complete${NC} || echo -e ${RED}Incomplete${NC})"
echo -e "Stage 2 (Overlaps):       ${stage2_complete}/4 tissues"
echo -e "Stage 3 (Coloc ABF):      ${stage3_complete}/4 tissues"
echo -e "Stage 4 (Aggregated):     $(check_file ${AGG_FILE} && echo -e ${GREEN}Complete${NC} || echo -e ${RED}Incomplete${NC})"
echo ""

# Overall status
if [[ -f "${AGG_FILE}" ]]; then
    echo -e "${GREEN}✓ Pipeline COMPLETE for ${TRAIT}${NC}"
    exit 0
elif [[ ${stage3_complete} -gt 0 ]]; then
    echo -e "${YELLOW}⚠ Pipeline IN PROGRESS for ${TRAIT}${NC}"
    exit 0
else
    echo -e "${RED}✗ Pipeline NOT STARTED or FAILED for ${TRAIT}${NC}"
    exit 1
fi
