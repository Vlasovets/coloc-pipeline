#!/bin/bash
# Check pipeline progress and show completed/pending stages

cd /home/itg/oleg.vlasovets/projects/coloc-pipeline

source /home/itg/oleg.vlasovets/miniconda3/etc/profile.d/conda.sh
conda activate coloc-pipeline

echo "============================================"
echo "Colocalization Pipeline Status"
echo "============================================"
echo ""

# Load traits and tissues from config
TRAITS=$(grep 'traits:' config.yaml | sed 's/.*\[//' | sed 's/\].*//' | tr ',' ' ')
TISSUES="high_grade_cartilage low_grade_cartilage synovium fat_pad"

echo "Configuration:"
echo "  Traits: $TRAITS"
echo "  Tissues: $(echo $TISSUES | wc -w) tissues"
echo ""

# Check Stage 1: GWAS VCF conversion
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Stage 1: GWAS VCF Conversion"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
for trait in $TRAITS; do
    FILE="results/gwas_vcf/${trait}_hg38.vcf.gz"
    if [ -f "$FILE" ]; then
        SIZE=$(du -h "$FILE" | cut -f1)
        echo "  ✓ $trait ($SIZE)"
    else
        echo "  ✗ $trait (pending)"
    fi
done
echo ""

# Check Stage 2: QTL-GWAS Overlaps
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Stage 2: QTL-GWAS Overlap Detection"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
TOTAL_OVERLAPS=0
COMPLETE_OVERLAPS=0
for trait in $TRAITS; do
    echo "  $trait:"
    for tissue in $TISSUES; do
        FILE="results/overlaps/${trait}.${tissue}.overlaps.rda"
        TOTAL_OVERLAPS=$((TOTAL_OVERLAPS + 1))
        if [ -f "$FILE" ]; then
            COMPLETE_OVERLAPS=$((COMPLETE_OVERLAPS + 1))
            GENES=$(ls results/overlaps/${trait}.${tissue}.qtl_subset.txt 2>/dev/null | xargs wc -l 2>/dev/null | awk '{print $1-1}' || echo "?")
            echo "    ✓ $(echo $tissue | sed 's/_/ /g') (${GENES} genes)"
        else
            echo "    ✗ $(echo $tissue | sed 's/_/ /g') (pending)"
        fi
    done
done
echo "  Progress: $COMPLETE_OVERLAPS/$TOTAL_OVERLAPS completed"
echo ""

# Check Stage 3: Coloc ABF
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Stage 3: Colocalization Analysis (ABF)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
TOTAL_COLOC=0
COMPLETE_COLOC=0
for trait in $TRAITS; do
    echo "  $trait:"
    for tissue in $TISSUES; do
        FILE="results/coloc_abf/${trait}.${tissue}.colocABF_results.txt"
        TOTAL_COLOC=$((TOTAL_COLOC + 1))
        if [ -f "$FILE" ]; then
            COMPLETE_COLOC=$((COMPLETE_COLOC + 1))
            TESTS=$(wc -l < "$FILE" 2>/dev/null | awk '{print $1-1}' || echo "?")
            echo "    ✓ $(echo $tissue | sed 's/_/ /g') (${TESTS} tests)"
        else
            echo "    ✗ $(echo $tissue | sed 's/_/ /g') (pending)"
        fi
    done
done
echo "  Progress: $COMPLETE_COLOC/$TOTAL_COLOC completed"
echo ""

# Check Stage 5: Aggregated results
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Stage 5: Result Aggregation"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
for trait in $TRAITS; do
    FILE="results/results/${trait}_coloc_aggregated.txt"
    if [ -f "$FILE" ]; then
        RESULTS=$(wc -l < "$FILE" 2>/dev/null | awk '{print $1-1}')
        echo "  ✓ $trait ($RESULTS results)"
    else
        echo "  ✗ $trait (pending)"
    fi
done
echo ""

# Check recent activity
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Recent Activity (last 5 log files):"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if [ -d "results/logs" ]; then
    ls -lt results/logs/*.log 2>/dev/null | head -5 | awk '{print "  " $9 " (" $6 " " $7 " " $8 ")"}'
else
    echo "  No logs found"
fi
echo ""

# Show incomplete jobs
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "To resume pipeline:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Cluster mode:  sbatch submit_pipeline.sh"
echo "  Local mode:    bash run_pipeline_local.sh"
echo "  Dry run:       snakemake --dry-run -n"
echo ""
echo "Snakemake will automatically:"
echo "  • Skip completed stages (checkpoints)"
echo "  • Resume from last successful step"
echo "  • Rerun incomplete/failed jobs only"
echo "============================================"
