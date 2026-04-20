# Colocalization Pipeline

A Snakemake pipeline for colocalization analysis between GWAS summary statistics and
tissue-specific eQTLs. Developed for osteoarthritis GWAS × joint-tissue eQTL analysis
at the Zeggini lab / ITG, Helmholtz Munich.

---

## Overview

The pipeline tests whether GWAS signals co-localize with eQTLs using two complementary methods:

- **Coloc ABF** (`coloc.abf`) — Approximate Bayes Factor method for genome-wide screening
- **Coloc SuSiE** (`coloc.susie`) — fine-mapping for ambiguous regions with potential multiple causal variants

A locus is flagged as co-localizing when PP4 ≥ 0.8 (strong) or PP4 > 0.6 AND PP4/PP3 > 2 (moderate).
Ambiguous ABF results (0.25 ≤ PP4 < 0.8) are resolved by SuSiE where credible sets are available.

### Applied to

| Trait | Cases | Controls | Description |
|-------|------:|--------:|-------------|
| KNEE | 172,256 | 1,144,244 | Knee osteoarthritis |
| TKR | 48,161 | 958,463 | Total knee replacement |
| ALLOA | 489,952 | 1,471,094 | All osteoarthritis |

| Tissue | eQTL N |
|--------|-------:|
| high_grade_cartilage | 115 |
| low_grade_cartilage | 113 |
| synovium | 109 |
| fat_pad | 97 |

### Results summary

| Trait | Colocalizations (ABF+SuSiE) |
|-------|--------------------------:|
| KNEE | 168 |
| TKR | 174 |
| ALLOA | 200 |

---

## Pipeline Stages

```
GWAS sumstats
     │
     ▼
Stage 1: GWAS VCF conversion          [~1–2 h]
  liftover hg19→hg38, allele harmonization
     │
     ▼  results/gwas_vcf/{trait}_hg38.vcf.bgz
     │
Stage 2: QTL-GWAS overlap detection   [~2 h per tissue]
  GenomicRanges ±1 Mb windows, VCF extraction
     │
     ▼  results/overlaps/{trait}.{tissue}.*
     │
Stage 3: Coloc ABF                    [~25–35 min per tissue]
  coloc.abf, mclapply × 4 cores
     │
     ▼  results/coloc_abf/{trait}.{tissue}.colocABF_results.txt
     │
Stage 4: Aggregate ABF results        [~1 min]
  combine tissues, filter by threshold
     │
     ├── results/results/{trait}_all_coloc_results.txt
     └── results/results/{trait}_significant_coloc_results.txt
     │
Stage 5: Coloc SuSiE fine-mapping     [~30 min per tissue]
  runsusie + coloc.susie for ambiguous ABF pairs (PP4 0.25–0.8)
     │
     ▼  results/coloc_susie/{trait}.{tissue}.colocSuSiE_results.txt
     │
Stage 7: Combine ABF + SuSiE          [~1 min]
  final decision: strong ABF → ABF; ambiguous → SuSiE; fallback → ABF
     │
     ├── results/results/{trait}_combined_coloc_results.txt
     └── results/results/{trait}_combined_significant.txt
```

---

## Requirements

- GWAS summary statistics (hg19 or hg38, gzipped TSV)
- GWAS independent signals file (CSV, hg38)
- eQTL permutation results (`egenes_df.txt`) per tissue
- eQTL nominal results (parquet files, one per chromosome) per tissue
- Variant annotation file (hg38, biallelic SNPs)
- Reference genotype files for LD calculation (plink bfile format)
- SLURM cluster with conda

---

## Quick Start

### 1. Configure

```bash
cp config.yaml.example config.yaml
# Edit paths, traits, tissues, and sample sizes
```

### 2. Submit full pipeline for a trait

```bash
cd /path/to/coloc-pipeline

# Stage 1 — GWAS VCF conversion
sbatch scripts/stage1_gwas_vcf.sh KNEE

# Stage 2 — QTL-GWAS overlaps (after Stage 1)
sbatch scripts/stage2_overlaps.sh KNEE

# Stage 3 — Coloc ABF, one job per tissue (after Stage 2)
sbatch --mem=64G scripts/stage3_coloc_abf.sh KNEE high_grade_cartilage
sbatch --mem=64G scripts/stage3_coloc_abf.sh KNEE low_grade_cartilage
sbatch --mem=64G scripts/stage3_coloc_abf.sh KNEE synovium
sbatch --mem=64G scripts/stage3_coloc_abf.sh KNEE fat_pad

# Stage 4 — Aggregate ABF results (after all Stage 3)
sbatch scripts/stage4_aggregate.sh KNEE

# Stage 5 — SuSiE fine-mapping, one job per tissue (after Stage 4)
sbatch scripts/stage5_susie.sh KNEE high_grade_cartilage
sbatch scripts/stage5_susie.sh KNEE synovium
sbatch scripts/stage5_susie.sh KNEE low_grade_cartilage
sbatch scripts/stage5_susie.sh KNEE fat_pad

# Stage 7 — Combine ABF + SuSiE (after all Stage 5)
sbatch scripts/stage7_combine.sh KNEE
```

### 3. Check status

```bash
squeue -u $USER
tail -f results/logs/coloc_abf_KNEE.high_grade_cartilage.log
```

---

## Repository Structure

```
coloc-pipeline/
├── config.yaml                  # All paths and parameters
├── Snakefile                    # Main Snakemake entry point
├── workflow/
│   ├── rules/                   # Per-stage Snakemake rules
│   └── scripts/                 # R analysis scripts
│       ├── 1_convert_gwas_to_vcf.R
│       ├── 2_find_qtl_gwas_overlaps.R
│       ├── 3_run_coloc_abf.R
│       ├── 4_run_coloc_susie.R
│       ├── 5_postprocess_coloc.R
│       ├── 7_combine_results.R
│       └── Coloc_helper_functions.R
├── scripts/                     # SLURM submission wrappers
│   ├── stage1_gwas_vcf.sh
│   ├── stage2_overlaps.sh
│   ├── stage3_coloc_abf.sh
│   ├── stage4_aggregate.sh
│   ├── stage5_susie.sh
│   └── stage7_combine.sh
├── tests/                       # Unit and integration tests
├── envs/                        # Conda environment definitions
└── markdown/
    ├── PIPELINE_GUIDE.md        # Detailed operational guide
    └── PIPELINE_STATUS.md       # Production run status
```

---

## Output Structure

```
results/
├── gwas_vcf/          # Stage 1: {trait}_hg38.vcf.bgz
├── overlaps/          # Stage 2: {trait}.{tissue}.overlaps.rda / gwas_data.rda
├── coloc_abf/         # Stage 3: {trait}.{tissue}.colocABF_results.txt
├── coloc_susie/       # Stage 5: {trait}.{tissue}.colocSuSiE_results.txt
├── results/           # Stages 4 & 7: aggregated and combined outputs
└── logs/              # Per-stage logs
```

---

## Key Parameters

| Parameter | Value |
|-----------|-------|
| Window size | ±1 Mb around GWAS signal |
| ABF strong threshold | PP4 ≥ 0.8 |
| ABF moderate threshold | PP4 > 0.6 AND PP4/PP3 > 2 |
| SuSiE candidate threshold | PP4 > 0.25 (from ABF) |
| Prior p12 | 1×10⁻⁵ |
| SuSiE coverage | 0.95 |
| SuSiE max signals (L) | 5 |
| overlap.min | 0 (disabled) |

---

## Tests

```bash
# Activate the Snakemake-managed conda env
RSCRIPT=$(find .snakemake/conda -name "Rscript" -path "*/bin/Rscript" | head -1)

${RSCRIPT} tests/test_stage1.R   # 13/13
${RSCRIPT} tests/test_stage2.R   # 16/16
${RSCRIPT} tests/test_stage3.R   # 5/5
${RSCRIPT} tests/test_stage4.R   # 12/12
```

---

## Documentation

- [PIPELINE_GUIDE.md](markdown/PIPELINE_GUIDE.md) — full operational reference, debugging tips, known issues
- [PIPELINE_STATUS.md](markdown/PIPELINE_STATUS.md) — production run history and results

---

© ITG / Zeggini Lab, Helmholtz Munich, 2026
