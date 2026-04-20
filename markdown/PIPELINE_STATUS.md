# Colocalization Pipeline - Implementation Status

## ✅ Completed & Tested

### Stage 1: GWAS VCF Conversion
- **Status**: ✓ Complete and tested
- **Script**: `workflow/scripts/1_convert_gwas_to_vcf.R`
- **Test**: `tests/test_stage1.R` — PASSES (13/13)
- **Output**: KNEE: 184M VCF, 11.3M SNPs, hg38

### Stage 2: QTL-GWAS Overlap Detection
- **Status**: ✓ Complete and tested
- **Script**: `workflow/scripts/2_find_qtl_gwas_overlaps.R`
- **Test**: `tests/test_stage2.R` — PASSES (16/16)
- **Output**: Overlaps for all 4 tissues (KNEE):

  | Tissue | Genes with overlapping GWAS signal |
  |--------|------------------------------------|
  | high_grade_cartilage | 1,158 |
  | low_grade_cartilage | 1,630 |
  | synovium | 1,448 |
  | fat_pad | 106 |

### Stage 3: Coloc ABF Analysis
- **Status**: ✓ Complete — all 4 tissues (SLURM jobs 34965770, 34984512–14)
- **Script**: `workflow/scripts/3_run_coloc_abf.R`
- **Rule**: `workflow/rules/stage3_coloc_abf.smk`
- **Test**: `tests/test_stage3.R` — PASSES (5/5 gene-signal pairs)
- **Production runs** (KNEE):

  | Tissue | Tests | Runtime | Job |
  |--------|-------|---------|-----|
  | high_grade_cartilage | 2,287 | ~25 min | 34965770 |
  | low_grade_cartilage | 3,175 | 33 min | 34984512 |
  | synovium | 2,828 | 28 min | 34984513 |
  | fat_pad | 225 | 10 min | 34984514 |
  | **Total** | **8,515** | | |

### Stage 4: Aggregate Results (ABF)
- **Status**: ✓ Complete — all 4 tissues (SLURM jobs 34982216, 34984785)
- **Script**: `workflow/scripts/5_postprocess_coloc.R`
- **Rule**: `workflow/rules/stage5_postprocess.smk`
- **Test**: `tests/test_stage4.R` — PASSES (Part A: 7/7, Part B: 5/5)
- **Production run**: KNEE, all 4 tissues — **211 significant colocalizations**

  | Tissue | Significant (PP4≥0.8) |
  |--------|----------------------|
  | high_grade_cartilage | 54 |
  | low_grade_cartilage | 81 |
  | synovium | 71 |
  | fat_pad | 5 |
  | **Total** | **211** |
  | Unique GWAS signals | 93 |
  | Unique genes | 78 |

- **Output**:
  - `results/KNEE_all_coloc_results.txt` (8,515 rows)
  - `results/KNEE_significant_coloc_results.txt` (211 rows)

### Stage 5: Coloc SuSiE Fine-mapping
- **Status**: ✓ Complete for KNEE (all 4 tissues)
- **Script**: `workflow/scripts/4_run_coloc_susie.R`
- **Rule**: `workflow/rules/stage4_coloc_susie.smk`
- **Key fixes required** (see `CLAUDE.md` for details):
  - Filter QTL to coloc_region before SuSiE (not full cis-window)
  - Rename lbf_variable columns to positional keys (`pos_NNNNN`) before `coloc.susie()`
  - Pass `overlap.min=0` to disable posterior trimming filter
- **Production runs** (KNEE):

  | Tissue | Pairs tested | Results | Colocalise=TRUE (PP.H4>0.8) | Jobs |
  |--------|-------------|---------|----------------------------|------|
  | synovium | 11/50 | 19 | 10 | 35053371 |
  | high_grade_cartilage | 9/72 | 19 | 1 | 35164360 |
  | low_grade_cartilage | 14/68 | 25 | 4 | 35164361 |
  | fat_pad | 0 candidates | — | — | 35165852 |
  | **Total** | | **63** | **15** | |

- **Notable KNEE SuSiE hits**:
  - `ENSG00000176371 / rs11633450` — colocates in all 3 tissues (PP.H4 ≥ 0.83)
  - `ENSG00000257815 / rs1152966` (synovium) — PP.H4 = 0.9999
  - `ENSG00000180155 / rs1114761` (low_grade_cartilage) — PP.H4 up to 0.995
  - `ENSG00000171368 / rs11749927` (synovium) — PP.H4 = 0.78

### Modular Snakemake Pipeline
- **Status**: ✓ Complete
- **Structure**: Main `Snakefile` + 5 stage-specific rule files in `workflow/rules/`

### SLURM Submission Scripts
- **Status**: ✓ Complete and tested
- `scripts/stage1_gwas_vcf.sh`, `stage2_overlaps.sh`, `stage3_coloc_abf.sh`,
  `stage3_only.sh`, `stage4_aggregate.sh`, `stage5_susie.sh`,
  `submit_pipeline.sh`, `check_status.sh`, `clean.sh`

---

## 📊 Current Pipeline Status (2026-04-14)

```
KNEE / all 4 tissues
  Stage 1 (GWAS VCF):    ✓ Complete (184M, 11.3M SNPs)
  Stage 2 (Overlaps):    ✓ Complete (4/4 tissues, 4,342 total genes)
  Stage 3 (Coloc ABF):   ✓ Complete (8,515 tests across 4 tissues)
  Stage 4 (Aggregation): ✓ Complete (211 significant, 93 signals, 78 genes)
  Stage 5 (SuSiE):       ✓ Complete (63 results, 15 Colocalise=TRUE)

TKR / all 4 tissues
  Stage 1–4 (pipeline):  ⏳ Running — SLURM job 35165853
  Stage 5 (SuSiE):       ⏳ Pending — jobs 35165940–35165943 (dependency on 35165853)

ALLOA / all 4 tissues
  Stage 1–4 (pipeline):  ⏳ Running — SLURM job 35165854
  Stage 5 (SuSiE):       ⏳ Pending — jobs 35165944–35165947 (dependency on 35165854)
```

---

## 📝 Next Steps

1. **Post-processing**: After TKR/ALLOA Stage 5 completes, create combined final output
   merging ABF and SuSiE results across all 3 traits and 4 tissues
2. **Merge branch**: `o-feature/stage-3-implementation` → `main` via PR
3. **Visualization**: Add locus zoom plots (`workflow/scripts/6_generate_plots.R`)
