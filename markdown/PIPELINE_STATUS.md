# Colocalization Pipeline - Implementation Status

## ✅ Completed & Tested

### Stage 1: GWAS VCF Conversion
- **Status**: ✓ Complete and tested
- **Script**: `workflow/scripts/1_convert_gwas_to_vcf.R`
- **Output**: VCF files with liftover to hg38
- **Test Result**: Success (KNEE: 184M VCF, 11.3M SNPs)

### Stage 2: QTL-GWAS Overlap Detection  
- **Status**: ✓ Complete and tested
- **Script**: `workflow/scripts/2_find_qtl_gwas_overlaps.R`
- **Output**: Overlap files for all 4 tissues
- **Test Result**: Success
  - high_grade_cartilage: 1,158 genes
  - low_grade_cartilage: 1,630 genes
  - synovium: 1,448 genes
  - fat_pad: 106 genes

### Stage 3: Coloc ABF Analysis
- **Status**: ✓ Complete and tested (SLURM job 34965770)
- **Script**: `workflow/scripts/3_run_coloc_abf.R`
- **Rule**: `workflow/rules/stage3_coloc_abf.smk`
- **Test**: `tests/test_stage3.R` — PASSES (5/5 gene-signal pairs, all assertions green)
- **Production run**: KNEE / high_grade_cartilage — 2,287 colocalization tests
- **Key fixes applied** (branch `o-feature/stage-3-implementation`):
  1. Load `overlaps`/`gwas_windows`/`qtl_windows` from Stage 2 RDA (not `overlap_df`)
  2. Map `overlaps` columns → `perform_coloc()` format (`cpg.id`, `signal_gwas.*`)
  3. Fix `fread(qtl_subset, header=TRUE)$gene_id` (was `header=FALSE`)
  4. Use `snakemake@params$qtl_nominal_dir` + `cis_qtl.cis_qtl_pairs.chr{N}.parquet`
  5. Remove non-existent args (`tissue`, `gwas_cc_ratio`, `GTEX_APP`) from `perform_coloc()` call
  6. Create intermediate sumstats dirs before calling `perform_coloc()`

### Stage 4: Aggregate Results
- **Status**: ✓ Complete and tested (SLURM job 34982216)
- **Script**: `workflow/scripts/5_postprocess_coloc.R`
- **Rule**: `workflow/rules/stage5_postprocess.smk`
- **Test**: `tests/test_stage4.R` — PASSES (Part A: 7/7 mock tests, Part B: 5/5 real data tests)
- **Production run**: KNEE — 2,287 total tests, **54 significant colocalizations** (40 signals, 32 genes)
- **Output**:
  - `results/KNEE_all_coloc_results.txt` (2,287 rows)
  - `results/KNEE_significant_coloc_results.txt` (54 rows)
- **Key fixes applied**:
  1. Fixed tissue extraction regex: `sub("^[^.]+\\.(.+)\\.colocABF_results\\.txt$", "\\1", basename(f))`
  2. Fixed GWAS_ID extraction regex: `sub("^([^.]+)\\..*", "\\1", basename(f))`
  3. Dynamic Snakemake input function (`_available_coloc_results`) to support partial tissue results

### Modular Snakemake Pipeline
- **Status**: ✓ Complete
- **Structure**: 
  - Main `Snakefile` (58 lines, imports modular rules)
  - 5 stage-specific rule files in `workflow/rules/`
  - Clean separation of concerns

### SLURM Submission Scripts
- **Status**: ✓ Complete and tested
- **Scripts**:
  - `scripts/submit_pipeline.sh` - Main orchestrator
  - `scripts/stage1_gwas_vcf.sh` - Stage 1 runner
  - `scripts/stage2_overlaps.sh` - Stage 2 runner  
  - `scripts/stage3_coloc_abf.sh` - Stage 3 runner
  - `scripts/stage3_only.sh` - Stage 3 runner (skips Stage 2 re-run via `--rerun-triggers mtime`)
  - `scripts/stage4_aggregate.sh` - Stage 4 runner
  - `scripts/check_status.sh` - Status checker (tested)
  - `scripts/clean.sh` - Cleanup utility

## 📊 Current Pipeline Status

```
Pipeline: KNEE / high_grade_cartilage

Stage 1 (GWAS VCF):      ✓ Complete (184M, 11.3M SNPs)
Stage 2 (Overlaps):      ✓ Complete (4/4 tissues, 4,342 total genes)
Stage 3 (Coloc ABF):     ✓ Complete & tested (2,287 tests, job 34965770)
Stage 4 (Aggregation):   ✓ Complete & tested (54 significant, job 34982216)
```

## ⚠️ Remaining / Optional

### Stage 5 (Optional): Coloc SuSiE Fine-mapping
- **Status**: Rule defined, script is template
- **Script**: `workflow/scripts/4_run_coloc_susie.R` (template)
- **Priority**: Low (can be added later for ambiguous cases)

## 📝 Next Steps

### Short-term (Extend to More Tissues):
1. Run Stage 3 for remaining tissues: `low_grade_cartilage`, `synovium`, `fat_pad`
   ```bash
   sbatch scripts/stage3_only.sh KNEE low_grade_cartilage
   sbatch scripts/stage3_only.sh KNEE synovium
   sbatch scripts/stage3_only.sh KNEE fat_pad
   ```
2. Re-run Stage 4 to aggregate all tissues
   ```bash
   sbatch scripts/stage4_aggregate.sh KNEE
   ```

### Long-term (Enhancement):
1. Add Stage 5 (SuSiE fine-mapping)
2. Add visualization/summary plots
3. Extend to additional GWAS traits
