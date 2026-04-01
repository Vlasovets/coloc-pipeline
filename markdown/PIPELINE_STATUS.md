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
  - `scripts/stage4_aggregate.sh` - Stage 4 runner
  - `scripts/check_status.sh` - Status checker (tested)
  - `scripts/clean.sh` - Cleanup utility

## ⚠️ Needs Implementation

### Stage 3: Coloc ABF Analysis
- **Status**: ✓ Complete and tested
- **Script**: `workflow/scripts/3_run_coloc_abf.R`
- **Rule**: `workflow/rules/stage3_coloc_abf.smk`
- **Test**: `tests/test_stage3.R` — **PASSES** (5/5 gene-signal pairs, all assertions green)
- **Key fixes applied** (branch `o-feature/stage-3-implementation`):
  1. Load `overlaps`/`gwas_windows`/`qtl_windows` from Stage 2 RDA (not `overlap_df`)
  2. Map `overlaps` columns → `perform_coloc()` format (`cpg.id`, `signal_gwas.*`)
  3. Fix `fread(qtl_subset, header=TRUE)$gene_id` (was `header=FALSE`)
  4. Use `snakemake@params$qtl_nominal_dir` + `cis_qtl.cis_qtl_pairs.chr{N}.parquet`
  5. Remove non-existent args (`tissue`, `gwas_cc_ratio`, `GTEX_APP`) from `perform_coloc()` call
  6. Create intermediate sumstats dirs before calling `perform_coloc()`

### Stage 4: Aggregate Results
- **Status**: Template exists, needs full implementation
- **Template**: `workflow/scripts/5_postprocess_coloc.R` (74 lines)
- **Working Code**: `src/5_postprocess_coloc.R`
- **Action Required**: Complete Snakemake integration

**Implementation Plan**:
```r
# Current template has:
- aggregate_coloc_results() function
- summarize_results() function
- CLI parsing

# Needs completion:
- Integration with snakemake@input
- Proper aggregation across tissues
- Summary statistics generation
- Output to snakemake@output files
```

### Stage 5 (Optional): Coloc SuSiE Fine-mapping
- **Status**: Rule defined, script is template
- **Script**: `workflow/scripts/4_run_coloc_susie.R` (template)
- **Priority**: Low (can be added later for ambiguous cases)

## 🚀 Quick Test Options

### Option 1: Test Modular Pipeline (Stages 1-2 Only)
```bash
# This works right now!
bash scripts/check_status.sh KNEE
sbatch scripts/stage2_overlaps.sh KNEE "high_grade_cartilage"
```

### Option 2: Manual Stage 3 Test (Using src/)
```bash
# Run the working src version manually
cd /home/itg/oleg.vlasovets/projects/coloc-pipeline/src
# Edit paths and run 3_run_coloc_abf.R
```

### Option 3: Complete Stage 3-4 Implementation
This requires:
1. Adapting src/3_run_coloc_abf.R to use Snakemake parameters
2. Testing with one tissue  
3. Adapting src/5_postprocess_coloc.R similarly
4. Full pipeline test

## 📊 Current Pipeline Status

```
Pipeline: KNEE

Stage 1 (GWAS VCF):      ✓ Complete (184M, 11.3M SNPs)
Stage 2 (Overlaps):      ✓ Complete (4/4 tissues, 4,342 total genes)
Stage 3 (Coloc ABF):     ✓ Implemented & tested (KNEE/high_grade_cartilage)
Stage 4 (Aggregation):   ⏳ Needs testing
```

## 📝 Next Steps

### Immediate (Testing What Works):
1. ✓ Test Stage 1-2 modular structure - DONE
2. ✓ Test SLURM scripts - DONE  
3. ✓ Verify status checker - DONE

### Short-term (Complete Pipeline):
1. Implement Stage 3 (3_run_coloc_abf.R)
   - Adapt src version to Snakemake
   - Add proper input/output handling
   - Test with one tissue first

2. Implement Stage 4 (5_postprocess_coloc.R)
   - Adapt aggregation logic
   - Add summary statistics
   - Test with Stage 3 output

3. End-to-end test
   - Run full pipeline for KNEE
   - Validate all outputs
   - Document results

### Long-term (Enhancement):
1. Add Stage 5 (SuSiE fine-mapping)
2. Add Stage 6 (Visualization/plots)
3. Optimize parallel processing
4. Add more comprehensive error handling

## 💡 Recommendation

**For immediate testing**: Use the working Stages 1-2 to validate the modular pipeline structure.

**For production**: Complete Stage 3-4 implementation by adapting the working src code to Snakemake interface. This should take 2-3 hours for someone familiar with R and Snakemake.

**Alternative**: Keep using src scripts directly for Stages 3-4 until Snakemake versions are ready.
