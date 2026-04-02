# Coloc Pipeline — Comprehensive Guide

A self-contained reference for running the colocalization pipeline on the HPC cluster.
A new colleague should be able to run the full pipeline independently using this document.

---

## 1. Overview

### What the pipeline does
The coloc pipeline tests whether GWAS signals for osteoarthritis traits (KNEE, TKR, ALLOA)
co-localize with tissue-specific eQTLs in four OA-relevant tissues. Co-localization is
assessed using the Approximate Bayes Factor (ABF) method from the `coloc` R package.
A shared causal variant is flagged when PP4 ≥ 0.8 (strong evidence) or
PP4 > 0.6 AND PP4/PP3 > 2 (moderate evidence).

### Inputs required
- GWAS summary statistics (one file per trait, gzipped TSV)
- eQTL permutation results (`egenes_df.txt`) per tissue
- eQTL nominal results (parquet files per chromosome) per tissue
- Variant annotation file (hg38, biallelic SNPs)
- GWAS independent signals file (CSV, hg38 coordinates)

### Outputs produced
- `results/gwas_vcf/{trait}_hg38.vcf.bgz` — GWAS in VCF format (hg38)
- `results/overlaps/{trait}.{tissue}.*` — QTL-GWAS overlap files
- `results/coloc_abf/{trait}.{tissue}.colocABF_results.txt` — per-tissue coloc results
- `results/results/{trait}_all_coloc_results.txt` — all tests aggregated
- `results/results/{trait}_significant_coloc_results.txt` — significant hits only

### Pipeline architecture

```
GWAS sumstats ─────────────────────────────────────────────┐
                                                            ▼
               Stage 1: GWAS VCF conversion (liftover hg19→hg38)
                              │
                              ▼  results/gwas_vcf/KNEE_hg38.vcf.bgz
                              │
eQTL permutations ────────────┤
GWAS signals (hg38) ──────────┤
                              ▼
               Stage 2: QTL-GWAS overlap detection  [per tissue]
                              │
                              ▼  results/overlaps/KNEE.{tissue}.*.rda/txt
                              │
eQTL nominal (parquet) ───────┤
Variant annotation ───────────┤
                              ▼
               Stage 3: Coloc ABF analysis  [per tissue, ~25 min]
                              │
                              ▼  results/coloc_abf/KNEE.{tissue}.colocABF_results.txt
                              │
                              ▼
               Stage 4: Aggregate + filter significant results
                              │
                              ├─▶ results/results/KNEE_all_coloc_results.txt
                              └─▶ results/results/KNEE_significant_coloc_results.txt
```

---

## 2. Repository Structure

```
coloc-pipeline/
├── Snakefile                    # Main pipeline entry point (imports modular rules)
├── config.yaml                  # All paths and parameters
├── environment.yml              # Conda env for Snakemake itself (coloc-pipeline)
├── envs/
│   └── r_coloc.yml              # R conda env (coloc, arrow, data.table, etc.)
├── workflow/
│   ├── rules/
│   │   ├── stage1_gwas_vcf.smk        # Stage 1 rule
│   │   ├── stage2_overlaps.smk        # Stage 2 rule
│   │   ├── stage3_coloc_abf.smk       # Stage 3 rule
│   │   ├── stage4_coloc_susie.smk     # Stage 4 rule (SuSiE, optional)
│   │   └── stage5_postprocess.smk     # Stage 5 aggregation rule
│   └── scripts/
│       ├── 1_convert_gwas_to_vcf.R    # Stage 1 script
│       ├── 2_find_qtl_gwas_overlaps.R # Stage 2 script
│       ├── 3_run_coloc_abf.R          # Stage 3 script (core coloc logic)
│       ├── 4_run_coloc_susie.R        # Stage 4 script (SuSiE, template)
│       ├── 5_postprocess_coloc.R      # Stage 5 aggregation script
│       ├── 6_generate_plots.R         # Optional plots script
│       ├── Coloc_helper_functions.R   # perform_coloc() function
│       ├── coloc_helpers.R            # make_vcf() and other helpers
│       ├── data_loader.R              # Loading functions (load_gwas_signals, etc.)
│       └── qtl_processor.R           # Window creation and overlap functions
├── scripts/
│   ├── stage1_gwas_vcf.sh       # SLURM script for Stage 1
│   ├── stage2_overlaps.sh       # SLURM script for Stage 2
│   ├── stage3_coloc_abf.sh      # SLURM script for Stage 3
│   ├── stage3_only.sh           # Stage 3 only (skips Stage 2 re-run)
│   ├── stage4_aggregate.sh      # SLURM script for Stage 4 aggregation
│   ├── check_status.sh          # Pipeline status checker
│   ├── clean.sh                 # Cleanup utility
│   └── submit_pipeline.sh       # Main orchestrator
├── tests/
│   ├── test_stage1.R            # Stage 1 unit tests (mock + real data)
│   ├── test_stage2.R            # Stage 2 unit tests (mock + real data)
│   ├── test_stage3.R            # Stage 3 integration test
│   └── test_stage4.R            # Stage 4 integration test
├── markdown/
│   ├── PIPELINE_GUIDE.md        # This file
│   ├── PIPELINE_STATUS.md       # Current implementation status
│   └── ...                      # Tutorial and reference documents
└── logs/                        # SLURM job logs (stage*_<jobid>.out/err)
```

---

## 3. Environment Setup

### Conda environments

There are two conda environments:

**1. `coloc-pipeline`** — for running Snakemake and submitting jobs
```bash
source /home/itg/oleg.vlasovets/miniconda3/etc/profile.d/conda.sh
conda activate coloc-pipeline
snakemake --version   # Should print 9.16.0 or similar
```

**2. R coloc environment** — managed automatically by Snakemake via `--use-conda`
- Located at `.snakemake/conda/<hash>_/` after first build
- Contains: R, coloc 5.2.3, arrow 21.0.0, data.table, GenomicRanges, gwasvcf, etc.
- Built automatically on first `snakemake --use-conda` run (~10 min)

### One-time setup
```bash
cd /home/itg/oleg.vlasovets/projects/coloc-pipeline
conda activate coloc-pipeline

# Build R conda env (first time only — takes ~10 min)
snakemake --use-conda --conda-create-envs-only

# Verify pipeline config parses without errors
snakemake -n --dry-run 2>&1 | head -20
```

---

## 4. Input Data

### Required files

| File | Location | Format |
|------|----------|--------|
| GWAS summary stats | `/lustre/groups/itg/teams/zeggini/projects/GO2/GO2SummStats/ALL.MAIN/GO.FILTER.GW.final.meta.results.ALL.{trait}.FULL.MAFless0.01.Nupdated.txt.gz` | gzipped TSV with columns: CHR, POS, NEA, EA, CPTID, EAF, BETA, SE, P |
| Variant annotation | `/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/variant_ann_hg38.txt.gz` | gzipped TSV with columns: chr, pos, rsid, ref, alt, position |
| GWAS signals | `/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/GO2_index_signals_b38.csv` | CSV with independent GWAS signals, hg38 coordinates |
| eQTL permutations | `/lustre/groups/itg/.../eQTL_results/{tissue}/cov_df_*/cis_permutation/egenes_df.txt` | TSV with columns: phenotype_id, variant_id, qval, ... |
| eQTL nominal | `/lustre/groups/itg/.../eQTL_results/{tissue}/cov_df_*/cis_nominal/cis_qtl.cis_qtl_pairs.chr{N}.parquet` | Parquet, one file per chromosome |

### Tissue-specific paths (from `config.yaml`)

| Tissue | eQTL permutations subdirectory | QTL n |
|--------|-------------------------------|-------|
| `high_grade_cartilage` | `cov_df_cart_high_K35` | 115 |
| `low_grade_cartilage` | `cov_df_cart_low_K45` | 113 |
| `synovium` | `cov_df_synovium_K30` | 109 |
| `fat_pad` | `cov_df_fat_pad_K06` | 97 |

All input files are read-only from shared storage — do not copy or modify them.

---

## 5. Configuration

All parameters are in `config.yaml` at the repository root.

### Key parameters

```yaml
output_dir: "/lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results"

traits:
  - KNEE
  - TKR
  - ALLOA

tissues:
  - high_grade_cartilage
  - low_grade_cartilage
  - synovium
  - fat_pad

window_size: 1000000       # 1 Mb window around GWAS signals
qval_threshold: 0.05       # eGene significance threshold
cores: 4                   # Parallel cores for Stage 3

# Significance thresholds for Stage 4
pp4_threshold: 0.8         # PP4 >= 0.8 = strong colocalization
pp4_pp3_ratio: 2.0         # PP4/PP3 > 2 for moderate evidence

sample_sizes:
  KNEE:
    cases: 172256
    controls: 1144244
```

### How to add a new tissue

1. Add the tissue name to `tissues:` list in `config.yaml`
2. Add the tissue to `qtl_permutation_files:` with path to `egenes_df.txt`
3. Add the tissue to `qtl_nominal_dirs:` with path to nominal parquet directory
4. Add the tissue to `qtl_sample_sizes:` with sample size
5. Run Stage 2 and Stage 3 for the new tissue:
   ```bash
   sbatch scripts/stage2_overlaps.sh KNEE <new_tissue>
   # After Stage 2 completes:
   sbatch scripts/stage3_only.sh KNEE <new_tissue>
   # After Stage 3 completes:
   sbatch scripts/stage4_aggregate.sh KNEE
   ```

---

## 6. Running the Full Pipeline (All Stages Together)

For a trait/tissue combination where no outputs exist yet, use the main orchestrator:

```bash
cd /home/itg/oleg.vlasovets/projects/coloc-pipeline
conda activate coloc-pipeline
sbatch scripts/submit_pipeline.sh KNEE
```

This submits a SLURM job that runs Snakemake with `snakemake rule all`, processing
all 4 tissues in sequence.

**Expected total runtime for KNEE (all 4 tissues):**
- Stage 1: ~1–2 hours
- Stage 2: ~2–3 hours per tissue (~8 hours total, 4 tissues sequential)
- Stage 3: ~25–30 min per tissue
- Stage 4: ~5 min

**Important:** Snakemake re-checks provenance by default. If scripts changed, it may
re-run earlier stages. Use `--rerun-triggers mtime` (already in all SLURM scripts)
to prevent this once outputs are present.

---

## 7. Running Individual Stages

### Stage 1: GWAS VCF Conversion

**What it does:** Reads GWAS summary statistics (hg19), lifts over to hg38, applies
allele harmonization against the variant annotation, and writes a bgzipped, tabix-indexed
VCF file.

```bash
sbatch scripts/stage1_gwas_vcf.sh KNEE
```

**SLURM header:**
```bash
#SBATCH --job-name=coloc_stage1
#SBATCH --output=logs/stage1_%j.out
#SBATCH --error=logs/stage1_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
```

**Inputs:**
- `GO.FILTER.GW.final.meta.results.ALL.KNEE.FULL.MAFless0.01.Nupdated.txt.gz`
- `variant_ann_hg38.txt.gz`

**Outputs:**
- `results/gwas_vcf/KNEE_hg38.vcf.bgz`
- `results/gwas_vcf/KNEE_hg38.vcf.bgz.tbi`

**Verify completion:**
```bash
ls -lh /lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results/gwas_vcf/KNEE_hg38.vcf.bgz*
# Expected: ~184M VCF + .tbi index
```

---

### Stage 2: QTL-GWAS Overlap Detection

**What it does:** Reads GWAS signals, eQTL eGenes, and creates 1 Mb windows around
each. Finds genomic overlaps using GenomicRanges. Extracts GWAS association data from
the VCF for all overlapping windows. Performs allele harmonization.

```bash
# Single tissue:
sbatch scripts/stage2_overlaps.sh KNEE high_grade_cartilage

# Multiple tissues in one job:
sbatch scripts/stage2_overlaps.sh KNEE "low_grade_cartilage synovium fat_pad"
```

**SLURM header:**
```bash
#SBATCH --job-name=coloc_stage2
#SBATCH --output=logs/stage2_%j.out
#SBATCH --error=logs/stage2_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
```

**Inputs:**
- `results/gwas_vcf/KNEE_hg38.vcf.bgz` (Stage 1 output)
- `egenes_df.txt` for the tissue
- `GO2_index_signals_b38.csv`

**Outputs (per tissue):**
- `results/overlaps/KNEE.{tissue}.overlaps.rda` — overlap table + window objects
- `results/overlaps/KNEE.{tissue}.qtl_subset.txt` — gene IDs to process in Stage 3
- `results/overlaps/KNEE.{tissue}.gwas_data.rda` — harmonized GWAS associations

**Verify completion:**
```bash
wc -l /lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results/overlaps/KNEE.high_grade_cartilage.qtl_subset.txt
# Expected: ~1,158 genes (plus header line)
```

**Note on runtime:** Stage 2 reads the full GWAS VCF for all GWAS signal windows.
This takes ~2–3 hours per tissue. A checkpoint file is saved mid-run, so interrupted
jobs can resume without repeating VCF extraction.

---

### Stage 3: Coloc ABF Analysis

**What it does:** For each gene-signal pair in `qtl_subset.txt`, loads eQTL nominal
statistics from parquet files, performs allele flipping, then runs `coloc.abf()` to
compute PP0–PP4. Results are written to a flat text file.

```bash
# Standard (may re-run Stage 2 if Snakemake detects code changes):
sbatch scripts/stage3_coloc_abf.sh KNEE high_grade_cartilage

# Recommended (skips Stage 2 re-run if its outputs already exist):
sbatch scripts/stage3_only.sh KNEE high_grade_cartilage
```

**SLURM header:**
```bash
#SBATCH --job-name=coloc_stage3_high_grade_cartilage
#SBATCH --output=logs/stage3_high_grade_cartilage_%j.out
#SBATCH --error=logs/stage3_high_grade_cartilage_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
```

**Inputs:**
- `results/overlaps/KNEE.{tissue}.overlaps.rda`
- `results/overlaps/KNEE.{tissue}.qtl_subset.txt`
- `results/overlaps/KNEE.{tissue}.gwas_data.rda`
- `results/gwas_vcf/KNEE_hg38.vcf.bgz`
- eQTL nominal parquet files (`cis_qtl.cis_qtl_pairs.chr{1..22}.parquet`)

**Output:**
- `results/coloc_abf/KNEE.{tissue}.colocABF_results.txt`
  Columns: `cpg`, `cpg_pos`, `gwas_signal`, `gwas_lead_snp`, `PP0`, `PP1`, `PP2`, `PP3`, `PP4`, `nvariants`, `tissue`, `GWAS_ID`

**Verify completion:**
```bash
wc -l /lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results/coloc_abf/KNEE.high_grade_cartilage.colocABF_results.txt
# Expected: 2,287 + 1 header line
```

**Parallelism:** Stage 3 uses `mclapply` with `cores=4` (set in `config.yaml`).
Each gene-signal pair is tested independently.

---

### Stage 4: Aggregate Results

**What it does:** Reads all available `colocABF_results.txt` files across tissues,
combines them, calculates summary statistics, filters for significant colocalizations
(PP4 ≥ 0.8 OR PP4 > 0.6 AND PP4/PP3 > 2.0), and writes two output files.

```bash
sbatch scripts/stage4_aggregate.sh KNEE
```

**SLURM header:**
```bash
#SBATCH --job-name=coloc_stage4
#SBATCH --output=logs/stage4_%j.out
#SBATCH --error=logs/stage4_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
```

**Inputs:** All available `results/coloc_abf/KNEE.*.colocABF_results.txt`
(dynamic input — runs even if only some tissues have completed Stage 3)

**Outputs:**
- `results/results/KNEE_all_coloc_results.txt` — all tests, all tissues
- `results/results/KNEE_significant_coloc_results.txt` — PP4 ≥ 0.8 or moderate evidence

**Verify completion:**
```bash
wc -l /lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results/results/KNEE_*coloc_results.txt
head -2 /lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results/results/KNEE_significant_coloc_results.txt
```

---

## 8. Running All Tissues

### Recommended sequence

Stage 2 can process multiple tissues in one SLURM job (sequential within job):

```bash
# Step 1: Run Stage 2 for all 4 tissues at once
sbatch --job-name=coloc_stage2_all \
       --time=24:00:00 --mem=64G --cpus-per-task=4 \
       --partition=cpu_p --qos=cpu_normal \
       scripts/stage2_overlaps.sh KNEE "high_grade_cartilage low_grade_cartilage synovium fat_pad"

# Step 2: After Stage 2 completes, submit Stage 3 for each tissue
# (can submit all 4 simultaneously — they run in parallel as separate SLURM jobs)
sbatch scripts/stage3_only.sh KNEE high_grade_cartilage
sbatch scripts/stage3_only.sh KNEE low_grade_cartilage
sbatch scripts/stage3_only.sh KNEE synovium
sbatch scripts/stage3_only.sh KNEE fat_pad

# Step 3: After all Stage 3 jobs complete, run Stage 4
sbatch scripts/stage4_aggregate.sh KNEE
```

### Checking progress

```bash
# Queue status
squeue -u oleg.vlasovets

# Pipeline output status
bash scripts/check_status.sh KNEE

# Check individual tissue Stage 3 output
ls -lh /lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results/coloc_abf/KNEE.*.colocABF_results.txt
```

---

## 9. Output Files

### File descriptions

| File | Description |
|------|-------------|
| `results/gwas_vcf/KNEE_hg38.vcf.bgz` | GWAS summary stats as bgzipped VCF, hg38 |
| `results/gwas_vcf/KNEE_hg38.vcf.bgz.tbi` | Tabix index for VCF |
| `results/overlaps/KNEE.{tissue}.overlaps.rda` | R object: `overlaps`, `gwas_windows`, `qtl_windows` |
| `results/overlaps/KNEE.{tissue}.qtl_subset.txt` | Tab-separated: `gene_id` column, one gene per row |
| `results/overlaps/KNEE.{tissue}.gwas_data.rda` | R object: `GWAS_associations` data.table |
| `results/coloc_abf/KNEE.{tissue}.colocABF_results.txt` | Tab-separated coloc results, one row per gene-signal pair |
| `results/results/KNEE_all_coloc_results.txt` | All coloc tests combined across tissues |
| `results/results/KNEE_significant_coloc_results.txt` | Filtered to significant colocalizations only |
| `results/logs/*.log` | Per-rule R script output logs |

### Interpreting significant colocalizations

The significant results file contains rows where:
- `PP4 >= 0.8` — strong evidence of a shared causal variant, OR
- `PP4 > 0.6 AND PP4/PP3 > 2.0` — moderate evidence (high PP4 relative to PP3)

Key columns:
- `cpg` — gene ID (Ensembl)
- `gwas_signal` — GWAS lead SNP rsID for the signal window
- `PP4` — posterior probability of shared causal variant (higher = stronger coloc)
- `PP3` — posterior probability of two distinct causal variants
- `tissue` — tissue in which colocalization was detected
- `GWAS_ID` — trait (KNEE, TKR, ALLOA)
- `nvariants` — number of variants in the coloc test window

---

## 10. Tests

### Running the test suite

Tests are R scripts (same framework as production scripts) with two parts each:
- **Part A**: mock/synthetic data (always runs, no production inputs needed)
- **Part B**: real data validation (skipped automatically if outputs not present)

```bash
# Run on interactive node (not login node):
salloc --nodes=1 --cpus-per-task=4 --mem=16G --time=1:00:00 \
       --partition=interactive_cpu_p --qos=interactive_cpu_short --nice=10000

source /home/itg/oleg.vlasovets/miniconda3/etc/profile.d/conda.sh

# Activate the R conda env (contains all required packages)
conda activate /home/itg/oleg.vlasovets/projects/coloc-pipeline/.snakemake/conda/<hash>_

cd /home/itg/oleg.vlasovets/projects/coloc-pipeline

Rscript tests/test_stage1.R
Rscript tests/test_stage2.R
Rscript tests/test_stage3.R
Rscript tests/test_stage4.R
```

### What each test file covers

| Test file | Part A (mock) | Part B (real) |
|-----------|---------------|---------------|
| `test_stage1.R` | Allele flipping, variant ID format, biallelic filter, column validation, p-value/EAF ranges, sample size parameterization, missing file error | VCF exists and non-empty, tabix index present, VCF header format, log file |
| `test_stage2.R` | QTL window creation, GWAS window creation, overlap detection, output columns, chromosome validity, gene list extraction, empty QTL handling, GWAS harmonization, qtl_subset.txt format | Load overlaps.rda, column check, chromosome ranges, qtl_subset.txt, gwas_data.rda, varid_for_coloc format |
| `test_stage3.R` | — | Source helpers, load overlap data, column mapping, load GWAS data, load QTL gene list, load variant annotation, load QTL nominal data, allele flipping, output directories, run perform_coloc() |
| `test_stage4.R` | Read/combine files, required columns, PP value range, significant filtering, summary statistics, empty file handling, output file writing | Read real Stage 3 output, column check, PP ranges, significant filter, write aggregated output |

### How to add new tests
Follow the existing pattern in any test file:
1. Add a `section("Name")` call to label the test
2. Write assertions using `fail()` (stops with exit code 1) or `pass()` (continues)
3. Use `skip()` for conditional tests that depend on real data
4. Put synthetic-data tests in Part A, real-data tests in Part B

---

## 11. Troubleshooting

### Job fails immediately at Snakemake step

**Symptom:** Log shows `CommandNotFoundError: snakemake not found`

**Fix:** Activate the correct conda environment before running:
```bash
source /home/itg/oleg.vlasovets/miniconda3/etc/profile.d/conda.sh
conda activate coloc-pipeline
```
All SLURM scripts do this automatically. If submitting via `sbatch`, the script
handles activation internally.

---

### Stage 3 re-triggers Stage 2

**Symptom:** Stage 3 SLURM job unexpectedly runs Stage 2 again (2+ hours extra)

**Cause:** Snakemake detects code changes (default trigger) and re-runs upstream rules.

**Fix:** Use `scripts/stage3_only.sh` which passes `--rerun-triggers mtime`:
```bash
sbatch scripts/stage3_only.sh KNEE high_grade_cartilage
```

---

### "Unlocked working directory" warning

**Symptom:** `snakemake --unlock` needed before running

**Cause:** Previous run exited uncleanly leaving a `.snakemake/locks/` file.

**Fix:** All SLURM scripts run `snakemake --unlock` automatically before each job.
If running manually:
```bash
cd /home/itg/oleg.vlasovets/projects/coloc-pipeline
conda activate coloc-pipeline
snakemake --unlock
```

---

### Stage 2 fails with "No GWAS data extracted"

**Symptom:** `gwas_data.rda` saved as empty, Stage 3 fails with no data

**Cause:** GWAS VCF does not contain any variants in the signal windows for this tissue.

**Diagnosis:** Check the Stage 2 log:
```bash
cat /lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results/logs/overlaps_KNEE.{tissue}.log | grep -E "ERROR|variants|region"
```

---

### Stage 3 fails with "object 'overlap_df' not found"

**Cause:** Old version of `3_run_coloc_abf.R` expected variable name `overlap_df`,
but Stage 2 saves `overlaps`. Fixed in branch `o-feature/stage-3-implementation`.

**Fix:** Ensure you're on the correct branch:
```bash
git branch   # Should show * o-feature/stage-3-implementation
```

---

### Stage 4 shows wrong tissue/GWAS_ID in output

**Symptom:** `tissue` column contains the full filename, `GWAS_ID` column contains a path

**Cause:** Old regex in `5_postprocess_coloc.R`. Fixed — the current version uses:
```r
tissue  <- sub("^[^.]+\\.(.+)\\.colocABF_results\\.txt$", "\\1", basename(f))
gwas_id <- sub("^([^.]+)\\..*", "\\1", basename(f))
```

---

### How to resume a failed run

Most SLURM scripts include `--rerun-incomplete` which resumes partial Snakemake outputs.
Stage 2 additionally writes a checkpoint file (`.checkpoint_raw.rds`) that lets it
skip the slow VCF extraction step if the job is re-submitted.

```bash
# Re-submit the failed stage — Snakemake will resume automatically
sbatch scripts/stage3_only.sh KNEE high_grade_cartilage
```

---

### Checking SLURM job logs

```bash
# List recent job logs
ls -lt /home/itg/oleg.vlasovets/projects/coloc-pipeline/logs/ | head -20

# View stdout for a specific job
cat logs/stage3_34965770.out

# View stderr (Snakemake progress + R stderr)
cat logs/stage3_34965770.err

# Check detailed R log from Snakemake (more verbose than SLURM logs)
cat /lustre/scratch/users/oleg.vlasovets/coloc-pipeline/results/logs/coloc_abf_KNEE.high_grade_cartilage.log
```

---

### Known limitations

- **Stage 2 takes ~2–3 hours per tissue** due to full VCF scanning. This is a known
  bottleneck; the checkpoint file mitigates re-runs.
- **Stage 3 produces no output for some tissues** if `qtl_subset.txt` is empty
  (no overlapping genes found). This is expected for `fat_pad` which has few eGenes.
- **Snakemake provenance warnings** (`Rules with missing metadata`) are expected
  for jobs run with `--rerun-triggers mtime`. They do not cause re-runs.
- **Stage 4 (SuSiE)** is a template only — the R script is not implemented.

---

## 12. Pipeline Status

See [PIPELINE_STATUS.md](PIPELINE_STATUS.md) for the current implementation state,
completed job IDs, and result counts.

**As of 2026-04-02:**
- Stage 1 ✓ (KNEE, 11.3M SNPs)
- Stage 2 ✓ (KNEE, all 4 tissues)
- Stage 3 ✓ (KNEE/high_grade_cartilage, 2,287 tests; others in progress)
- Stage 4 ✓ (KNEE/high_grade_cartilage: 54 significant colocalizations)
