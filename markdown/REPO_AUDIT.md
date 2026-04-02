# Repo Audit Report
**Date:** 2026-04-02  
**Branch audited:** `o-feature/stage-3-implementation`  
**Purpose:** Identify redundant files, branches, large binaries, untracked docs, and gitignore gaps before cleanup.

---

## 1. Branch Inventory

| Branch | Tip commit | Tracked remote | Status |
|--------|-----------|---------------|--------|
| `o-feature/stage-3-implementation` *(current)* | `7aa2223` | no | **Active work — not yet merged** |
| `main` | `f5c132d` | `origin/main` | Production baseline |
| `dev` | `52882b1` | `origin/dev` | Stale — 6 commits behind current branch |
| `o-feature/stage-3` | `5b75b51` | no | Superseded (see below) |
| `o-test-rule` | `52882b1` | no | Superseded (see below) |

### Commits ahead of main (current branch)
```
7aa2223 docs: replace all personal paths with portable variables
4d68c81 docs: add portable environment setup (conda + apptainer), runtimes
2255359 test: add Stage 1 and Stage 2 tests; docs: comprehensive pipeline guide
7cf7935 feat: fix tissue/GWAS_ID regex and mark Stage 3+4 complete
084030b feat: implement Stage 4 aggregate results
64e9573 feat: implement and test stage 3 - coloc ABF analysis
52882b1 feat(stage2): Fix VCF extraction and add local scratch optimization
47abfa6 debug stage 2
648407d Combine tests for comprehensive R environment tests (#20)
02a80a1 Merge pull request #16 from Vlasovets/o-feature/github-workflow
b692445 update the workflow
```

### Branch cleanup candidates

**`o-feature/stage-3`** — has 2 unique commits not in `o-feature/stage-3-implementation`:
- `d2651dc add tests for stage 2` — older, superseded by commit `2255359` in current branch (proper test_stage1.R / test_stage2.R with all 29 passing tests)
- `5b75b51 update .gitignore` — may contain useful .gitignore changes not yet merged; **verify before deleting**

**`o-test-rule`** — tip is `52882b1`, same commit as `dev`. No unique content.  
**Recommendation:** safe to delete `o-test-rule`. Review `o-feature/stage-3` .gitignore diff first.

---

## 2. Untracked Files

### Large binaries (should be gitignored)

| Path | Size | Action |
|------|------|--------|
| `coloc-pipeline.sif` | 946 MB | Add to `.gitignore` (`*.sif`) |
| `containers/r_coloc/Miniforge3-Linux-x86_64.sh` | ~100 MB | Add to `.gitignore` (`**/Miniforge3*.sh`) |
| `containers/r_coloc_reuse/Miniforge3-Linux-x86_64.sh` | ~100 MB | Same |
| `containers/r_coloc_reuse/conda/` | part of 101 MB | Installed conda env — gitignore the whole `containers/` dir |

### Build artifacts / container work

| Path | Size | Action |
|------|------|--------|
| `containers/r_coloc/` | 101 MB | Intermediate build dir — gitignore |
| `containers/r_coloc_reuse/` | 101 MB | Intermediate build dir — gitignore |
| `coloc-pipeline.def` | small | Definition file — **consider tracking** (reproducible container builds) |
| `scripts/build_container.sh` | small | Build helper — **consider tracking** alongside `.def` |

### Untracked documentation (markdown/)

These docs are untracked but **not** gitignored. They will be committed if someone runs `git add markdown/`.

| File | Content | Recommendation |
|------|---------|---------------|
| `markdown/BUILD_CONTAINER.md` | Container build instructions | Keep untracked or add to .gitignore |
| `markdown/BUILD_QUICK_START.md` | Quick-start for container | May overlap with PIPELINE_GUIDE.md; review |
| `markdown/GIT_SNAKEMAKE_CONTAINER_TUTORIAL.md` | Tutorial | Tutorial content, not pipeline ops — keep local |
| `markdown/README.md` | Unknown | Review contents |
| `markdown/TUTORIAL_QUICK_REFERENCE.md` | Tutorial reference | Tutorial content — keep local |

**Tracked docs (correctly):** `PIPELINE_GUIDE.md`, `PIPELINE_STATUS.md`  
**Gitignored (correctly):** `CONTINUE_PIPELINE.md`, `REFACTORING_SUMMARY.md`

### Other untracked files

| Path | Action |
|------|--------|
| `Snakefile_tutorial` | Tutorial artifact — add to .gitignore or delete |
| `workflow/rules/test.smk` | Development test rule — review; add to .gitignore or track if useful |

---

## 3. Tracked Files — Redundancy Analysis

### `src/` vs `workflow/scripts/`

`src/` (9 files, tracked) contains the **original pre-Snakemake scripts** from before the pipeline was built:

| `src/` file | Relationship to `workflow/scripts/` |
|-------------|-------------------------------------|
| `src/convert_gwas_to_vcf_hg38.R` | Predecessor of `workflow/scripts/1_convert_gwas_to_vcf.R` |
| `src/3_run_coloc_abf.R` | Predecessor of `workflow/scripts/3_run_coloc_abf.R` |
| `src/5_postprocess_coloc.R` | Predecessor of `workflow/scripts/5_postprocess_coloc.R` |
| `src/Coloc_helper_functions.R` | Same file now in `workflow/scripts/Coloc_helper_functions.R` |
| `src/1_make_GO2_b38_vcf.R`, `src/1_make_vcf.ipynb` | Older GWAS prep approaches |
| `src/2_run_coloc_abf.ipynb`, `src/4_run_coloc_susie.R` | Exploratory notebooks/scripts |
| `src/README.md` | `src/`-specific readme |

`src/` is retained per `# !src/` comment in `.gitignore`. Its role is historical reference.  
**Recommendation:** Keep `src/` as-is (it is small and documents provenance). Optionally rename to `src_legacy/` to make intent clearer.

### `logs/` and `jobs/` directories

Both exist on disk but are **correctly gitignored** — `git ls-files` returns nothing for them.  
`logs/` has 30 files (SLURM `.out`/`.err` from past runs). These are local only.

### `workflow/scripts/` — potential internal redundancy

| File | Note |
|------|------|
| `Coloc_helper_functions.R` (37 KB) | Main coloc helper library |
| `coloc_helpers.R` (38 KB) | Nearly same size — possible duplicate or older version |
| `data_loader.R` | May overlap with helpers in `Coloc_helper_functions.R` |

**Recommendation:** Compare `coloc_helpers.R` vs `Coloc_helper_functions.R` — if one is unused, remove it.

---

## 4. Gitignore Gaps

The following untracked paths would be accidentally staged by `git add .`:

```
# Add these to .gitignore:
*.sif                          # Apptainer/Singularity images
containers/                    # Intermediate container builds
Miniforge3-Linux-x86_64.sh    # (redundant if containers/ is ignored)
Snakefile_tutorial             # Tutorial artifact
workflow/rules/test.smk        # Dev-only rule (or track it)

# Consider (if keeping docs local):
markdown/BUILD_CONTAINER.md
markdown/BUILD_QUICK_START.md
markdown/GIT_SNAKEMAKE_CONTAINER_TUTORIAL.md
markdown/README.md
markdown/TUTORIAL_QUICK_REFERENCE.md
```

---

## 5. Summary Table

| Category | Count | Action needed |
|----------|-------|---------------|
| Branches to delete | 2 (`o-test-rule`, `o-feature/stage-3`) | Verify .gitignore diff on stage-3 first |
| Large untracked binaries to gitignore | 3+ | Add `*.sif`, `containers/` to .gitignore |
| Build artifacts to gitignore | 2 dirs | `containers/r_coloc`, `r_coloc_reuse` |
| Untracked docs to decide on | 5 | Add to .gitignore or track |
| Untracked scripts to decide on | 2 (`coloc-pipeline.def`, `scripts/build_container.sh`) | Track them |
| Possibly redundant scripts | 1 (`coloc_helpers.R`) | Verify vs `Coloc_helper_functions.R` |
| Tracked legacy code | `src/` (9 files) | Keep as-is or clarify with comment |

---

## 6. Recommended Cleanup Order

1. **Check `o-feature/stage-3` .gitignore diff** — extract any useful changes, then delete branch
2. **Delete `o-test-rule`** — no unique content
3. **Update `.gitignore`** — add `*.sif`, `containers/`, `Snakefile_tutorial`, untracked markdown docs
4. **Track `coloc-pipeline.def` + `scripts/build_container.sh`** — reproducibility
5. **Investigate `coloc_helpers.R`** — delete if unused duplicate of `Coloc_helper_functions.R`
6. **Merge current branch to main** — 11 commits of stage 3/4 work + docs + tests pending merge
