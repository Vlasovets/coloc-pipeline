## Colocalization Pipeline - Quick Reference

### Run Complete Pipeline

```bash
# Cluster mode (recommended)
sbatch submit_pipeline.sh

# Local mode (testing)
bash run_pipeline_local.sh
```

### Check Progress

```bash
bash check_pipeline_status.sh
```

### Resume After Failure

```bash
# Just rerun the same command - checkpoints are automatic!
sbatch submit_pipeline.sh
```

### Common Commands

```bash
# Preview what will run (dry run)
snakemake --dry-run -n

# Force rerun specific stage
snakemake --forcerun run_coloc_abf --cores 8

# Force rerun everything
snakemake --forceall --cores 8

# Clean incomplete jobs
snakemake --rerun-incomplete

# View logs
tail -f logs/snakemake_*.out
tail -f results/logs/coloc_abf_KNEE.high_grade_cartilage.log
```

### File Locations

```
results/gwas_vcf/       → Stage 1 outputs (GWAS VCFs)
results/overlaps/       → Stage 2 outputs (QTL-GWAS overlaps)
results/coloc_abf/      → Stage 3 outputs (Coloc results)
results/results/        → Stage 5 outputs (Aggregated)
logs/                   → Snakemake logs
results/logs/           → Individual stage logs
```

### Checkpoints Explained

✓ **Automatic**: Snakemake tracks completed stages via output files
✓ **Resume**: Rerun after failure picks up where it stopped
✓ **Skip**: Completed stages with valid outputs are skipped
✓ **Parallel**: Independent jobs (tissues) run simultaneously

### Need Help?

- Full docs: `less CHECKPOINT_SYSTEM.md`
- Workflow: `less WORKFLOW_DOCUMENTATION.md`
- Config: `vim config.yaml`
