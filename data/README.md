# Data Directory

This folder contains both bundled demo inputs and intermediate files produced by the pipeline.

## Bundled Demo Inputs

The repository includes a small `chr22` demo dataset so the pipeline can be smoke-tested without first downloading the full benchmark inputs:

- `raw/meta/panel.tsv`: sample-to-population metadata
- `raw/1000G_chr22_pruned.vcf.gz`
- `raw/1000G_chr22_pruned.vcf.gz.tbi`

These files are enough to run the default `chr22` configuration in [`../configs/dataset.yaml`](../configs/dataset.yaml).

## Generated During Runtime

- `sample_lists/`: reference and test sample lists
- `prepared/`: prepared chromosome-specific VCFs used by Gnomix and RFMix

## If You Want to Use Your Own Data

Point the pipeline to another raw data root:

```bash
python scripts/run_pipeline.py --data-root /path/to/raw --chromosomes 22
```

Expected layout:

```text
<data-root>/
├── meta/panel.tsv
├── 1000G_chr22_pruned.vcf.gz
└── 1000G_chr22_pruned.vcf.gz.tbi
```
