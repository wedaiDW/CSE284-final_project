# Local Ancestry Benchmark Pipeline

This repository benchmarks two local ancestry inference (LAI) tools, Gnomix and RFMix, on the same input data and summarizes their behavior with agreement, ancestry composition, tract-level statistics, and runtime logs.

## Motivation

The main goal of this project is to make side-by-side LAI benchmarking easy to reproduce. Instead of only running one tool and manually inspecting outputs, this pipeline standardizes:

- sample selection from a population panel
- chromosome-specific input preparation
- tool execution for Gnomix and RFMix
- postprocessing into a common per-site format
- metric computation and figure generation

This makes it easier to answer practical questions such as:

- Do the two tools agree on local ancestry labels?
- Do they produce similar ancestry proportions per sample?
- Are the inferred tracts similarly smooth or fragmented?
- What are the runtime and operational tradeoffs?

## Included Demo Assets

The repository already includes a small `chr22` demo so reviewers can run something immediately:

- input demo data: [`data/README.md`](data/README.md)
- sample outputs: [`results/README.md`](results/README.md)

The default config [`configs/dataset.yaml`](configs/dataset.yaml) is already set to chromosome 22, so the fastest smoke test is the bundled demo run.

## Quickstart

### 1. Create the environment

```bash
conda env create -f environment.yml
conda activate lai-benchmark
```

If you prefer a specific Python binary, all thin shell wrappers support `PYTHON=/path/to/python`.

### 2. Run the bundled `chr22` demo

```bash
bash scripts/run_chr22_demo.sh
```

Or explicitly:

```bash
python scripts/run_pipeline.py \
  --config configs/dataset.yaml \
  --tools configs/tools.yaml \
  --chromosomes 22 \
  --resolved-config configs/dataset.runtime.chr22.yaml
```

### 3. Point the pipeline to your own raw data

```bash
python scripts/run_pipeline.py \
  --data-root /path/to/raw \
  --chromosomes 22
```

If your panel or VCFs are not under the default layout, you can override them directly:

```bash
python scripts/run_pipeline.py \
  --panel-tsv /path/to/panel.tsv \
  --vcf-pattern '/path/to/1000G_chr{chrom}_pruned.vcf.gz' \
  --chromosomes 22
```

## Expected Raw Data Layout

The default layout is:

```text
data/raw/
├── meta/panel.tsv
├── 1000G_chr22_pruned.vcf.gz
└── 1000G_chr22_pruned.vcf.gz.tbi
```

For multiple chromosomes:

```text
<data-root>/
├── meta/panel.tsv
├── 1000G_chr1_pruned.vcf.gz
├── 1000G_chr1_pruned.vcf.gz.tbi
├── ...
├── 1000G_chr22_pruned.vcf.gz
└── 1000G_chr22_pruned.vcf.gz.tbi
```

If data are missing, [`scripts/fetch_or_point_data.py`](scripts/fetch_or_point_data.py) now prints the expected layout and example commands instead of only saying “place data at the data root”.

## Thin Shell Wrappers

The repository includes thin wrappers again for reviewers who want a short entry point:

- [`scripts/run_chr22_demo.sh`](scripts/run_chr22_demo.sh): full end-to-end `chr22` demo run
- [`scripts/benchmark_chr22.sh`](scripts/benchmark_chr22.sh): prepare inputs and run benchmark stages
- [`scripts/rebuild_reports.sh`](scripts/rebuild_reports.sh): recompute metrics and figures from existing tool outputs

These wrappers intentionally stay thin and pass through to the Python modules, so the repo still preserves a clean `configs/` + `scripts/` + `src/` structure.

## Step-by-Step Commands

```bash
bash scripts/check_env.sh
python scripts/fetch_or_point_data.py --config configs/dataset.yaml
python scripts/make_splits.py --config configs/dataset.yaml --overwrite
python scripts/prepare_inputs.py --config configs/dataset.yaml
bash scripts/benchmark_all.sh --config configs/dataset.yaml --tools configs/tools.yaml
python scripts/make_figures.py --config configs/dataset.yaml
```

## Project Structure

- `configs/`: dataset and tool YAML configuration
- `scripts/`: orchestration, wrappers, and command-line entrypoints
- `src/`: reusable loading, metrics, tract, and plotting logic
- `data/`: bundled demo inputs plus prepared intermediates
- `results/`: benchmark logs, tool outputs, metrics, and figures

## Outputs

After a successful run, the key artifacts are:

- `results/gnomix/chr*/per_site.tsv`
- `results/rfmix/chr*/per_site.tsv`
- `results/metrics/comparison_dashboard.csv`
- `results/metrics/tool_comparison_summary.csv`
- `results/figures/*.png`
- `results/benchmark/*.time.txt`

## Analysis Snapshot

The current bundled `chr22` example already demonstrates what this project measures.

From the included `chr22` run:

- `50` PEL test samples were compared against `CEU / CHB / YRI` reference panels
- overall Gnomix vs RFMix agreement was `0.9521`
- Gnomix produced `88` tracts with mean tract length `8.04 Mb`
- RFMix produced `93` tracts with mean tract length `7.61 Mb`

This project intentionally evaluates the tools from multiple angles:

- site-level agreement
- ancestry proportion per sample
- tract count and tract length distribution
- ancestry label confusion matrix
- runtime logs

## Important Limitation

For the bundled 1000 Genomes `PEL` benchmark, there is no official per-site local ancestry ground truth. That means the default evaluation is primarily agreement-based rather than accuracy-based.

This is still useful for operational benchmarking, but it is not the same as measuring biological accuracy against truth labels.

## Remaining Work

- add a simulated admixed benchmark with known ground-truth labels
- report chromosome-by-chromosome and ancestry-specific breakdowns
- add a tiny automated smoke test around the bundled `chr22` demo
- expand the report section with more biological interpretation of discordant samples

## Notes

- `environment.yml` is the recommended setup path for reproducibility
- if your local default Python has environment issues, run wrappers as `PYTHON=/path/to/python bash scripts/run_chr22_demo.sh`
- the repo includes both code and example outputs so reviewers can inspect results even before rerunning the pipeline
