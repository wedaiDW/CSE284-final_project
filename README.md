# lai-benchmark

A benchmark framework for comparing Local Ancestry Inference tools Gnomix vs RFMix. This repo evaluates and generates visualation for both tool on agreement and runtime.

You will have tool manually install RFMix and Gnomix tool to run this benchmark.

```tsv
sample_id\tchrom\tpos\tanc_label
NAxxxx\t22\t16050075\tAFR
```

## Quickstart

```bash
# If conda is available
conda env create -f environment.yml
conda activate lai-benchmark

# or: pip install -r requirements.txt

python scripts/run_pipeline.py

# data root is optional
python scripts/run_pipeline.py \
  --data-root /default/directory/data/raw/if/leave/blank \
  --chromosomes 22
```

## Or Step-by-Step

```bash
bash scripts/check_env.sh
python scripts/fetch_or_point_data.py --config configs/dataset.yaml
python scripts/make_splits.py --config configs/dataset.yaml --overwrite
python scripts/prepare_inputs.py --config configs/dataset.yaml
bash scripts/benchmark_all.sh --config configs/dataset.yaml --tools configs/tools.yaml
python scripts/make_figures.py --config configs/dataset.yaml
```

## Config Files

- `configs/dataset.yaml`
  - chromosomes, reference populations, admixed test populations
  - sample sizes and data paths
- `configs/tools.yaml`
  - gnomix / rfmix command templates, thread count, and options

## Outputs

- `results/benchmark/`: `/usr/bin/time -v` logs (runtime/memory)
- `results/metrics/`: metric CSV files
- `results/figures/`: final plots for report/presentation

## Implementation Note

- Thin shell wrappers were removed.
- Keep using `bash scripts/check_env.sh` and `bash scripts/benchmark_all.sh`.
- Use Python modules directly for data check/prep/tool runs.
