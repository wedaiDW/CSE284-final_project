# Results Directory

This folder stores example outputs and regenerated artifacts from the benchmark pipeline.

## What Reviewers Should Look At First

- `metrics/comparison_dashboard.csv`: one-row summary of the current benchmark
- `metrics/tool_comparison_summary.csv`: per-tool tract and ancestry summaries
- `figures/tool_comparison_overview.png`: compact overview figure
- `figures/ancestry_proportion_stacked.png`: ancestry composition summary

## Example Included Output

The bundled `chr22` demo has been run and the repository includes example outputs for:

- `results/gnomix/chr22/`
- `results/rfmix/chr22/`
- `results/metrics/`
- `results/figures/`

This means reviewers can inspect both raw tool outputs and postprocessed summaries even before rerunning the code.

## Rebuilding the Reports

If `per_site.tsv` files already exist and you only want to refresh the analysis products:

```bash
bash scripts/rebuild_reports.sh
```
