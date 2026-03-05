from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.plotting import (
    plot_agreement_distribution,
    plot_ancestry_stacked_bar,
    plot_label_confusion_heatmap,
    plot_tract_length_distribution,
    plot_tool_comparison_overview,
)
from src.utils import ensure_dir


def parse_args():
    parser = argparse.ArgumentParser(description="Generate benchmark figures")
    parser.add_argument("--config", default="configs/dataset.yaml", help="Unused currently, kept for interface")
    parser.add_argument("--metrics-dir", default="results/metrics", help="Input metrics directory")
    parser.add_argument("--outdir", default="results/figures", help="Output figure directory")
    return parser.parse_args()


def main():
    args = parse_args()
    metrics_dir = Path(args.metrics_dir)
    outdir = Path(args.outdir)
    ensure_dir(outdir)

    site_agreement = pd.read_csv(metrics_dir / "site_agreement.csv")
    tract_segments = pd.read_csv(metrics_dir / "tract_segments.tsv", sep="\t")
    ancestry_props = pd.read_csv(metrics_dir / "ancestry_props.csv")
    label_confusion_rates = pd.read_csv(metrics_dir / "label_confusion_matrix_rates.csv", index_col=0)
    tool_summary = pd.read_csv(metrics_dir / "tool_comparison_summary.csv")
    comparison_dashboard = pd.read_csv(metrics_dir / "comparison_dashboard.csv")

    plot_agreement_distribution(site_agreement, outdir / "agreement_distribution.png")
    plot_tract_length_distribution(tract_segments, outdir / "tract_length_distribution.png")
    plot_ancestry_stacked_bar(ancestry_props, outdir / "ancestry_proportion_stacked.png")
    plot_label_confusion_heatmap(label_confusion_rates, outdir / "label_confusion_heatmap.png")
    plot_tool_comparison_overview(tool_summary, comparison_dashboard, outdir / "tool_comparison_overview.png")

    print(f"figures saved")


if __name__ == "__main__":
    main()
