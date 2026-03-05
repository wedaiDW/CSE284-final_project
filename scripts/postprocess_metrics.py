from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.io_gnomix import load_gnomix_per_site
from src.io_rfmix import load_rfmix_per_site
from src.metrics import compute_ancestry_proportions, compute_site_agreement
from src.tracts import compute_tract_segments, summarize_tracts
from src.utils import apply_label_map, ensure_dir, load_yaml


def parse_args():
    parser = argparse.ArgumentParser(description="Postprocess LAI outputs and compute metrics")
    parser.add_argument("--config", default="configs/dataset.yaml", help="Path to dataset config")
    parser.add_argument("--gnomix-dir", default="results/gnomix", help="Gnomix output root")
    parser.add_argument("--rfmix-dir", default="results/rfmix", help="RFMix output root")
    parser.add_argument("--outdir", default="results/metrics", help="Metrics output directory")
    return parser.parse_args()


def maybe_attach_population(df, panel_tsv):
    if not panel_tsv.exists() or df.empty:
        return df

    panel = pd.read_csv(panel_tsv, sep="\t")
    if "sample" not in panel.columns or "pop" not in panel.columns:
        return df

    pop_map = panel[["sample", "pop"]].drop_duplicates()
    return df.merge(pop_map, left_on="sample_id", right_on="sample", how="left").drop(columns=["sample"])


def build_tool_comparison_summary(tracts_df, ancestry_props):
    if tracts_df.empty:
        return pd.DataFrame(
            columns=[
                "tool",
                "n_samples",
                "n_tracts",
                "mean_tract_length_bp",
                "median_tract_length_bp",
                "mean_tracts_per_sample",
            ]
        )

    all_labels = sorted(ancestry_props["anc_label"].dropna().astype(str).unique().tolist())
    rows: list[dict[str, float | int | str]] = []
    for tool in sorted(tracts_df["tool"].dropna().astype(str).unique()):
        tdf = tracts_df[tracts_df["tool"] == tool].copy()
        tracts_per_sample = tdf.groupby("sample_id").size()
        row: dict[str, float | int | str] = {
            "tool": tool,
            "n_samples": int(tdf["sample_id"].nunique()),
            "n_tracts": int(len(tdf)),
            "mean_tract_length_bp": float(tdf["tract_length_bp"].mean()),
            "median_tract_length_bp": float(tdf["tract_length_bp"].median()),
            "mean_tracts_per_sample": float(tracts_per_sample.mean()),
        }

        prop_df = ancestry_props[ancestry_props["tool"] == tool].copy()
        if not prop_df.empty:
            pivot = prop_df.pivot_table(
                index="sample_id",
                columns="anc_label",
                values="prop_sites",
                aggfunc="sum",
                fill_value=0.0,
            )
            for anc in all_labels:
                row[f"mean_prop_{anc}"] = float(pivot[anc].mean()) if anc in pivot.columns else 0.0
        rows.append(row)

    return pd.DataFrame(rows).sort_values("tool").reset_index(drop=True)


def build_comparison_dashboard(site_summary, tool_summary):
    sample_rows = site_summary[site_summary["sample_id"] != "__overall__"].copy()
    overall_rows = site_summary[site_summary["sample_id"] == "__overall__"].copy()

    row: dict[str, float | int] = {
        "n_samples": int(len(sample_rows)),
        "overall_agreement_rate": float(overall_rows["agreement_rate"].iloc[0]) if not overall_rows.empty else float("nan"),
        "mean_sample_agreement_rate": float(sample_rows["agreement_rate"].mean()) if not sample_rows.empty else float("nan"),
        "median_sample_agreement_rate": float(sample_rows["agreement_rate"].median()) if not sample_rows.empty else float("nan"),
        "p10_sample_agreement_rate": float(sample_rows["agreement_rate"].quantile(0.10))
        if not sample_rows.empty
        else float("nan"),
        "p90_sample_agreement_rate": float(sample_rows["agreement_rate"].quantile(0.90))
        if not sample_rows.empty
        else float("nan"),
    }

    if set(tool_summary["tool"].tolist()) >= {"gnomix", "rfmix"}:
        g = tool_summary[tool_summary["tool"] == "gnomix"].iloc[0]
        r = tool_summary[tool_summary["tool"] == "rfmix"].iloc[0]
        row["delta_mean_tract_length_bp_gnomix_minus_rfmix"] = float(
            g["mean_tract_length_bp"] - r["mean_tract_length_bp"]
        )
        row["delta_median_tract_length_bp_gnomix_minus_rfmix"] = float(
            g["median_tract_length_bp"] - r["median_tract_length_bp"]
        )

        anc_cols = sorted([c for c in tool_summary.columns if c.startswith("mean_prop_")])
        for col in anc_cols:
            row[f"delta_{col}_gnomix_minus_rfmix"] = float(g.get(col, 0.0) - r.get(col, 0.0))

    return pd.DataFrame([row])


def main():
    args = parse_args()
    cfg = load_yaml(args.config)

    label_map = cfg.get("ancestry_label_map", {})
    panel_tsv = Path(cfg.get("panel_tsv", "./data/raw/panel.tsv"))

    outdir = Path(args.outdir)
    ensure_dir(outdir)

    gnomix = load_gnomix_per_site(Path(args.gnomix_dir))
    rfmix = load_rfmix_per_site(Path(args.rfmix_dir))

    gnomix["anc_label"] = apply_label_map(gnomix["anc_label"], label_map)
    rfmix["anc_label"] = apply_label_map(rfmix["anc_label"], label_map)

    gnomix = maybe_attach_population(gnomix, panel_tsv)
    rfmix = maybe_attach_population(rfmix, panel_tsv)

    gnomix.to_csv(outdir / "gnomix_per_site.tsv", sep="\t", index=False)
    rfmix.to_csv(outdir / "rfmix_per_site.tsv", sep="\t", index=False)

    merged, site_summary = compute_site_agreement(gnomix, rfmix)
    merged.to_csv(outdir / "site_agreement_detailed.tsv", sep="\t", index=False)
    site_summary.to_csv(outdir / "site_agreement.csv", index=False)

    label_conf_counts = pd.crosstab(
        merged["anc_label_gnomix"],
        merged["anc_label_rfmix"],
        rownames=["gnomix_label"],
        colnames=["rfmix_label"],
    )
    label_conf_counts.to_csv(outdir / "label_confusion_matrix_counts.csv")
    total_pairs = float(label_conf_counts.to_numpy().sum())
    label_conf_rates = label_conf_counts / total_pairs if total_pairs > 0 else label_conf_counts
    label_conf_rates.to_csv(outdir / "label_confusion_matrix_rates.csv")

    gnomix_tracts = compute_tract_segments(gnomix)
    gnomix_tracts["tool"] = "gnomix"
    rfmix_tracts = compute_tract_segments(rfmix)
    rfmix_tracts["tool"] = "rfmix"

    all_tracts = pd.concat([gnomix_tracts, rfmix_tracts], ignore_index=True)
    all_tracts.to_csv(outdir / "tract_segments.tsv", sep="\t", index=False)

    tract_stats = summarize_tracts(all_tracts)
    tract_stats.to_csv(outdir / "tract_stats.csv", index=False)

    gnomix_prop = compute_ancestry_proportions(gnomix)
    gnomix_prop["tool"] = "gnomix"
    rfmix_prop = compute_ancestry_proportions(rfmix)
    rfmix_prop["tool"] = "rfmix"

    ancestry_props = pd.concat([gnomix_prop, rfmix_prop], ignore_index=True)
    ancestry_props.to_csv(outdir / "ancestry_props.csv", index=False)

    tool_summary = build_tool_comparison_summary(all_tracts, ancestry_props)
    tool_summary.to_csv(outdir / "tool_comparison_summary.csv", index=False)

    dashboard = build_comparison_dashboard(site_summary, tool_summary)
    dashboard.to_csv(outdir / "comparison_dashboard.csv", index=False)

    print(f"[ok] wrote metrics to {outdir}")
    print(f"[ok] site agreement rows: {len(site_summary)}")
    print(f"[ok] tract segments rows: {len(all_tracts)}")


if __name__ == "__main__":
    main()
