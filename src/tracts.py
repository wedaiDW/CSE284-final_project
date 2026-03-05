from __future__ import annotations

import pandas as pd


def compute_tract_segments(df):
    required = {"sample_id", "chrom", "pos", "anc_label"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"missing columns for tract computation: {sorted(missing)}")

    rows = []

    for (sample_id, chrom), sub in df.sort_values(["sample_id", "chrom", "pos"]).groupby(
        ["sample_id", "chrom"],
        sort=False,
    ):
        positions = sub["pos"].tolist()
        labels = sub["anc_label"].tolist()

        if not positions:
            continue

        start_idx = 0
        for i in range(1, len(labels) + 1):
            is_break = i == len(labels) or labels[i] != labels[start_idx]
            if not is_break:
                continue

            tract_start = int(positions[start_idx])
            tract_end = int(positions[i - 1])
            n_sites = i - start_idx
            rows.append(
                {
                    "sample_id": sample_id,
                    "chrom": chrom,
                    "anc_label": labels[start_idx],
                    "tract_start": tract_start,
                    "tract_end": tract_end,
                    "tract_length_bp": tract_end - tract_start + 1,
                    "n_sites": n_sites,
                }
            )
            start_idx = i

    return pd.DataFrame(rows)


def summarize_tracts(tracts_df):
    if tracts_df.empty:
        return pd.DataFrame(
            columns=[
                "tool",
                "anc_label",
                "n_tracts",
                "mean_length_bp",
                "median_length_bp",
                "mean_sites",
                "median_sites",
            ]
        )

    grouping = ["anc_label"]
    if "tool" in tracts_df.columns:
        grouping = ["tool", "anc_label"]

    summary = (
        tracts_df.groupby(grouping, as_index=False)
        .agg(
            n_tracts=("tract_length_bp", "size"),
            mean_length_bp=("tract_length_bp", "mean"),
            median_length_bp=("tract_length_bp", "median"),
            mean_sites=("n_sites", "mean"),
            median_sites=("n_sites", "median"),
        )
        .sort_values(grouping)
    )
    return summary
