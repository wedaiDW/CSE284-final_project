from __future__ import annotations

import pandas as pd


def compute_site_agreement(gnomix, rfmix):
    key_cols = ["sample_id", "chrom", "pos"]

    g = gnomix[key_cols + ["anc_label"]].drop_duplicates(key_cols)
    r = rfmix[key_cols + ["anc_label"]].drop_duplicates(key_cols)

    merged = g.merge(r, on=key_cols, suffixes=("_gnomix", "_rfmix"), how="inner")
    if merged.empty:
        raise ValueError("no overlapping (sample_id, chrom, pos) sites between gnomix and rfmix")

    merged["agree"] = merged["anc_label_gnomix"] == merged["anc_label_rfmix"]

    by_sample = (
        merged.groupby("sample_id", as_index=False)
        .agg(n_sites=("agree", "size"), agreement_rate=("agree", "mean"))
        .sort_values("sample_id")
    )

    overall = pd.DataFrame(
        {
            "sample_id": ["__overall__"],
            "n_sites": [int(merged["agree"].size)],
            "agreement_rate": [float(merged["agree"].mean())],
        }
    )

    summary = pd.concat([by_sample, overall], ignore_index=True)
    return merged, summary


def compute_ancestry_proportions(df):
    if df.empty:
        return pd.DataFrame(columns=["sample_id", "anc_label", "n_sites", "prop_sites"])

    counts = (
        df.groupby(["sample_id", "anc_label"], as_index=False)
        .size()
        .rename(columns={"size": "n_sites"})
    )

    totals = counts.groupby("sample_id", as_index=False)["n_sites"].sum().rename(columns={"n_sites": "n_total"})
    out = counts.merge(totals, on="sample_id", how="left")
    out["prop_sites"] = out["n_sites"] / out["n_total"]
    out = out.drop(columns=["n_total"]).sort_values(["sample_id", "anc_label"])
    return out.reset_index(drop=True)
