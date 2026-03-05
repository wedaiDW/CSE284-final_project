from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def _save(fig, path):
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_agreement_distribution(site_agreement, out_path):
    data = site_agreement[site_agreement["sample_id"] != "__overall__"].copy()

    fig, ax = plt.subplots(figsize=(8, 4.5))
    if "agreement_rate" in data.columns and not data.empty:
        sns.histplot(data["agreement_rate"], bins=20, ax=ax, color="#1f77b4")
        ax.set_xlabel("Per-sample agreement rate")
        ax.set_ylabel("Count")
        ax.set_title("Gnomix vs RFMix agreement distribution")
    else:
        ax.text(0.5, 0.5, "No agreement data", ha="center", va="center")
        ax.set_axis_off()

    _save(fig, out_path)


def plot_tract_length_distribution(tract_segments, out_path):
    fig, ax = plt.subplots(figsize=(8, 4.5))

    if not tract_segments.empty and "tract_length_bp" in tract_segments.columns:
        plot_df = tract_segments.copy()
        plot_df["tract_length_bp"] = pd.to_numeric(plot_df["tract_length_bp"], errors="coerce")
        plot_df = plot_df.dropna(subset=["tract_length_bp"])
        if "tool" in plot_df.columns:
            sns.histplot(
                data=plot_df,
                x="tract_length_bp",
                hue="tool",
                bins=30,
                element="step",
                stat="density",
                common_norm=False,
                ax=ax,
            )
        else:
            sns.histplot(plot_df["tract_length_bp"], bins=30, ax=ax)
        ax.set_xlabel("Tract length (bp)")
        ax.set_ylabel("Density")
        ax.set_title("Tract length distribution")
    else:
        ax.text(0.5, 0.5, "No tract data", ha="center", va="center")
        ax.set_axis_off()

    _save(fig, out_path)


def plot_ancestry_stacked_bar(ancestry_props, out_path):
    fig, ax = plt.subplots(figsize=(10, 5))

    if ancestry_props.empty:
        ax.text(0.5, 0.5, "No ancestry proportion data", ha="center", va="center")
        ax.set_axis_off()
        _save(fig, out_path)
        return

    plot_df = ancestry_props.copy()
    if "tool" in plot_df.columns:
        plot_df["group"] = plot_df["tool"] + "|" + plot_df["sample_id"]
    else:
        plot_df["group"] = plot_df["sample_id"]

    groups = plot_df["group"].drop_duplicates().tolist()
    if len(groups) > 30:
        keep = set(groups[:30])
        plot_df = plot_df[plot_df["group"].isin(keep)]

    pivot = (
        plot_df.pivot_table(index="group", columns="anc_label", values="prop_sites", aggfunc="mean", fill_value=0)
        .sort_index()
    )

    bottom = None
    for anc in pivot.columns:
        values = pivot[anc].values
        ax.bar(pivot.index, values, label=anc, bottom=bottom)
        if bottom is None:
            bottom = values
        else:
            bottom = bottom + values

    ax.set_ylabel("Proportion")
    ax.set_xlabel("Sample (tool|sample_id)")
    ax.set_title("Ancestry proportion summary")
    ax.tick_params(axis="x", labelrotation=90)
    ax.legend(title="Ancestry", bbox_to_anchor=(1.01, 1), loc="upper left")

    _save(fig, out_path)


def plot_label_confusion_heatmap(label_confusion_rates, out_path):
    fig, ax = plt.subplots(figsize=(6, 5))

    if label_confusion_rates.empty:
        ax.text(0.5, 0.5, "No confusion matrix data", ha="center", va="center")
        ax.set_axis_off()
        _save(fig, out_path)
        return

    heat_df = label_confusion_rates.copy() * 100.0
    sns.heatmap(
        heat_df,
        annot=True,
        fmt=".1f",
        cmap="YlOrRd",
        cbar_kws={"label": "Percent of all shared sites (%)"},
        ax=ax,
    )
    ax.set_title("Gnomix vs RFMix label confusion")
    ax.set_xlabel("RFMix label")
    ax.set_ylabel("Gnomix label")

    _save(fig, out_path)


def plot_tool_comparison_overview(
    tool_summary,
    comparison_dashboard,
    out_path,
):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    if tool_summary.empty:
        ax1.text(0.5, 0.5, "No tool summary data", ha="center", va="center")
        ax1.set_axis_off()
        ax2.set_axis_off()
        _save(fig, out_path)
        return

    plot_df = tool_summary.copy()
    x = plot_df["tool"].astype(str)

    sns.barplot(
        data=plot_df,
        x="tool",
        y="mean_tract_length_bp",
        hue="tool",
        dodge=False,
        legend=False,
        ax=ax1,
        palette="Set2",
    )
    ax1.set_title("Mean tract length by tool")
    ax1.set_xlabel("Tool")
    ax1.set_ylabel("Mean tract length (bp)")

    anc_cols = sorted([c for c in plot_df.columns if c.startswith("mean_prop_")])
    if anc_cols:
        bottom = pd.Series([0.0] * len(plot_df))
        for col in anc_cols:
            anc = col.replace("mean_prop_", "")
            values = pd.to_numeric(plot_df[col], errors="coerce").fillna(0.0)
            ax2.bar(x, values, bottom=bottom, label=anc)
            bottom = bottom + values
        ax2.set_ylim(0, 1.0)
        ax2.set_title("Mean ancestry composition by tool")
        ax2.set_xlabel("Tool")
        ax2.set_ylabel("Proportion")
        ax2.legend(title="Ancestry", loc="upper right")
    else:
        ax2.text(0.5, 0.5, "No ancestry proportion data", ha="center", va="center")
        ax2.set_axis_off()

    if not comparison_dashboard.empty and "overall_agreement_rate" in comparison_dashboard.columns:
        overall = float(comparison_dashboard["overall_agreement_rate"].iloc[0]) * 100.0
        fig.suptitle(f"Head-to-head overview (overall agreement: {overall:.2f}%)")
    else:
        fig.suptitle("Head-to-head overview")

    _save(fig, out_path)
