from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import yaml


def load_yaml(path):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"config not found: {p}")
    with p.open("r") as f:
        return yaml.safe_load(f) or {}


def ensure_dir(path):
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def read_sample_list(path):
    p = Path(path)
    if not p.exists():
        return []
    return [line.strip() for line in p.read_text().splitlines() if line.strip()]


def normalize_chrom_value(value):
    text = str(value)
    text = text.replace("chr", "").replace("CHR", "")
    return text


def apply_label_map(series, label_map):
    if not label_map:
        return series.astype(str)

    normalized_map = {str(k): str(v) for k, v in label_map.items()}
    return series.astype(str).map(lambda x: normalized_map.get(x, x))


def standardize_per_site(df):
    col_aliases = {
        "sample": "sample_id",
        "id": "sample_id",
        "ind": "sample_id",
        "individual": "sample_id",
        "chr": "chrom",
        "chromosome": "chrom",
        "position": "pos",
        "bp": "pos",
        "label": "anc_label",
        "ancestry": "anc_label",
    }

    renamed = {}
    for c in df.columns:
        key = c.strip()
        renamed[c] = col_aliases.get(key, key)

    df = df.rename(columns=renamed)

    required = ["sample_id", "chrom", "pos", "anc_label"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"per-site table missing required columns: {missing}")

    out = df.copy()
    out["sample_id"] = out["sample_id"].astype(str)
    out["chrom"] = out["chrom"].map(normalize_chrom_value)
    out["pos"] = pd.to_numeric(out["pos"], errors="coerce")
    out = out.dropna(subset=["pos"])
    out["pos"] = out["pos"].astype(int)
    out["anc_label"] = out["anc_label"].astype(str)
    out = out.sort_values(["sample_id", "chrom", "pos"]).reset_index(drop=True)
    return out
