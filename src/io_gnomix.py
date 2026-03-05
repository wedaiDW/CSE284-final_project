from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.utils import standardize_per_site


def _load_tsv_files(files):
    frames = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t", comment="#")
        except Exception as e:
            raise RuntimeError(f"failed reading {f}: {e}") from e
        if df.empty:
            continue
        df["source_file"] = str(f)
        frames.append(df)

    if not frames:
        raise ValueError("no non-empty gnomix per-site files were found")

    return standardize_per_site(pd.concat(frames, ignore_index=True))


def load_gnomix_per_site(root):
    root = Path(root)
    if not root.exists():
        raise FileNotFoundError(f"gnomix output root not found: {root}")

    files = sorted(root.rglob("per_site.tsv"))
    if not files:
        files = sorted(root.rglob("*.tsv"))

    if not files:
        raise FileNotFoundError(
            f"no per-site TSV found under {root}. expected results/gnomix/chr*/per_site.tsv"
        )

    return _load_tsv_files(files)
