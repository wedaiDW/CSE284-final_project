from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.utils import load_yaml, read_sample_list


def parse_args():
    parser = argparse.ArgumentParser(description="Create reference/test sample lists.")
    parser.add_argument("--config", default="configs/dataset.yaml", help="Path to dataset yaml")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing split files")
    return parser.parse_args()


def require_columns(df, required):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"panel file missing columns: {missing}. expected at least {required}")


def sample_by_pop(
    df,
    pop_col,
    sample_col,
    pops,
    random_seed,
    n_per_pop = None,
):
    chunks = []
    for pop in pops:
        sub = df[df[pop_col] == pop].copy()
        if sub.empty:
            print(f"[warn] no samples found for pop={pop}")
            continue
        if n_per_pop is None or n_per_pop >= len(sub):
            picked = sub
        else:
            picked = sub.sample(n=n_per_pop, random_state=random_seed)
        chunks.append(picked)

    if not chunks:
        return pd.DataFrame(columns=df.columns)

    out = pd.concat(chunks, axis=0, ignore_index=True)
    out = out.drop_duplicates(subset=[sample_col])
    return out


def main():
    args = parse_args()
    cfg = load_yaml(args.config)

    panel_path = Path(cfg["panel_tsv"])
    if not panel_path.exists():
        raise SystemExit(f"[error] panel file not found: {panel_path}")

    sample_col = cfg.get("panel_sample_col", "sample")
    pop_col = cfg.get("panel_pop_col", "pop")

    df = pd.read_csv(panel_path, sep="\t")
    require_columns(df, [sample_col, pop_col])

    reference_pops = list(cfg.get("reference_pops", []))
    test_pops = list(cfg.get("admixed_test_pops", []))
    n_ref_per_pop = cfg.get("n_ref_per_pop")
    n_test = cfg.get("n_test")
    random_seed = int(cfg.get("random_seed", 42))

    ref_df = sample_by_pop(
        df=df,
        pop_col=pop_col,
        sample_col=sample_col,
        pops=reference_pops,
        n_per_pop=n_ref_per_pop,
        random_seed=random_seed,
    )

    test_pool = df[df[pop_col].isin(test_pops)].copy()
    test_pool = test_pool[~test_pool[sample_col].isin(ref_df[sample_col])].copy()
    if test_pool.empty:
        raise SystemExit("[error] no test samples remain after filtering and overlap removal")

    if n_test is None or n_test >= len(test_pool):
        test_df = test_pool
    else:
        test_df = test_pool.sample(n=n_test, random_state=random_seed)

    ref_path = Path(cfg.get("reference_list", "./data/sample_lists/reference.txt"))
    test_path = Path(cfg.get("test_list", "./data/sample_lists/test.txt"))
    out_dir = ref_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    if (ref_path.exists() or test_path.exists()) and not args.overwrite:
        raise SystemExit(
            "[error] split files already exist. use --overwrite to replace. "
            f"existing: {ref_path}, {test_path}"
        )

    ref_df[[sample_col]].to_csv(ref_path, sep="\t", index=False, header=False)
    test_df[[sample_col]].to_csv(test_path, sep="\t", index=False, header=False)

    ref_detail = out_dir / "reference_with_pop.tsv"
    test_detail = out_dir / "test_with_pop.tsv"
    ref_df[[sample_col, pop_col]].rename(columns={sample_col: "sample", pop_col: "pop"}).to_csv(
        ref_detail,
        sep="\t",
        index=False,
    )
    test_df[[sample_col, pop_col]].rename(columns={sample_col: "sample", pop_col: "pop"}).to_csv(
        test_detail,
        sep="\t",
        index=False,
    )

    summary = pd.DataFrame(
        {
            "group": ["reference", "test"],
            "n_samples": [len(ref_df), len(test_df)],
            "pops": [",".join(sorted(ref_df[pop_col].unique())), ",".join(sorted(test_df[pop_col].unique()))],
        }
    )
    summary_path = out_dir / "split_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    _ = read_sample_list(ref_path)
    _ = read_sample_list(test_path)

    print(f"reference list: {ref_path} ({len(ref_df)} samples)")
    print(f"test list: {test_path} ({len(test_df)} samples)")
    print(f"summary: {summary_path}")


if __name__ == "__main__":
    main()
