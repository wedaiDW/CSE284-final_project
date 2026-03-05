#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml

ROOT_DIR = Path(__file__).resolve().parents[1]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Unified LAI benchmark entrypoint with optional path/chromosome overrides."
    )
    parser.add_argument("--config", default="configs/dataset.yaml", help="Base dataset config")
    parser.add_argument("--tools", default="configs/tools.yaml", help="Tools config")
    parser.add_argument(
        "--resolved-config",
        default="configs/dataset.runtime.yaml",
        help="Path to write the resolved runtime config",
    )
    parser.add_argument("--data-root", help="Override data_root, e.g. /abs/path/to/data")
    parser.add_argument("--panel-tsv", help="Override panel_tsv file path")
    parser.add_argument(
        "--vcf-pattern",
        help="Override vcf_pattern, e.g. /abs/path/1000G_chr{chrom}_pruned.vcf.gz",
    )
    parser.add_argument(
        "--chromosomes",
        help="Override chromosomes, comma-separated (e.g. 22 or 21,22)",
    )
    parser.add_argument(
        "--overwrite-splits",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Whether to pass --overwrite to scripts/make_splits.py",
    )
    parser.add_argument(
        "--env-check",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Whether to run scripts/check_env.sh before pipeline",
    )
    parser.add_argument(
        "--make-figures",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Whether to run scripts/make_figures.py",
    )
    parser.add_argument(
        "--write-config-only",
        action="store_true",
        help="Only write resolved config and exit without running pipeline",
    )
    return parser.parse_args()


def resolve_path(path_like: str, *, allow_braces: bool = False) -> str:
    text = str(path_like).strip()
    if allow_braces:
        p = Path(text.replace("{chrom}", "__CHROM__"))
        if p.is_absolute():
            return text
        return str((ROOT_DIR / text).resolve())

    p = Path(text).expanduser()
    if p.is_absolute():
        return str(p)
    return str((ROOT_DIR / p).resolve())


def parse_chromosomes(value: str) -> list[str]:
    items = [x.strip() for x in value.replace(" ", ",").split(",") if x.strip()]
    if not items:
        raise ValueError("--chromosomes is empty")
    out = []
    for chrom in items:
        c = chrom.lower().replace("chr", "")
        if not c:
            continue
        out.append(c)
    if not out:
        raise ValueError("no valid chromosome values were parsed")
    return out


def load_config(path):
    p = Path(resolve_path(path))
    if not p.exists():
        raise FileNotFoundError(f"config not found: {p}")
    with p.open("r") as f:
        return yaml.safe_load(f) or {}


def guess_panel_path(data_root):
    candidates = [data_root / "meta" / "panel.tsv", data_root / "panel.tsv"]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def run_cmd(cmd):
    printable = " ".join(shlex.quote(x) for x in cmd)
    print(f"[run] {printable}")
    subprocess.run(cmd, cwd=ROOT_DIR, check=True)


def main():
    args = parse_args()

    cfg = load_config(args.config)
    tools = resolve_path(args.tools)
    resolved_cfg_path = Path(resolve_path(args.resolved_config))
    resolved_cfg_path.parent.mkdir(parents=True, exist_ok=True)

    if args.data_root:
        data_root = Path(resolve_path(args.data_root))
        cfg["data_root"] = str(data_root)
        if not args.panel_tsv:
            cfg["panel_tsv"] = str(guess_panel_path(data_root))
        if not args.vcf_pattern:
            cfg["vcf_pattern"] = str(data_root / "1000G_chr{chrom}_pruned.vcf.gz")

    if args.panel_tsv:
        cfg["panel_tsv"] = resolve_path(args.panel_tsv)
    if args.vcf_pattern:
        cfg["vcf_pattern"] = resolve_path(args.vcf_pattern, allow_braces=True)
    if args.chromosomes:
        cfg["chromosomes"] = parse_chromosomes(args.chromosomes)

    resolved_cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

    print(f"[info] resolved config: {resolved_cfg_path}")
    print(f"[info] data_root: {cfg.get('data_root')}")
    print(f"[info] panel_tsv: {cfg.get('panel_tsv')}")
    print(f"[info] vcf_pattern: {cfg.get('vcf_pattern')}")
    print(f"[info] chromosomes: {cfg.get('chromosomes')}")

    if args.write_config_only:
        print("[ok] config written only")
        return

    if args.env_check:
        run_cmd(["bash", "scripts/check_env.sh"])

    run_cmd([sys.executable, "scripts/fetch_or_point_data.py", "--config", str(resolved_cfg_path)])

    split_cmd = [sys.executable, "scripts/make_splits.py", "--config", str(resolved_cfg_path)]
    if args.overwrite_splits:
        split_cmd.append("--overwrite")
    run_cmd(split_cmd)

    run_cmd(
        [
            "bash",
            "scripts/benchmark_all.sh",
            "--config",
            str(resolved_cfg_path),
            "--tools",
            str(tools),
        ]
    )

    if args.make_figures:
        run_cmd([sys.executable, "scripts/make_figures.py", "--config", str(resolved_cfg_path)])

    print("Run complete. Check results directory.")


if __name__ == "__main__":
    main()
