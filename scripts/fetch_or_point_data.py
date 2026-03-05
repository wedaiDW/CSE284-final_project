from __future__ import annotations

import argparse
from pathlib import Path

import yaml

ROOT_DIR = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate configured panel/VCF paths.")
    parser.add_argument("--config", default="configs/dataset.yaml", help="Path to dataset yaml")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        raise SystemExit(f"[error] config not found: {config_path}")

    cfg = yaml.safe_load(config_path.read_text()) or {}

    data_root = Path(cfg.get("data_root", "./data/raw"))
    panel_tsv = Path(cfg.get("panel_tsv", data_root / "panel.tsv"))
    vcf_pattern = cfg.get("vcf_pattern", "./data/raw/chr{chrom}.vcf.gz")
    chromosomes = [str(c) for c in cfg.get("chromosomes", ["22"])]

    data_root.mkdir(parents=True, exist_ok=True)

    print(f"data_root: {data_root}")
    print(f"panel_tsv: {panel_tsv}")
    print(f"vcf_pattern: {vcf_pattern}")
    print(f"chromosomes: {', '.join(chromosomes)}")

    missing = []
    if not panel_tsv.exists():
        missing.append(str(panel_tsv))

    for chrom in chromosomes:
        vcf = Path(vcf_pattern.format(chrom=chrom))
        if not vcf.exists():
            missing.append(str(vcf))

    if missing:
        print("\nSome expected files are missing:")
        for path in missing:
            print(f"  - {path}")
        print("\nPlace data at the data root.")
    else:
        print("Sanity check OK")


if __name__ == "__main__":
    main()
