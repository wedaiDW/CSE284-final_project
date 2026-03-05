#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import yaml

ROOT_DIR = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare subset benchmark VCFs.")
    parser.add_argument("--config", default="configs/dataset.yaml", help="Path to dataset yaml")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text()) or {}

    chromosomes = [str(c) for c in cfg.get("chromosomes", ["22"])]
    vcf_pattern = cfg.get("vcf_pattern", "./data/raw/chr{chrom}.vcf.gz")
    prepared_root = Path(cfg.get("prepared_root", "./data/prepared"))
    reference_list = Path(cfg.get("reference_list", "./data/sample_lists/reference.txt"))
    test_list = Path(cfg.get("test_list", "./data/sample_lists/test.txt"))

    if not reference_list.exists() or not test_list.exists():
        raise SystemExit(
            "[error] sample list files are missing. run scripts/make_splits.py first."
        )

    prepared_root.mkdir(parents=True, exist_ok=True)
    manifest_path = prepared_root / "manifest.tsv"

    with reference_list.open() as f:
        ref_samples = [line.strip() for line in f if line.strip()]
    with test_list.open() as f:
        test_samples = [line.strip() for line in f if line.strip()]
    all_samples = sorted(set(ref_samples + test_samples))

    if not all_samples:
        raise SystemExit("[error] no samples found in reference/test lists")

    bcftools = shutil.which("bcftools")
    rows: list[str] = ["chrom\tprepared_vcf\treference_list\ttest_list"]

    for chrom in chromosomes:
        source_vcf = Path(vcf_pattern.format(chrom=chrom))
        chrom_dir = prepared_root / f"chr{chrom}"
        chrom_dir.mkdir(parents=True, exist_ok=True)

        sample_file = chrom_dir / "all_samples.txt"
        sample_file.write_text("\n".join(all_samples) + "\n")

        out_vcf = chrom_dir / f"benchmark_input.chr{chrom}.vcf.gz"

        if source_vcf.exists() and bcftools:
            cmd = [
                bcftools,
                "view",
                "-r",
                str(chrom),
                "-S",
                str(sample_file),
                "-Oz",
                "-o",
                str(out_vcf),
                str(source_vcf),
            ]
            print("[run]", " ".join(cmd))
            subprocess.run(cmd, check=True)
            subprocess.run([bcftools, "index", "-f", str(out_vcf)], check=True)
        elif source_vcf.exists() and not bcftools:
            out_vcf.write_text(
                "# bcftools missing; placeholder created. install bcftools to prepare real subset VCF.\n"
            )
            print("[warn] bcftools not found. wrote placeholder VCF:", out_vcf)
        else:
            out_vcf.write_text(
                "# source VCF missing; placeholder created. update configs/dataset.yaml::vcf_pattern\n"
            )
            print("[warn] source VCF missing for chrom", chrom, "expected", source_vcf)

        rows.append(f"{chrom}\t{out_vcf}\t{reference_list}\t{test_list}")

    manifest_path.write_text("\n".join(rows) + "\n")
    print(f"[ok] wrote manifest: {manifest_path}")


if __name__ == "__main__":
    main()
