#!/usr/bin/env python3
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import shutil
import subprocess
from pathlib import Path

import yaml

ROOT_DIR = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run RFMix benchmark (real or mock fallback).")
    parser.add_argument("--config", default="configs/dataset.yaml", help="Path to dataset yaml")
    parser.add_argument("--tools", default="configs/tools.yaml", help="Path to tools yaml")
    return parser.parse_args()


def read_positions(vcf_path, max_sites=5000):
    if not vcf_path.exists():
        return []

    opener = gzip.open if vcf_path.suffix == ".gz" else open
    positions = []
    try:
        with opener(vcf_path, "rt") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                positions.append(int(parts[1]))
                if max_sites is not None and len(positions) >= max_sites:
                    break
    except Exception:
        return []

    return positions


def mock_per_site(out_file, samples, chrom, positions):
    labels = ["AFR", "EUR", "EAS"]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w") as f:
        f.write("sample_id\tchrom\tpos\tanc_label\n")
        for sample in samples:
            for pos in positions:
                key = f"rfmix::{sample}::{chrom}::{pos}".encode()
                idx = int(hashlib.sha1(key).hexdigest(), 16) % len(labels)
                f.write(f"{sample}\t{chrom}\t{pos}\t{labels[idx]}\n")


def run_cmd(cmd):
    subprocess.run(cmd, check=True)


def subset_vcf(vcf_in, samples_file, vcf_out, bcftools_bin):
    run_cmd(
        [
            bcftools_bin,
            "view",
            "-S",
            str(samples_file),
            "-Oz",
            "-o",
            str(vcf_out),
            str(vcf_in),
        ]
    )
    run_cmd([bcftools_bin, "index", "-f", str(vcf_out)])


def load_panel_map(panel_tsv):
    with panel_tsv.open("r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "sample" not in reader.fieldnames or "pop" not in reader.fieldnames:
            raise ValueError("panel_tsv must contain columns: sample, pop")
        return {row["sample"]: row["pop"] for row in reader if row.get("sample") and row.get("pop")}


def write_sample_map(reference_samples, panel_map, out_path):
    missing = [s for s in reference_samples if s not in panel_map]
    if missing:
        raise ValueError(f"sample_map missing {len(missing)} reference samples in panel (example: {missing[:5]})")

    with out_path.open("w") as f:
        for sample in reference_samples:
            f.write(f"{sample}\t{panel_map[sample]}\n")


def write_synthetic_genetic_map(chrom, positions, out_path):
    if not positions:
        raise ValueError("cannot build genetic map from empty position list")
    start = positions[0]
    with out_path.open("w") as f:
        for pos in positions:
            cm = (pos - start) / 1_000_000.0
            f.write(f"{chrom}\t{pos}\t{cm:.6f}\n")


def parse_code_to_label(line):
    after_colon = line.split(":", 1)[-1].replace("\t", " ").strip()
    code_to_label: dict[int, str] = {}
    for token in after_colon.split():
        if "=" not in token:
            continue
        label, code = token.split("=", 1)
        code_to_label[int(code)] = label
    if not code_to_label:
        raise ValueError("failed parsing subpopulation codes from RFMix MSP header")
    return code_to_label


def sample_haplotype_columns(header_cols):
    pairs = {}
    for idx, col in enumerate(header_cols[6:], start=6):
        sample = col.rsplit(".", 1)[0]
        pairs.setdefault(sample, []).append(idx)

    out = {}
    for sample, idxs in pairs.items():
        if len(idxs) >= 2:
            out[sample] = (idxs[0], idxs[1])
    return out


def write_per_site_from_msp(
    msp_path,
    out_file,
    samples,
    positions,
    chrom,
):
    if not positions:
        raise ValueError("cannot render per_site.tsv from empty position list")

    with msp_path.open("r") as f:
        first = f.readline().strip()
        second = f.readline().strip()
        if not first.startswith("#Subpopulation order/codes:"):
            raise ValueError("unexpected MSP header: missing subpopulation code line")
        if not second.startswith("#chm"):
            raise ValueError("unexpected MSP header: missing column header line")

        code_to_label = parse_code_to_label(first)
        header_cols = second.lstrip("#").split("\t")
        hap_pairs = sample_haplotype_columns(header_cols)

        use_samples = [s for s in samples if s in hap_pairs]
        if not use_samples:
            raise ValueError("none of test samples were found in MSP output")

        out_file.parent.mkdir(parents=True, exist_ok=True)
        with out_file.open("w") as w:
            w.write("sample_id\tchrom\tpos\tanc_label\n")
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                seg_start = int(parts[1])
                seg_end = int(parts[2])

                lo = bisect.bisect_left(positions, seg_start)
                hi = bisect.bisect_right(positions, seg_end)
                seg_positions = positions[lo:hi]
                if not seg_positions:
                    continue

                for sample in use_samples:
                    i0, i1 = hap_pairs[sample]
                    c0 = int(float(parts[i0]))
                    c1 = int(float(parts[i1]))
                    chosen_code = c0 if c0 == c1 else c0
                    anc_label = code_to_label.get(chosen_code, str(chosen_code))
                    for pos in seg_positions:
                        w.write(f"{sample}\t{chrom}\t{pos}\t{anc_label}\n")


def main():
    args = parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text()) or {}
    tools = yaml.safe_load(Path(args.tools).read_text()) or {}

    rfmix_cfg = tools.get("rfmix", {})
    threads = int(tools.get("threads", 4))
    chromosomes = [str(c) for c in cfg.get("chromosomes", ["22"])]
    prepared_root = Path(cfg.get("prepared_root", "./data/prepared"))
    reference_list = Path(cfg.get("reference_list", "./data/sample_lists/reference.txt"))
    test_list = Path(cfg.get("test_list", "./data/sample_lists/test.txt"))
    panel_tsv = Path(cfg.get("panel_tsv", "./data/raw/panel.tsv"))
    genetic_map_pattern = cfg.get("rfmix_genetic_map_pattern")

    if not test_list.exists():
        raise SystemExit("[error] test list missing. run scripts/make_splits.py first")
    if not reference_list.exists():
        raise SystemExit("[error] reference list missing. run scripts/make_splits.py first")

    samples = [s.strip() for s in test_list.read_text().splitlines() if s.strip()]
    reference_samples = [s.strip() for s in reference_list.read_text().splitlines() if s.strip()]
    if not samples:
        raise SystemExit("[error] no test samples found in test list")
    if not reference_samples:
        raise SystemExit("[error] no reference samples found in reference list")

    rfmix_bin = shutil.which("rfmix")
    bcftools_bin = shutil.which("bcftools")
    enabled = bool(rfmix_cfg.get("enabled", True))
    cmd_template = rfmix_cfg.get("cmd_template", "")
    panel_map = load_panel_map(panel_tsv) if panel_tsv.exists() else {}

    for chrom in chromosomes:
        out_dir = Path("results/rfmix") / f"chr{chrom}"
        out_dir.mkdir(parents=True, exist_ok=True)

        prepared_vcf = prepared_root / f"chr{chrom}" / f"benchmark_input.chr{chrom}.vcf.gz"
        query_haps = out_dir / f"query.chr{chrom}.vcf.gz"
        ref_haps = out_dir / f"ref.chr{chrom}.vcf.gz"
        sample_map = out_dir / f"sample_map.chr{chrom}.tsv"
        synthetic_map = out_dir / f"genetic_map.chr{chrom}.synthetic.tsv"
        out_prefix = out_dir / "rfmix"
        msp_path = Path(str(out_prefix) + ".msp.tsv")
        log_path = out_dir / "run.log"
        per_site_path = out_dir / "per_site.tsv"

        positions = read_positions(prepared_vcf, max_sites=5000)
        if not positions:
            positions = [10_000 + i * 100 for i in range(1000)]

        ran_real = False
        if enabled and rfmix_bin and cmd_template and bcftools_bin and panel_map and prepared_vcf.exists():
            try:
                subset_vcf(prepared_vcf, reference_list, ref_haps, bcftools_bin)
                subset_vcf(prepared_vcf, test_list, query_haps, bcftools_bin)
                write_sample_map(reference_samples, panel_map, sample_map)

                if genetic_map_pattern:
                    genetic_map = Path(str(genetic_map_pattern).format(chrom=chrom))
                    if not genetic_map.exists():
                        all_positions = read_positions(ref_haps, max_sites=None)
                        write_synthetic_genetic_map(chrom, all_positions, synthetic_map)
                        genetic_map = synthetic_map
                else:
                    all_positions = read_positions(ref_haps, max_sites=None)
                    write_synthetic_genetic_map(chrom, all_positions, synthetic_map)
                    genetic_map = synthetic_map

                format_dict = {
                    "query_haps": query_haps,
                    "ref_haps": ref_haps,
                    "sample_map": sample_map,
                    "genetic_map": genetic_map,
                    "out_prefix": out_prefix,
                    "threads": threads,
                    "chrom": chrom,
                }
                cmd = cmd_template.format(**format_dict)
                if rfmix_cfg.get("extra_args"):
                    cmd = f"{cmd} {rfmix_cfg['extra_args']}"

                print("[run]", cmd)
                with log_path.open("w") as logf:
                    proc = subprocess.run(cmd, shell=True, stdout=logf, stderr=subprocess.STDOUT)

                if proc.returncode == 0 and msp_path.exists():
                    write_per_site_from_msp(msp_path, per_site_path, samples=samples, positions=positions, chrom=chrom)
                    ran_real = True
                else:
                    print(f"[warn] rfmix command failed on chr{chrom}. fallback to mock per-site output")
            except Exception as e:
                print(f"[warn] real rfmix setup failed on chr{chrom}: {e}. fallback to mock per-site output")
                log_path.write_text(f"real rfmix setup failed: {e}\n")

        if not ran_real:
            if not log_path.exists():
                log_path.write_text(
                    "rfmix executable/template unavailable or failed. mock per-site output generated.\n"
                )
            mock_per_site(per_site_path, samples=samples, chrom=chrom, positions=positions)

        mode = "real" if ran_real else "mock"
        print(f"result wrote {per_site_path} ({mode})")


if __name__ == "__main__":
    main()
