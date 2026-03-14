#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON:-python}"
CONFIG="${CONFIG:-configs/dataset.yaml}"

"$PYTHON_BIN" scripts/postprocess_metrics.py --config "$CONFIG"
"$PYTHON_BIN" scripts/make_figures.py --config "$CONFIG"
