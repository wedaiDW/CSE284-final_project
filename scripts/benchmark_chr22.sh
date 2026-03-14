#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON:-python}"
CONFIG="${CONFIG:-configs/dataset.yaml}"
TOOLS_CONFIG="${TOOLS_CONFIG:-configs/tools.yaml}"

"$PYTHON_BIN" scripts/prepare_inputs.py --config "$CONFIG"
bash scripts/benchmark_all.sh --config "$CONFIG" --tools "$TOOLS_CONFIG"
