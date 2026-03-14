#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON:-python}"
CONFIG="${CONFIG:-configs/dataset.yaml}"
TOOLS_CONFIG="${TOOLS_CONFIG:-configs/tools.yaml}"
RESOLVED_CONFIG="${RESOLVED_CONFIG:-configs/dataset.runtime.chr22.yaml}"

"$PYTHON_BIN" scripts/run_pipeline.py \
  --config "$CONFIG" \
  --tools "$TOOLS_CONFIG" \
  --chromosomes 22 \
  --resolved-config "$RESOLVED_CONFIG" \
  "$@"
