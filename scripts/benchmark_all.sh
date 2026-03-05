#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

CONFIG="configs/dataset.yaml"
TOOLS_CONFIG="configs/tools.yaml"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG="$2"
      shift 2
      ;;
    --tools)
      TOOLS_CONFIG="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

mkdir -p results/benchmark
TIME_BIN="/usr/bin/time"
if [[ ! -x "$TIME_BIN" ]]; then
  TIME_BIN="time"
fi

echo "[benchmark] using timer: $TIME_BIN"
TIME_ARGS="-v"
if [[ "$(uname -s)" == "Darwin" ]]; then
  TIME_ARGS="-l"
fi

"$TIME_BIN" "$TIME_ARGS" -o results/benchmark/prepare.time.txt \
  python scripts/prepare_inputs.py --config "$CONFIG"

"$TIME_BIN" "$TIME_ARGS" -o results/benchmark/gnomix.time.txt \
  python scripts/run_gnomix.py --config "$CONFIG" --tools "$TOOLS_CONFIG"

"$TIME_BIN" "$TIME_ARGS" -o results/benchmark/rfmix.time.txt \
  python scripts/run_rfmix.py --config "$CONFIG" --tools "$TOOLS_CONFIG"

"$TIME_BIN" "$TIME_ARGS" -o results/benchmark/postprocess.time.txt \
  python scripts/postprocess_metrics.py --config "$CONFIG"

echo "[ok] benchmark logs saved to results/benchmark/"
