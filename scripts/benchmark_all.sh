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
TIME_ARGS="-v"
if [[ "$(uname -s)" == "Darwin" ]]; then
  TIME_BIN=""
  TIME_ARGS=""
fi

if [[ -n "$TIME_BIN" ]]; then
  echo "[benchmark] using timer: $TIME_BIN"
else
  echo "[benchmark] timing disabled on Darwin"
fi

run_timed() {
  local outfile="$1"
  shift

  if [[ -z "$TIME_BIN" ]]; then
    echo "[warn] timing skipped on Darwin to avoid /usr/bin/time permission issues" >&2
    printf "timer_skipped=1\n" >"$outfile"
    "$@"
    return 0
  fi

  if "$TIME_BIN" "$TIME_ARGS" -o "$outfile" "$@"; then
    return 0
  fi

  echo "[warn] timer failed for: $*" >&2
  echo "[warn] falling back to untimed execution" >&2
  printf "timer_failed=1\n" >"$outfile"
  "$@"
}

run_timed results/benchmark/prepare.time.txt \
  python scripts/prepare_inputs.py --config "$CONFIG"

run_timed results/benchmark/gnomix.time.txt \
  python scripts/run_gnomix.py --config "$CONFIG" --tools "$TOOLS_CONFIG"

run_timed results/benchmark/rfmix.time.txt \
  python scripts/run_rfmix.py --config "$CONFIG" --tools "$TOOLS_CONFIG"

run_timed results/benchmark/postprocess.time.txt \
  python scripts/postprocess_metrics.py --config "$CONFIG"

echo "[ok] benchmark logs saved to results/benchmark/"
