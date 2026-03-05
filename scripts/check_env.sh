set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "[check] Repository root: $ROOT_DIR"

required_cmds=(python bash)
optional_cmds=(bcftools plink gnomix rfmix)

missing_required=0
for cmd in "${required_cmds[@]}"; do
  if command -v "$cmd" >/dev/null 2>&1; then
    echo "[ok] command found: $cmd"
  else
    echo "[error] missing required command: $cmd"
    missing_required=1
  fi
done

if [[ -x "/usr/bin/time" ]]; then
  echo "[ok] command found: /usr/bin/time"
else
  echo "[warn] /usr/bin/time not found (benchmark script may fail)"
fi

for cmd in "${optional_cmds[@]}"; do
  if command -v "$cmd" >/dev/null 2>&1; then
    echo "[ok] optional command found: $cmd"
  else
    echo "[warn] optional command missing: $cmd"
  fi
done

python - <<'PY'
import importlib.util

packages = ["pandas", "numpy", "matplotlib", "yaml", "seaborn"]
missing = [p for p in packages if importlib.util.find_spec(p) is None]
if missing:
    print(f"[warn] missing Python packages: {', '.join(missing)}")
else:
    print("[ok] Python packages available: pandas, numpy, matplotlib, yaml, seaborn")
PY

mkdir -p data/sample_lists data/prepared
mkdir -p results/gnomix results/rfmix results/benchmark results/metrics results/figures

echo "[check] done"

if [[ "$missing_required" -ne 0 ]]; then
  exit 1
fi
