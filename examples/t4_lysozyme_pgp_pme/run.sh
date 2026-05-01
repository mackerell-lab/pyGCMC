#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
RUNS="$HERE/runs"
STEPS="${EXAMPLE_STEPS:-1}"
SEED="${EXAMPLE_SEED:-24681357}"
MODES_CSV="${EXAMPLE_MODES:-pgp_full,pme,pgp_full_pme}"

if [[ -n "${GCMC_CPU:-}" ]]; then
  GCMC_BIN="$GCMC_CPU"
elif command -v gcmc_cpu >/dev/null 2>&1; then
  GCMC_BIN="$(command -v gcmc_cpu)"
else
  GCMC_BIN="$ROOT/build/bin/gcmc_cpu"
fi

if [[ ! -x "$GCMC_BIN" ]]; then
  echo "gcmc_cpu not found or not executable: $GCMC_BIN" >&2
  exit 1
fi

rm -rf "$RUNS"
mkdir -p "$RUNS"
: > "$RUNS/wall_time.tsv"

IFS=',' read -r -a MODES <<< "$MODES_CSV"
for mode in "${MODES[@]}"; do
  mode="$(echo "$mode" | xargs)"
  [[ -z "$mode" ]] && continue
  run_dir="$RUNS/$mode"
  mkdir -p "$run_dir"
  sed \
    -e "s|__ROOT__|$ROOT|g" \
    -e "s|__RUN_DIR__|$run_dir|g" \
    -e "s|__MODE__|$mode|g" \
    -e "s|__SEED__|$SEED|g" \
    -e "s|__STEPS__|$STEPS|g" \
    "$HERE/base.inp.in" > "$run_dir/run.inp"
  cat > "$run_dir/command.txt" <<EOF
$GCMC_BIN --inp $run_dir/run.inp --prefix $run_dir/gcmc --dump-accept $run_dir/acceptance.jsonl --dump-params $run_dir/params.json --print-freq 1000000000 --traj-freq 1000000000 --checkpoint-freq 0 --no-stats
EOF
  start_ns="$(date +%s%N)"
  "$GCMC_BIN" \
    --inp "$run_dir/run.inp" \
    --prefix "$run_dir/gcmc" \
    --dump-accept "$run_dir/acceptance.jsonl" \
    --dump-params "$run_dir/params.json" \
    --print-freq 1000000000 \
    --traj-freq 1000000000 \
    --checkpoint-freq 0 \
    --no-stats \
    > "$run_dir/stdout.log" \
    2> "$run_dir/stderr.log"
  end_ns="$(date +%s%N)"
  python3 - "$mode" "$start_ns" "$end_ns" "$RUNS/wall_time.tsv" <<'PY'
import sys
mode, start, end, path = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
elapsed = (end - start) / 1_000_000_000.0
with open(path, "a", encoding="utf-8") as handle:
    handle.write(f"{mode}\t{elapsed:.6f}\n")
PY
done

python3 "$ROOT/examples/common/analyze_pgp_pme.py" \
  "$RUNS" \
  --system "T4 lysozyme, 37,871 atoms, 64^3 PME mesh" \
  --json "$RUNS/summary.json" \
  --log "$RUNS/paper_summary.log"
