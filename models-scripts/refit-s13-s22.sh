#!/usr/bin/env bash
# Refit models s13 and s18 (harmonised metabolic syndrome definition) and the
# moderation models s21 (baseline relative fat mass) and s22 (baseline
# omega-3 index), three at a time (4 chains each, 12 cores). Output from each
# script goes to models-scripts/<model>-run.log. Run from the project root:
#
#   bash models-scripts/refit-s13-s22.sh

set -u
cd "$(dirname "$0")/.."

lanes=${LANES:-3}
models=(s21 s22 s13)

run_one() {
  m=$1
  start=$(date +%s)
  echo "$(date '+%F %T') start $m"
  if Rscript "models-scripts/model-$m.R" > "models-scripts/$m-run.log" 2>&1; then
    status=ok
  else
    status=FAILED
  fi
  echo "$(date '+%F %T') $status $m ($(( ($(date +%s) - start) / 60 )) min)"
}
export -f run_one

printf "%s\n" "${models[@]}" | xargs -P "$lanes" -I{} bash -c 'run_one {}'
