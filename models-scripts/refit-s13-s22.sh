#!/usr/bin/env bash
# Refit models s13 and s18 (harmonised metabolic syndrome definition) and the
# moderation models s19 (own baseline value), s20 (baseline burden score),
# s21 (baseline relative fat mass) and s22 (baseline omega-3 index), three at
# a time (4 chains each, 12 cores). Diagnosis thresholds are passed as data
# (diag_thresholds() in R/data-prep.R). Output from each script goes to
# models-scripts/<model>-run.log. Run from the project root:
#
#   bash models-scripts/refit-s13-s22.sh

set -u
cd "$(dirname "$0")/.."

lanes=${LANES:-3}
models=(s13 s18 s19 s20 s21 s22)

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
