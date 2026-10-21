#!/usr/bin/env bash
# Refit models s7 to s20 by running their model-s*.R scripts, three at a time
# (4 chains each, 12 cores). Models are started longest first so the
# lanes finish at about the same time. Output from each script goes to
# models-scripts/<model>-run.log. Run from the project root:
#
#   bash models-scripts/refit-s7-s20.sh
#
# Approximate sampling times at K = 11 (minutes): s14 84, s15 83, s17 66,
# s16 39, s12 37, s20 20, s10 19, s19 19, s13 13, s8 6, s11 6, s18 5, s7 3,
# s9 3.

set -u
cd "$(dirname "$0")/.."

lanes=${LANES:-3}
models=(s14 s15 s17 s16 s12 s20 s10 s19 s13 s8 s11 s18 s7 s9)

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
