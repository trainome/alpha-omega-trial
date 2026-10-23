#!/usr/bin/env bash
# Refit models s14 to s17 with D = 6 latent factors by running their
# model-s*-D6.R scripts, three at a time (4 chains each, 12 cores). Models are
# started longest first. Output from each script goes to
# models-scripts/<model>-D6-run.log. Run from the project root:
#
#   bash models-scripts/refit-s14-s17-d6.sh
#
# Approximate sampling times at D = 4, K = 11 (minutes): s14 84, s15 83,
# s17 66, s16 39. Expect longer run times with D = 6.

set -u
cd "$(dirname "$0")/.."

lanes=${LANES:-3}
models=(s14 s15 s17 s16)

run_one() {
  m=$1
  start=$(date +%s)
  echo "$(date '+%F %T') start $m-D6"
  if Rscript "models-scripts/model-$m-D6.R" > "models-scripts/$m-D6-run.log" 2>&1; then
    status=ok
  else
    status=FAILED
  fi
  echo "$(date '+%F %T') $status $m-D6 ($(( ($(date +%s) - start) / 60 )) min)"
}
export -f run_one

printf "%s\n" "${models[@]}" | xargs -P "$lanes" -I{} bash -c 'run_one {}'
