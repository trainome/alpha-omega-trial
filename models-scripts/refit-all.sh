#!/usr/bin/env bash
# Full refit of all observed-data models by running their model-*.R scripts,
# two at a time (4 chains each):
#   * s7 to s22, including the D = 6 variants of s14 to s17 (metabolic syndrome)
#   * m1 (drop-out/missingness), m2 to m6 (fatty acids)
#   * mo1 to mo4 (OGTT)
# The m and mo scripts also refit their prior predictive fits. Models fitted to
# synthetic data (04-metsyn-models.qmd, 06-metsyn-simulations.qmd) are not
# included. Output from each script goes to
# models-scripts/<model>-run.log, progress is printed to the console. Run from
# the project root:
#
#   bash models-scripts/refit-all.sh 2>&1 | tee models-scripts/refit-all.log
#
# The models are split into two fixed lanes that run sequentially and in
# parallel to each other. The D = 6 scripts compile to the same executable as
# their D = 4 counterparts (e.g. model-s14.R and model-s14-D6.R both write
# models/s14.exe), each pair is therefore kept in the same lane so they never
# run at the same time. Lanes are balanced on approximate sampling times
# (minutes, three lanes, K = 12): s14 83, s15 81, s12 83, s17 71, s20 54,
# s10 54, s19 45, s16 42, s13 32, s8 22, s11 12, s18 10, s7 9, s9 8. D = 6
# variants and s21/s22 are longer. The m and mo models are short (sampling,
# posterior only: mo3 10, mo4 10, mo2 5, mo1 2, m1-m6 < 1) and run last.
# Expect roughly 9-10 hours in total.
#
# The m and mo scripts can be smoke tested (20 iterations, fits saved as
# fits/paper1/smoke.*.RDS) with: REFIT_SMOKE=1 Rscript models-scripts/model-m1.R

set -u
cd "$(dirname "$0")/.."

# Rscript may not be on PATH in a plain Git Bash session (R is not added to the
# Windows PATH by default), fall back to the newest per-user R installation.
if ! command -v Rscript > /dev/null 2>&1; then
  r_bin=$(ls -d "$LOCALAPPDATA"/Programs/R/R-*/bin/x64 2>/dev/null | sort -V | tail -n 1)
  if [ -z "$r_bin" ]; then
    echo "Rscript not found, add R's bin/x64 directory to PATH" >&2
    exit 1
  fi
  export PATH="$(cygpath -u "$r_bin"):$PATH"
fi
echo "Using $(command -v Rscript)"

lane_a=(s14 s14-D6 s17 s17-D6 s12 s19 s13 s11 s7 s21 mo3 mo1 m1 m3 m5)
lane_b=(s15 s15-D6 s16 s16-D6 s20 s10 s22 s8 s18 s9 mo4 mo2 m2 m4 m6)

failed_file=$(mktemp)

run_one() {
  m=$1
  start=$(date +%s)
  echo "$(date '+%F %T') start $m"
  if Rscript "models-scripts/model-$m.R" > "models-scripts/$m-run.log" 2>&1; then
    status=ok
  else
    status=FAILED
    echo "$m" >> "$failed_file"
  fi
  echo "$(date '+%F %T') $status $m ($(( ($(date +%s) - start) / 60 )) min)"
}

run_lane() {
  for m in "$@"; do
    run_one "$m"
  done
}

# Record the data version used for this refit (participants in the Stan data).
Rscript -e 'suppressMessages({source("./R/packages.R"); source("./R/data-prep.R")});
  cat("AOData", as.character(packageVersion("AOData")),
      "| N_P =", data_load_s7s10()$d$N_P, "\n")' 2>/dev/null

total_start=$(date +%s)

run_lane "${lane_a[@]}" &
run_lane "${lane_b[@]}" &
wait

echo "$(date '+%F %T') done ($(( ($(date +%s) - total_start) / 60 )) min)"
if [ -s "$failed_file" ]; then
  echo "Failed: $(tr '\n' ' ' < "$failed_file")"
  rm -f "$failed_file"
  exit 1
fi
rm -f "$failed_file"
