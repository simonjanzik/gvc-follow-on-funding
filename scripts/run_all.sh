#!/usr/bin/env bash
# run_all.sh: the whole pipeline in order, from any working directory. Usage: ./run_all.sh [--from NN]
# Stages 40 and 60 run for hours (bootstrap, simulation); on macOS caffeinate keeps the machine awake, both resume
# from their draw files after an interruption. Python is pinned to the interpreter the results were made with.
set -euo pipefail
cd "$(dirname "$0")"
PY=/usr/bin/python3
OUT=../outputs            # same folder as OUTPUT_DIR in config.py
KEEP=""; command -v caffeinate >/dev/null && KEEP="caffeinate -i"
FROM=10; [[ "${1:-}" == "--from" ]] && FROM="${2:?usage: ./run_all.sh [--from NN]}"
run() { local nn=$1; shift; (( nn >= FROM )) || return 0; echo "== $* =="; "$@"; }
$PY config.py > /dev/null      # config.R follows config.py before any R stage reads it
run 10 $PY 10_build_data.py
run 20 $PY 20_descriptives.py
run 30 $KEEP Rscript 30_estimators.R
run 30 $PY _values.py merge "$OUT/data/30_values.csv" 30_estimators.R
run 40 $KEEP Rscript 40_inference.R
run 40 $PY _values.py merge "$OUT/data/40_values.csv" 40_inference.R
run 50 $PY 50_sensitivity.py
run 60 $KEEP Rscript 60_simulation.R
run 60 $PY _values.py merge "$OUT/data/60_values.csv" 60_simulation.R
$PY _env.py > /dev/null
echo "done"
