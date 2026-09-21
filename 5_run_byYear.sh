#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Drive CommunityShape_CustomNulls_ByYear.R across every YEAR x (LEVEL, POOL)
# pairing, concurrently. No augmentation dimension and no year-averaging stage:
# each run is one Rscript process writing a uniquely named set of outputs
# (<LEVEL>_by_<POOL>_<YEAR>_{Pool,Individual,SwapMeans}Null.csv), so runs never
# collide -- they only share read-only inputs. Concurrency is capped at $JOBS.
#
# Prereqs: run CalculateAbundance_ByYear.R once first (the per-year abundance
# files must exist), and keep community_depth.R beside the R script (it is
# sourced).
#
# Usage:  ./run_byYear.sh                 # full grid, $JOBS at a time
#         JOBS=4 ./run_byYear.sh          # cap concurrency at 4
#         YEARS="2018" ./run_byYear.sh    # only one year
# ---------------------------------------------------------------------------
set -euo pipefail

REPO=/home/aly/Beetles/BeetleBodySizeVariation
JOBS=${JOBS:-16}                       # max concurrent runs (override: JOBS=N ./run_byYear.sh)

YEARS=(${YEARS:-2018 2019})            # override: YEARS="2018" ./run_byYear.sh
# Valid focal->pool pairings only: a pool must sit above the focal level, so
# site->site is degenerate and omitted. Edit this list to taste.
PAIRS=(
  "plot site"
  "plot domain"
  "plot all"
  "site domain"
  "site all"
)

cd "$REPO"
mkdir -p logs

#### run the grid, up to $JOBS at a time ####
# Emit one "<LEVEL> <POOL> <YEAR>" line per job; xargs keeps $JOBS alive.
echo ">> ${#PAIRS[@]} pairings x ${#YEARS[@]} years, $JOBS at a time"
for YEAR in "${YEARS[@]}"; do
  for PAIR in "${PAIRS[@]}"; do
    echo "$PAIR $YEAR"
  done
done | xargs -P "$JOBS" -L1 bash -c '
  LEVEL=$1; POOL=$2; YEAR=$3
  tag="${LEVEL}_by_${POOL}_${YEAR}"
  echo "[start] $tag"
  if Rscript CommunityShape_CustomNulls_ByYear.R "$LEVEL" "$POOL" "$YEAR" > "logs/${tag}.log" 2>&1; then
    echo "[ done] $tag"
  else
    echo "[FAIL ] $tag  (see logs/${tag}.log)"
  fi
' _

echo ">> all done. Outputs in ./Outputs, logs in ./logs"
