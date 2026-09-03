#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Drive Overlap_CustomNulls_ByYear.R across every AUGMENT x YEAR x (LEVEL,POOL)
# pairing, concurrently, then average the two years for each pairing.
#
# Safe to parallelize: each run is an independent Rscript process that writes a
# uniquely named output (<LEVEL>_by_<POOL>_<aug|noaug>_<YEAR>_*Null.csv), so runs
# never collide -- they only share read-only inputs. The only limit is RAM/CPU,
# so concurrency is capped at $JOBS.
#
# Prereqs: run CalculateAbundance_ByYear.R once first (the per-year abundance
# files must exist), and apply the commandArgs / AUGMENT edits to both R scripts.
#
# Usage:  ./run_ByYear.sh                    # full grid, $JOBS at a time
#         JOBS=4 ./run_ByYear.sh             # cap concurrency at 4
#         AUGS="FALSE" ./run_ByYear.sh       # only the no-augmentation runs
# ---------------------------------------------------------------------------
set -euo pipefail

REPO=/home/aly/Beetles/BeetleBodySizeVariation
JOBS=${JOBS:-6}                       # max concurrent runs (override: JOBS=N ./run_ByYear.sh)

YEARS=(2018 2019)
AUGS=(${AUGS:-TRUE FALSE})            # augmentation on and/or off (override: AUGS="FALSE" ...)
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

#### STAGE 1: run the grid, up to $JOBS at a time ####
# Emit one "<LEVEL> <POOL> <YEAR> <AUGMENT>" line per job; xargs keeps $JOBS alive.
echo ">> stage 1: ${#PAIRS[@]} pairings x ${#YEARS[@]} years x ${#AUGS[@]} aug settings, $JOBS at a time"
for AUG in "${AUGS[@]}"; do
  for YEAR in "${YEARS[@]}"; do
    for PAIR in "${PAIRS[@]}"; do
      echo "$PAIR $YEAR $AUG"
    done
  done
done | xargs -P "$JOBS" -L1 bash -c '
  LEVEL=$1; POOL=$2; YEAR=$3; AUG=$4
  if [ "$AUG" = TRUE ]; then AUG_TAG=aug; else AUG_TAG=noaug; fi
  tag="${LEVEL}_by_${POOL}_${AUG_TAG}_${YEAR}"
  echo "[start] $tag"
  if Rscript Overlap_CustomNulls_ByYear.R "$LEVEL" "$POOL" "$YEAR" "$AUG" > "logs/${tag}.log" 2>&1; then
    echo "[ done] $tag"
  else
    echo "[FAIL ] $tag  (see logs/${tag}.log)"
  fi
' _

#### STAGE 2: average the two years for each (pairing x aug) (fast; serial) ####
echo ">> stage 2: averaging years per pairing"
for PAIR in "${PAIRS[@]}"; do
  set -- $PAIR; LEVEL=$1; POOL=$2
  for AUG in "${AUGS[@]}"; do
    if [ "$AUG" = TRUE ]; then AUG_TAG=aug; else AUG_TAG=noaug; fi
    tag="${LEVEL}_by_${POOL}_${AUG_TAG}_avg"
    echo "[avg  ] $LEVEL $POOL $AUG_TAG"
    Rscript Overlap_CustomNulls_ByYearAverage.R "$LEVEL" "$POOL" "$AUG" > "logs/${tag}.log" 2>&1 \
      || echo "[FAIL ] avg $LEVEL $POOL $AUG_TAG  (see logs/${tag}.log)"
  done
done

echo ">> all done. Outputs in ./Outputs, logs in ./logs"
