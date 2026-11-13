#!/usr/bin/env bash
# Run Stage 1 as four INDEPENDENT ensembles for one column under one setting.
#
#   ./run_stage1_setting.sh <setting> <tag> [nsamples] [nchains]
#
# The chains are unseeded (tools/aem_HYP.jl documents why), so four repeats of
# one command are four independent ensembles. That is what makes a between-
# ensemble R-hat possible at all.
#
# Sampler knobs come from the environment and are read by tools/aem_HYP.jl:
#   SCT_NMAX           max nodes k      (production 40)
#   SCT_TMAX           max temperature  (production 2.5)
#   SCT_NCHAINSATONE   chains at T=1    (production 1)
#
# Output: data/aem<tag>_<setting>e<N>/ -- never data/aem<tag>/, so no
# production posterior can be overwritten.
set -u -o pipefail

SETTING="${1:?setting name, e.g. b1}"
TAG="${2:?column tag, e.g. 2P}"
NSAMPLES="${3:-200001}"
NCHAINS="${4:-4}"

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
LOGDIR="$HERE/logs"
mkdir -p "$LOGDIR"

# Two ensembles at a time: each takes NCHAINS workers and this machine has 8
# threads. Five concurrent trans-D runs exhaust RAM on this box.
CONC=2
if [ "$NCHAINS" -gt 4 ]; then CONC=1; fi

echo "setting=$SETTING tag=$TAG nsamples=$NSAMPLES nchains=$NCHAINS conc=$CONC"
echo "  SCT_NMAX=${SCT_NMAX:-40} SCT_TMAX=${SCT_TMAX:-2.5} SCT_NCHAINSATONE=${SCT_NCHAINSATONE:-1}"

cd "$REPO"
t_all=$(date +%s)
pids=""
for e in 1 2 3 4; do
    LABEL="_${SETTING}e${e}"
    LOG="$LOGDIR/stage1_${TAG}_${SETTING}e${e}.log"
    echo "  -> ensemble $e  label=$LABEL  log=$LOG"
    ( t0=$(date +%s)
      julia --project=. tools/aem_HYP.jl "$NSAMPLES" "$TAG" "$NCHAINS" "$LABEL" > "$LOG" 2>&1
      rc=$?
      t1=$(date +%s)
      echo "ENSEMBLE_DONE e=$e rc=$rc minutes=$(( (t1-t0)/60 ))" >> "$LOG"
    ) &
    pids="$pids $!"
    # throttle
    while [ "$(jobs -r | wc -l)" -ge "$CONC" ]; do wait -n 2>/dev/null || break; done
done
wait
t_end=$(date +%s)
echo "ALL_DONE setting=$SETTING tag=$TAG total_minutes=$(( (t_end-t_all)/60 ))"
for e in 1 2 3 4; do
    grep -h "ENSEMBLE_DONE" "$LOGDIR/stage1_${TAG}_${SETTING}e${e}.log" 2>/dev/null || echo "e=$e NO RESULT"
done
