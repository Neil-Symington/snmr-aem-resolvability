#!/usr/bin/env bash
# stage2_run.sh -- Stage 2 for one column, full length, kernel from --pool.
#
#   bash tools/stage2_run.sh <TAG> [--control] [--smoke]
#
# TAG is the plain column, e.g. 3S. Stage 2 is invoked with <TAG>T, which is the
# data file tools/smr_data.jl writes.
#
# --pool gives the inversion kernel the median of the POOLED post-burn-in T=1
# draws across that column's four Stage 1 ensembles, instead of one folder's
# aem_posterior.jld2. That is the inversion kernel only; the DATA are generated
# on the true conductivity (fixed 18 Nov, tag snmr-data-truth-20261118).
#
# --control inverts the SAME data with the kernel on the TRUE conductivity, to
# separate what the SNMR inversion loses to an imperfect kernel from what it
# loses to noise and non-uniqueness. Nothing downstream uses the control.
set -u -o pipefail
cd "$(dirname "$0")/.." || exit 1
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TAG="${1:?usage: $0 <TAG> [--control] [--smoke]}"; shift || true
CONTROL=0; SMOKE=0
for a in "$@"; do
    case "$a" in
        --control) CONTROL=1 ;;
        --smoke)   SMOKE=1 ;;
        *) echo "unknown argument $a" >&2; exit 2 ;;
    esac
done
NS=100001; [ "$SMOKE" = 1 ] && NS=2001

case "$TAG" in 3P) COL=1 ;; 2P) COL=2 ;; 3S) COL=3 ;; 2S) COL=4 ;;
  *) echo "unknown tag $TAG" >&2; exit 2 ;; esac

# Stage 2 resolves a relative --pool entry against data/, so the entries are
# BARE FOLDER NAMES. Passing "data/aem..." gives data/data/aem... and fails.
POOL=""
for E in 1 2 3 4; do
    d="aem${TAG}_m4c${COL}e${E}"
    [ -s "data/$d/aem_posterior.jld2" ] || { echo "missing data/$d/aem_posterior.jld2" >&2; exit 1; }
    POOL="${POOL:+$POOL,}$d"
done

[ -s "data/smr_data_murray${TAG}T.jld2" ] || {
    echo "missing data/smr_data_murray${TAG}T.jld2 -- run tools/smr_data.jl --tag $TAG" >&2; exit 1; }

echo "=== $(date '+%H:%M') Stage 2, $TAG, control=$CONTROL, nsamples=$NS"
echo "    pool: $POOL"
t0=$(date +%s)
if [ "$CONTROL" = 1 ]; then
    julia --project=. stages/stage2_snmr.jl --tier typical --conds C2 \
          --tag "${TAG}T" --truekernel --nsamples "$NS"
else
    julia --project=. stages/stage2_snmr.jl --tier typical --conds C2 \
          --tag "${TAG}T" --pool "$POOL" --nsamples "$NS"
fi
rc=$?
echo "=== $(date '+%H:%M') Stage 2 $TAG control=$CONTROL rc=$rc elapsed_min=$(( ($(date +%s)-t0)/60 ))"
exit $rc
