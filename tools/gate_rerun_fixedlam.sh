#!/usr/bin/env bash
# gate_rerun_fixedlam.sh -- re-run the release gate under the FIXED lambda2 set
# {0.01, 0.1}, both sweep directions, then the acceptance test.
#
#   bash tools/gate_rerun_fixedlam.sh [--smoke]
#
# 2P goes first, because its lambda-hat run has already finished and is archived
# for the record. 3P follows once its lambda-hat run exits -- that run is kept
# too, reported, and then superseded.
set -u -o pipefail
cd "$(dirname "$0")/.." || exit 1
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

SMOKE=0; [ "${1:-}" = "--smoke" ] && SMOKE=1
EXTRA=""; SUF=""
if [ "$SMOKE" = 1 ]; then EXTRA="--grid single_pass --npts 3"; SUF="_smoke"; fi

win_3P="--top 7.5 --base 23.5"
win_2P="--top 17.5 --base 33.5"

run_col() {
    local col="$1" win
    eval "win=\$win_${col}"
    echo "=== $(date '+%H:%M') Stage 4 fixed-lambda run, $col ${EXTRA:+(smoke)}"
    local t0; t0=$(date +%s)
    julia --project=. stages/stage4_gate_v2.jl --config "configs/synth_${col}.toml" \
        $win $EXTRA > "logs/stage4_${col}_fixedlam${SUF}.log" 2>&1
    local rc=$?
    echo "=== $(date '+%H:%M') $col rc=$rc elapsed_min=$(( ($(date +%s)-t0)/60 ))"
    grep -E "WINDOW OVERRIDDEN|absolute rule|lower edge|upper edge" \
        "logs/stage4_${col}_fixedlam${SUF}.log" || true
    return $rc
}

run_col 2P || exit $?

if [ "$SMOKE" = 0 ]; then
    echo "=== $(date '+%H:%M') waiting for the lambda-hat 3P run to exit"
    while pgrep -f "stage4_gate_v2.jl --config configs/synth_3P.toml" > /dev/null; do sleep 60; done
    echo "=== $(date '+%H:%M') it has exited; its edges, for the record:"
    grep -E "lower edge|upper edge" logs/stage4_3P.log || echo "  (none recorded)"
    cp -f outputs/synth_3P/stage4_synth_3P.jld2 \
          outputs/record_lamhat/stage4_synth_3P_lamhat.jld2 2>/dev/null || true
fi

run_col 3P || exit $?

if [ "$SMOKE" = 0 ]; then
    echo "=== $(date '+%H:%M') acceptance test, targets unchanged"
    julia --project=. test/test_gate_acceptance.jl
    echo "=== acceptance rc=$?"
fi
