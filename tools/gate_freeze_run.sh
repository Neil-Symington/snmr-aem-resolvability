#!/usr/bin/env bash
# gate_freeze_run.sh -- the release Stage 4 on 3P and 2P, then the acceptance test.
#
#   bash tools/gate_freeze_run.sh [--smoke]
#
# Runs the RELEASE entry point, not a gate_check driver: this is what the freeze
# is a freeze of. --smoke uses a tiny node count to exercise every save path.
set -u -o pipefail
cd "$(dirname "$0")/.." || exit 1
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

SMOKE=0; [ "${1:-}" = "--smoke" ] && SMOKE=1
EXTRA=""
if [ "$SMOKE" = 1 ]; then EXTRA="--grid single_pass --npts 3"; fi

# THE WINDOW THE TARGETS WERE ESTABLISHED ON, forced explicitly. Stage 3's
# window moved after the 18 Nov SNMR fix (saturated water criterion -> the
# window comes from a tie-break), so an acceptance run must force the window its
# target numbers came from or it compares two different problems. Production
# takes Stage 3's window; this path says so in the output and the test asserts it.
win_3P="--top 7.5 --base 23.5"
win_2P="--top 17.5 --base 33.5"

# The two columns run CONCURRENTLY. Each is a single-threaded chain of occam
# solves, so two processes sit comfortably on six cores and the wall clock
# halves. Logs are per column; a failure in either fails the script.
t0=$(date +%s)
pids=""
for col in 3P 2P; do
    echo "=== $(date '+%H:%M') Stage 4 release run, $col ${EXTRA:+(smoke)}"
    eval "win=\$win_${col}"
    julia --project=. stages/stage4_gate_v2.jl --config "configs/synth_${col}.toml" $win $EXTRA \
        > "logs/stage4_${col}${EXTRA:+_smoke}.log" 2>&1 &
    pids="$pids $!"
done
rc=0
for p in $pids; do wait "$p" || rc=1; done
echo "=== $(date '+%H:%M') both columns rc=$rc elapsed_min=$(( ($(date +%s)-t0)/60 ))"
for col in 3P 2P; do
    echo "--- $col, last lines ---"
    tail -12 "logs/stage4_${col}${EXTRA:+_smoke}.log"
done
[ "$rc" -eq 0 ] || exit $rc

if [ "$SMOKE" = 0 ]; then
    echo "=== $(date '+%H:%M') acceptance test"
    julia --project=. test/test_gate_acceptance.jl
    echo "=== acceptance rc=$?"
fi
