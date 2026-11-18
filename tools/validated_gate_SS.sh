#!/usr/bin/env bash
# validated_gate_SS.sh -- the VALIDATED gate (gate_check, tag
# gate-of-record-20261119) on 3S and 2S, with the settings that produced the
# 3P and 2P results: lambda2 {0.01, 0.1}, STATED limits 0.001-2 S/m
# (mbounds -3.0,0.30103), warm starts both sweep directions, windows
# 7.5-23.5 m (3S) and 17.5-33.5 m (2S).
#
# Two steps, as for 3P and 2P: a coarse forced-conductance profile bound to
# bound, then the warm-started bisection on its stored nodes.
set -u -o pipefail
cd "$(dirname "$0")/.." || exit 1
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

win_3S="--top 7.5  --base 23.5"
win_2S="--top 17.5 --base 33.5"
MB="-3.0,0.30103"
LAM="0.01,0.1"

echo "=== $(date +%H:%M) waiting for Stage 0 on 3S and 2S"
while :; do
    n=0
    for t in 3S 2S; do [ -s "outputs/synth_${t}/stage0_synth_${t}.jld2" ] && n=$((n+1)); done
    [ "$n" = 2 ] && break
    sleep 20
done
echo "=== $(date +%H:%M) Stage 0 present"

one() {
    local t="$1" win
    eval "win=\$win_${t}"
    local s0="outputs/synth_${t}/stage0_synth_${t}.jld2"
    echo "=== $(date +%H:%M) $t coarse profile, bound to bound"
    julia --project=. gate_check/profile_minchi2.jl \
        --stage0 "$s0" $win --lo -3.0 --hi 0.30103 --step 0.25 \
        --mbounds "$MB" --lam "$LAM" --starts warm \
        --truth "data/truth_murray${t}.jld2" --tag "${t}_sed" \
        > "logs/vgate_${t}_profile.log" 2>&1
    echo "=== $(date +%H:%M) $t profile rc=$?"
    echo "=== $(date +%H:%M) $t edge bisection"
    julia --project=. gate_check/edge_bisect.jl \
        --profile "gate_check/results/profile_${t}_sed.jld2" \
        --stage0 "$s0" $win --mbounds "$MB" --lam "$LAM" \
        --maxrefine 6 --seed 20261116 --tag "${t}_sed" \
        > "logs/vgate_${t}_bisect.log" 2>&1
    echo "=== $(date +%H:%M) $t bisect rc=$?"
}

one 3S & one 2S & wait
echo "=== $(date +%H:%M) VALIDATED_GATE_SS_DONE"
for t in 3S 2S; do
    echo "--- $t profile ---"; tail -14 "logs/vgate_${t}_profile.log"
    echo "--- $t bisect ---";  tail -12 "logs/vgate_${t}_bisect.log"
done
