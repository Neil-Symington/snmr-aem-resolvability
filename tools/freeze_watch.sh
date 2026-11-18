#!/usr/bin/env bash
# freeze_watch.sh -- wait for 3P and 2P Stage 2, repoint their configs at the
# run folders, then run the gate freeze through the launcher.
#
# Stage 2 output now lives in data/stage2/<tag>_<kernel>_<timestamp>/, so each
# column's [data].snmr_ensemble has to name the run that was actually used. The
# config validator checks that file exists, which is what blocked the freeze
# before Stage 2 on these columns had been run -- Stage 4 itself never reads it.
set -u
cd "$(dirname "$0")/.." || exit 1
LOG() { echo "$(date '+%H:%M') $*"; }

LOG "waiting for Stage 2 on 3P and 2P"
while :; do
    n=0
    for t in 3P 2P; do
        d=$(ls -d data/stage2/${t}T_pooled_* 2>/dev/null | sort | tail -1)
        [ -n "$d" ] && [ -s "$d/cond_C2_typical${t}T.jld2" ] && n=$((n+1))
    done
    [ "$n" = 2 ] && break
    sleep 30
done
LOG "both Stage 2 runs present"

for t in 3P 2P; do
    d=$(ls -d data/stage2/${t}T_pooled_* | sort | tail -1)
    f="$d/cond_C2_typical${t}T.jld2"
    sed -i "s|^snmr_ensemble = .*|snmr_ensemble = \"${f}\"   # eight-ensemble pooled run, 18 Nov|" \
        "configs/synth_${t}.toml"
    LOG "configs/synth_${t}.toml -> $f"
done

LOG "launching the gate freeze"
exec bash tools/safe_launch.sh \
    logs/gate_freeze_smoke.log logs/gate_freeze.log \
    'bash tools/gate_freeze_run.sh --smoke' \
    'bash tools/gate_freeze_run.sh'
