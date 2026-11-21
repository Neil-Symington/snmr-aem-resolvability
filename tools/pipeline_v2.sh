#!/usr/bin/env bash
# pipeline_v2.sh -- the whole v2 sequence, detached, resumable and watched.
#
# WHY THIS EXISTS. On 21 Nov a session-bound waiter ran pair 1 of 4, then died
# without printing anything, and the machine sat idle for 22 hours before
# anyone looked. This script fixes the class:
#
#   RESUMABLE  every step writes a marker in logs/chain_state/ when it
#              succeeds, and a step whose marker exists is skipped. Relaunching
#              costs nothing and picks up where it stopped.
#   DETACHED   launched with setsid+nohup, so it outlives the session.
#   WATCHED    a heartbeat file is touched every 5 minutes while the chain
#              lives. tools/pipeline_watchdog.sh, run from cron every 15
#              minutes, relaunches the chain if that file goes stale.
#
# Run it any number of times; concurrent copies are prevented by a lock.
set -u
cd "$(dirname "$0")/.." || exit 1

STATE=logs/chain_state;  mkdir -p "$STATE" logs
LOG="logs/pipeline_v2.log"
HB="$STATE/heartbeat"
LOCK="$STATE/lock"

say() { echo "- $(date '+%Y-%m-%d %H:%M') $*" | tee -a "$LOG" >/dev/null; }

# one copy at a time; a stale lock from a killed run is taken over
exec 9>"$LOCK"
flock -n 9 || { say "another copy holds the lock; exiting"; exit 0; }

touch "$HB"
( while kill -0 $$ 2>/dev/null; do touch "$HB"; sleep 300; done ) &
HBPID=$!
trap 'kill $HBPID 2>/dev/null' EXIT

# step <name> <command...>   -- run unless already marked done
step() {
    local name="$1"; shift
    if [ -f "$STATE/done_$name" ]; then say "skip $name (already done)"; return 0; fi
    say "START $name"
    if "$@" >> "logs/step_$name.log" 2>&1; then
        touch "$STATE/done_$name"; say "DONE  $name"
    else
        say "FAILED $name (exit $?) -- chain stops; fix and relaunch to resume"
        exit 1
    fi
}

POOL5=aem2Pv2_m4c5e1,aem2Pv2_m4c5e2,aem2Pv2_m4c5e3,aem2Pv2_m4c5e4
POOL6=aem2Sv2_m4c6e1,aem2Sv2_m4c6e2,aem2Sv2_m4c6e3,aem2Sv2_m4c6e4
J="julia --project=."

say "=== pipeline start ==="

# ---- Stage 1, two ensembles at a time -------------------------------------
step pair_c5_e12  bash tools/m4_pair.sh 5 1 2
step pair_c5_e34  bash tools/m4_pair.sh 5 3 4
step pair_c6_e12  bash tools/m4_pair.sh 6 1 2
step pair_c6_e34  bash tools/m4_pair.sh 6 3 4

# ---- Stage 2, pooled and true-kernel, eight ensembles each ----------------
step stage2_2Pv2_pooled      $J stages/stage2_snmr.jl --tier typical --conds C2 \
                                --tag 2Pv2T --nens 8 --pool "$POOL5"
step stage2_2Pv2_truekernel  $J stages/stage2_snmr.jl --tier typical --conds C2 \
                                --tag 2Pv2T --nens 8 --pool "$POOL5" --truekernel
step stage2_2Sv2_pooled      $J stages/stage2_snmr.jl --tier typical --conds C2 \
                                --tag 2Sv2T --nens 8 --pool "$POOL6"
step stage2_2Sv2_truekernel  $J stages/stage2_snmr.jl --tier typical --conds C2 \
                                --tag 2Sv2T --nens 8 --pool "$POOL6" --truekernel

say "=== Stage 2 complete; Stage 3 and the gate follow ==="
say "NOTE: Stage 3 and gate_check steps are appended by the operator once the"
say "      Stage 2 run folders are known -- they are named in the invocation."
say "=== pipeline reached the end of its defined steps ==="
