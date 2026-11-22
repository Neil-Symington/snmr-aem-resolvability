#!/usr/bin/env bash
# pipeline_v2b.sh -- the whole v2 sequence, detached, resumable and watched.
#
# Supersedes pipeline_v2.sh, which stopped after Stage 2. Same markers and
# same lock, so it resumes whatever that one completed.
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
# 9>&- so the heartbeat does NOT inherit the lock file descriptor. Without it
# the subshell keeps the lock after the pipeline exits -- its `sleep` holds it
# for up to five minutes -- and the next launch is refused with "another copy
# holds the lock". Seen 22 Nov 2026.
( while kill -0 $$ 2>/dev/null; do touch "$HB"; sleep 300; done ) 9>&- &
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
# --truekernel builds the kernel from the TRUE conductivity, so --pool is
# not just unnecessary, it is refused: the two are alternatives
# (stage2_snmr.jl:318). The ensemble count is set by --nens either way.
step stage2_2Pv2_truekernel  $J stages/stage2_snmr.jl --tier typical --conds C2 \
                                --truekernel --tag 2Pv2T --nens 8
step stage2_2Sv2_pooled      $J stages/stage2_snmr.jl --tier typical --conds C2 \
                                --tag 2Sv2T --nens 8 --pool "$POOL6"
# --truekernel builds the kernel from the TRUE conductivity, so --pool is
# not just unnecessary, it is refused: the two are alternatives
# (stage2_snmr.jl:318). The ensemble count is set by --nens either way.
step stage2_2Sv2_truekernel  $J stages/stage2_snmr.jl --tier typical --conds C2 \
                                --truekernel --tag 2Sv2T --nens 8

# ---- point each config at the Stage 2 run that just finished ------------
# The run folder carries a timestamp, so it cannot be written into the config
# in advance. newest_pooled picks the most recent pooled folder for the tag;
# the marker makes the edit happen once.
newest_pooled() { ls -d data/stage2/$1_pooled_* 2>/dev/null | sort | tail -1; }

point_config() {
    local tag=$1 cfg=$2 run
    run=$(newest_pooled "$tag")
    [ -n "$run" ] || { echo "no pooled run for $tag"; return 1; }
    local ens="$run/cond_C2_typical$tag.jld2"
    [ -s "$ens" ] || { echo "missing ensemble file: $ens"; return 1; }
    python3 - "$cfg" "$ens" <<'PY'
import sys, re
cfg, ens = sys.argv[1], sys.argv[2]
s = open(cfg).read()
s2 = re.sub(r'^snmr_ensemble = .*$', f'snmr_ensemble = "{ens}"', s, count=1, flags=re.M)
assert s2 != s, "snmr_ensemble line not found in " + cfg
open(cfg, "w").write(s2)
print("pointed", cfg, "->", ens)
PY
}

step point_2Pv2  point_config 2Pv2T configs/synth_2Pv2.toml
step point_2Sv2  point_config 2Sv2T configs/synth_2Sv2.toml

# ---- Stage 3 ------------------------------------------------------------
step stage3_2Pv2  $J stages/stage3_hypothesis.jl --config configs/synth_2Pv2.toml
step stage3_2Sv2  $J stages/stage3_hypothesis.jl --config configs/synth_2Sv2.toml

# ---- the validated gate_check drivers -----------------------------------
# Settings copied from the 2S run of record (profile_2S_sed.jld2): a uniform
# grid -3.0 to 0.25 at 0.25 dex, stated limits, the fixed lambda2 set. The
# window is test two's, 17.5-33.5 m, unchanged in v2.
GC_COMMON="--top 17.5 --base 33.5 --mbounds -3.0,0.30103 --lam 0.01,0.1"

step profile_2Pv2  $J gate_check/profile_minchi2.jl --stage0 outputs/synth_2Pv2/stage0_synth_2Pv2.jld2 \
                      --lo -3.0 --hi 0.25 --step 0.25 $GC_COMMON --tag 2Pv2_sed
step profile_2Sv2  $J gate_check/profile_minchi2.jl --stage0 outputs/synth_2Sv2/stage0_synth_2Sv2.jld2 \
                      --lo -3.0 --hi 0.25 --step 0.25 $GC_COMMON --tag 2Sv2_sed

step edges_2Pv2    $J gate_check/edge_bisect.jl --profile gate_check/results/profile_2Pv2_sed.jld2 \
                      --stage0 outputs/synth_2Pv2/stage0_synth_2Pv2.jld2 $GC_COMMON \
                      --maxrefine 6 --seed 20261116 --tag 2Pv2_sed
step edges_2Sv2    $J gate_check/edge_bisect.jl --profile gate_check/results/profile_2Sv2_sed.jld2 \
                      --stage0 outputs/synth_2Sv2/stage0_synth_2Sv2.jld2 $GC_COMMON \
                      --maxrefine 6 --seed 20261116 --tag 2Sv2_sed

# ---- the release gate, refined near each edge ---------------------------
# 15-node two_pass grid, refined to 0.05 dex within 0.3 dex of each edge
# (ruling of 21 Nov). 3P runs beside 2Pv2 so the two are comparable.
step relgate_2Pv2  $J stages/stage4_gate_v2.jl --config configs/synth_2Pv2.toml \
                      --top 17.5 --base 33.5 --grid two_pass --refinedex 0.05 --refinehalf 0.3
step relgate_3P    $J stages/stage4_gate_v2.jl --config configs/synth_3P.toml \
                      --top 7.5 --base 23.5 --grid two_pass --refinedex 0.05 --refinehalf 0.3

touch "$STATE/done_ALL"
say "=== PIPELINE COMPLETE: Stage 1-3, gate_check and the release gate ==="
