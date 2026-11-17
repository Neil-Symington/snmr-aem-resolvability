#!/usr/bin/env bash
# m4_pair_watch.sh -- hand the M4 inversions over from serial to two-at-a-time,
# decide whether that was worth it, and carry the decision into column 4.
#
# Neil, 17 Nov: "after tonight's 2P edge job finishes, run two ensembles at a
# time. Time the first pair. Keep it only if total throughput improves, i.e.
# each ensemble takes less than about 1.7x the single-run time. Otherwise go
# back to one. Report which."
#
# The switch happens at an ENSEMBLE BOUNDARY, not the moment the edge job ends.
# The edge job finishes around 19:30 and column 3's e2 around 22:00, so killing
# the serial runner when the edge job ends would throw away about an hour of e2.
# Waiting for e2 to land costs nothing and gives two clean single-run baselines
# instead of one.
#
# A serial runner killed mid-ensemble leaves a partial Stage 1 directory. HiQGA
# appends to its history files, so re-running into one would corrupt it. Partial
# directories are MOVED ASIDE, never deleted.
set -u
cd "$(dirname "$0")/.." || exit 1
REPO="$(pwd)"
LOG() { echo "$(date '+%H:%M') $*"; }
NOTE() { echo "- $(date '+%Y-%m-%d %H:%M') $*" >> "$REPO/M4_STATUS.md"; }

EDGE_LOG=logs/gatework_20261117b.log
MANIFEST=data/m4/manifest.toml
DECISION=data/m4/CONCURRENCY_DECISION.txt

mins_of() {   # mins_of <column> <ensemble>
    awk -v s="[column_$1.stage1_e$2]" '$0==s{f=1} f&&/^minutes/{print $3; exit}' "$MANIFEST"
}
done_ens() {  # done_ens <column> <ensemble>
    grep -q "^\[column_$1.stage1_e$2\]" "$MANIFEST" 2>/dev/null
}
stop_serial() {
    local p c
    p=$(ps -eo pid,args | grep -E "bash tools/m4_run\.sh" | grep -v grep | awk '{print $1}')
    if [ -n "$p" ]; then LOG "stopping the serial runner (pid $p)"; kill $p 2>/dev/null; sleep 2; fi
    c=$(ps -eo pid,args | grep "aem_HYP.jl" | grep -v grep | awk '{print $1}')
    if [ -n "$c" ]; then LOG "stopping its ensemble (pid $c)"; kill $c 2>/dev/null; sleep 5; fi
}
park_partial() {  # park_partial <tag> <column> <ensemble...>
    local tag="$1" col="$2"; shift 2
    local E D f ok NEW
    for E in "$@"; do
        D="data/aem${tag}_m4c${col}e${E}"
        [ -d "$D" ] || continue
        ok=1
        for f in aem_setup.jld2 aem_posterior.jld2 misfits_hypsynth_aem_s.bin \
                 models_hypsynth_aem_s.bin points_hypsynth_aem_s.bin; do
            [ -s "$D/$f" ] || ok=0
        done
        if [ "$ok" = 1 ]; then LOG "$D is complete; it will be skipped"
        else NEW="${D}_partial_killed_$(date +%H%M)"; mv "$D" "$NEW"
             LOG "moved partial $D aside to $NEW (not deleted)"; fi
    done
}

# ── 1. wait for the 2P edge job ──────────────────────────────────────────
LOG "waiting for the 2P edge job to finish ($EDGE_LOG)"
while ! grep -q "gate work finished" "$EDGE_LOG" 2>/dev/null; do sleep 20; done
LOG "2P edge job finished"

# ── 2. wait for column 3 e2, the second single-run baseline ──────────────
LOG "waiting for column 3 e2 to land in the manifest (switching at an ensemble boundary)"
while ! done_ens 3 2; do sleep 20; done
B1=$(mins_of 3 1); B2=$(mins_of 3 2)
LOG "single-run baselines: e1 = ${B1} min, e2 = ${B2} min"

stop_serial
park_partial 3S 3 3 4

# ── 3. the pair: column 3, e3 + e4 ───────────────────────────────────────
LOG "launching column 3 e3 + e4 concurrently through safe_launch.sh"
bash tools/safe_launch.sh \
    logs/m4_c3_pair_smoke.log logs/m4_c3_pair_real.log \
    'bash tools/m4_pair.sh --smoke 3 3 4' \
    'bash tools/m4_pair.sh 3 3 4'
# safe_launch nohups the real run; wait for it to land in the manifest
LOG "waiting for the pair to finish"
while ! done_ens 3 4; do sleep 20; done
PAIRW=$(mins_of 3 4)       # wall clock of the whole group

# ── 4. the decision ──────────────────────────────────────────────────────
BASE=$(python3 -c "print(($B1+$B2)/2.0)")
RATIO=$(python3 -c "print(f'{$PAIRW/$BASE:.2f}')")
KEEP=$(python3 -c "print('pair' if $PAIRW/$BASE < 1.7 else 'serial')")
{
    echo "Concurrency experiment, column 3, $(date '+%Y-%m-%d %H:%M')"
    echo "  single-run baselines: e1 = $B1 min, e2 = $B2 min, mean = $BASE min"
    echo "  two at a time:        e3 + e4 in $PAIRW min wall clock"
    echo "  ratio per ensemble:   ${RATIO}x the single-run time (keep if < 1.7)"
    echo "  DECISION:             $KEEP"
    echo
    echo "Each ensemble in the pair occupies the wall clock of the whole group,"
    echo "so a ratio below 2.0 already means more ensembles per hour, and 1.7 is"
    echo "the margin Neil set. Caveat: e1 ran alongside the 2P edge job, so the"
    echo "baselines are slightly pessimistic, which biases this toward 'pair'."
} > "$DECISION"
LOG "decision: $KEEP (ratio ${RATIO}x)"
NOTE "concurrency decision: $KEEP, ${RATIO}x the single-run time (baselines $B1/$B2 min, pair $PAIRW min)"

# ── 5. column 4, at the winning concurrency ──────────────────────────────
LOG "column 4: Stage 0"
M4_STAGE0_ONLY=1 bash tools/m4_run.sh 4 >> logs/m4_c4_stage0.log 2>&1 || {
    LOG "column 4 Stage 0 FAILED, see logs/m4_c4_stage0.log"; exit 1; }
if [ "$KEEP" = "pair" ]; then
    LOG "column 4: two pairs, e1+e2 then e3+e4"
    bash tools/m4_pair.sh 4 1 2 >> logs/m4_c4_pair.log 2>&1
    bash tools/m4_pair.sh 4 3 4 >> logs/m4_c4_pair.log 2>&1
else
    # One at a time, but still through m4_pair.sh with a group of one, so that
    # Stage 0 is not run a second time. Calling m4_run.sh here would re-run it
    # and append a SECOND [column_4] table to the manifest, which would make the
    # file invalid TOML.
    LOG "column 4: one at a time"
    for E in 1 2 3 4; do
        bash tools/m4_pair.sh 4 "$E" >> logs/m4_c4_serial.log 2>&1
    done
fi
LOG "column 4 finished"
NOTE "column 4 finished at concurrency: $KEEP"
