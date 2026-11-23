#!/usr/bin/env bash
# pipeline_uniform.sh -- the uniform-window gate on all four columns.
#
# Detached, resumable, watched: markers in logs/chain_state_uniform/, a flock,
# a heartbeat every 5 minutes. Outputs go to a NEW timestamped folder under
# gate_check/results/uniform_<stamp>/. No record artefact is opened for writing.
#
# Phase 1 per column: the record's conductance nodes plus the true value,
# sharded across up to five concurrent processes. Phase 2 (edges, repeat
# starts, the report) is gate_check/uniform_edges.jl, run once per column.
set -u
cd "$(dirname "$0")/.." || exit 1
STATE=logs/chain_state_uniform; mkdir -p "$STATE" logs
LOG=logs/pipeline_uniform.log; HB="$STATE/heartbeat"; LOCK="$STATE/lock"
say() { echo "- $(date '+%Y-%m-%d %H:%M') $*" >> "$LOG"; }

exec 9>"$LOCK"
flock -n 9 || { say "another copy holds the lock; exiting"; exit 0; }
touch "$HB"
( while kill -0 $$ 2>/dev/null; do touch "$HB"; sleep 300; done ) 9>&- &
trap 'kill %1 2>/dev/null' EXIT

# One run folder for the whole job; reused on a resume.
if [ -f "$STATE/RUNDIR" ]; then OUT=$(cat "$STATE/RUNDIR"); else
    OUT="gate_check/results/uniform_$(date '+%Y%m%d-%H%M%S')"
    mkdir -p "$OUT"; echo "$OUT" > "$STATE/RUNDIR"
fi
say "=== uniform gate start, output $OUT ==="

NODES_3P="0.0016,0.00505964,0.016,0.0505964,0.1426,0.179523,0.2051,0.226006,0.284525,0.358195,0.450941,0.567701,0.714694"
NODES_3S="0.016,0.0284525,0.0505964,0.0899746,0.16,0.284525,0.505964,0.899746,1.1246,1.6,2.84525,5.05964,8.99746,16.0,28.4525"
NODES_2Pv2="0.016,0.0284525,0.0505964,0.0899746,0.16,0.284525,0.505964,0.7449,0.899746,1.6,2.84525,5.05964,8.99746,16.0,28.4525"
NODES_2Sv2="0.016,0.0284525,0.0505964,0.0899746,0.16,0.284525,0.505964,0.899746,1.3962,1.6,2.84525,5.05964,8.99746,16.0,28.4525"

# shard a comma list into N roughly equal pieces, round-robin so the hard
# (high-conductance) values are spread rather than piled on one process
shard() {
    python3 - "$1" "$2" "$3" <<'PY'
import sys
vals = sys.argv[1].split(","); n = int(sys.argv[2]); k = int(sys.argv[3])
print(",".join(vals[k::n]))
PY
}

NSHARD=5
for col in 3P 3S 2Pv2 2Sv2; do
    if [ -f "$STATE/done_profile_$col" ]; then say "skip profile $col"; continue; fi
    eval "NODES=\$NODES_$col"
    say "START profile $col ($NSHARD shards)"
    pids=(); ok=1
    for k in $(seq 0 $((NSHARD-1))); do
        cs=$(shard "$NODES" "$NSHARD" "$k")
        [ -z "$cs" ] && continue
        # shard-level resume: a finished shard is not redone on a relaunch
        if [ -s "$OUT/shard_${col}_${k}.jld2" ]; then
            say "  shard $col/$k already written, skipping"; continue
        fi
        julia --project=. gate_check/uniform_gate.jl --col "$col" --cs "$cs" \
              --out "$OUT" --shard "$k" > "logs/ug_${col}_$k.log" 2>&1 &
        pids+=($!)
    done
    for p in "${pids[@]}"; do wait "$p" || ok=0; done
    if [ "$ok" = 1 ]; then touch "$STATE/done_profile_$col"; say "DONE  profile $col"
    else say "FAILED profile $col -- fix and relaunch to resume"; exit 1; fi
done

for col in 3P 3S 2Pv2 2Sv2; do
    if [ -f "$STATE/done_edges_$col" ]; then say "skip edges $col"; continue; fi
    say "START edges $col"
    if julia --project=. gate_check/uniform_edges.jl --col "$col" --out "$OUT" \
         > "logs/ue_$col.log" 2>&1; then
        touch "$STATE/done_edges_$col"; say "DONE  edges $col"
    else
        say "FAILED edges $col -- fix and relaunch to resume"; exit 1
    fi
done

touch "$STATE/done_ALL"
say "=== UNIFORM GATE COMPLETE, output $OUT ==="
