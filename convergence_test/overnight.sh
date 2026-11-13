#!/usr/bin/env bash
# Overnight queue, serialised behind B1.
#
# B1 is the control and its ensembles are the one thing here that cannot be
# cheaply redone, so nothing new starts until it is finished. Then the B1
# diagnostics run (that is the actual Part A deliverable), then the remaining
# Occam site runs.
#
# The Occam 2P site run already in flight is left alone: it is deep into the
# gate grid and is one process against B1's eight workers.
set -u
REPO=/home/symo/Documents/GA/AEM_SMR_Resolvability_Quickrun/release
cd "$REPO" || exit 1
B1LOG="$REPO/convergence_test/logs/driver_b1_all.log"
STATUS="$REPO/MORNING_STATUS.md"

say() { echo "[$(date -Is)] $*"; }

{
  echo "# Overnight status"
  echo
  echo "Started $(date -Is). This file is rewritten as each step finishes."
  echo
} > "$STATUS"

note() { echo "$*" >> "$STATUS"; }

say "waiting for B1 to finish (marker B1_ALL_DONE in driver_b1_all.log)"
while ! grep -q "B1_ALL_DONE" "$B1LOG" 2>/dev/null; do sleep 60; done
say "B1 finished"
note "## B1 (the control) — FINISHED"
note ""
note '```'
grep -h "ENSEMBLE_DONE" "$REPO"/convergence_test/logs/stage1_2P_b1e*.log 2>/dev/null | sed 's/^/2P /' >> "$STATUS"
grep -h "ENSEMBLE_DONE" "$REPO"/convergence_test/logs/stage1_3P_b1e*.log 2>/dev/null | sed 's/^/3P /' >> "$STATUS"
note '```'
note ""

# ── the Part A / B1 deliverable ────────────────────────────────────────────
for TAG in 2P 3P; do
    say "diagnosing B1 $TAG"
    L="$REPO/convergence_test/logs/diagnose_${TAG}_b1.log"
    julia --project=. convergence_test/diagnose_stage1.jl "$TAG" b1 4 > "$L" 2>&1
    rc=$?
    say "diagnose $TAG rc=$rc"
    note "## B1 diagnostics, $TAG (exit $rc)"
    note ""
    note '```'
    grep -vE "^\[ Info|^┌|^└|^│" "$L" | tail -45 >> "$STATUS"
    note '```'
    note ""
done

# ── the remaining Occam site runs ──────────────────────────────────────────
# 2P is already running or done; do not restart it.
for TAG in 2S 3P 3S; do
    # never clobber a finished run
    if [ -f "$REPO/outputs/occam_${TAG}/verdict_occam_${TAG}.csv" ]; then
        say "occam $TAG already complete, skipping"; continue
    fi
    say "occam site run $TAG"
    rm -rf "$REPO/outputs/occam_${TAG}"
    L="$REPO/occam_test/logs/site_${TAG}.log"
    t0=$(date +%s)
    julia --project=. run_site.jl "configs/occam_${TAG}.toml" > "$L" 2>&1
    rc=$?
    say "occam $TAG rc=$rc minutes=$(( ($(date +%s)-t0)/60 ))"
    note "## Occam site run, $TAG (exit $rc, $(( ($(date +%s)-t0)/60 )) min)"
    note ""
    if [ -f "$REPO/outputs/occam_${TAG}/verdict_occam_${TAG}.csv" ]; then
        note '```'
        grep -E "Verdict|Interval|Admissible|floor|Water criterion" \
            "$REPO/outputs/occam_${TAG}/verdict_occam_${TAG}.csv" 2>/dev/null | cut -c1-160 >> "$STATUS"
        note '```'
    else
        note "No verdict table written -- see occam_test/logs/site_${TAG}.log"
    fi
    note ""
done

say "ALL DONE"
note "---"
note "Queue finished $(date -Is)."
echo "OVERNIGHT_ALL_DONE $(date -Is)"
