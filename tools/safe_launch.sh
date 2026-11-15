#!/usr/bin/env bash
# safe_launch.sh -- never start a long run on an unexercised script.
#
#   tools/safe_launch.sh <smoke-log> <real-log> <smoke-cmd> <real-cmd>
#
# 1. Runs <smoke-cmd> to completion, including its save, and waits.
# 2. Refuses to start <real-cmd> if the smoke exited non-zero, or if its log
#    contains a Julia warning, the words "soft scope", or an ERROR.
# 3. Only then launches <real-cmd> with nohup, in the background.
#
# WHY. Four failures in M1 came from Julia's soft scope: a top-level loop that
# rebinds a variable gets a NEW local each iteration. Julia prints
# "Assignment to `x` in soft scope is ambiguous" and carries on, so the run can
# complete and silently lose the state it was accumulating -- which is how a
# warm-start chain was broken without erroring. The warning is the only signal,
# so here it is an error. Two more came from names used before definition in a
# save block, which only executes at the end of a run; a smoke that runs to
# completion is the only thing that catches those.
#
# Written for bash 3.2 (macOS ships it), so no mapfile, no associative arrays.
set -u

if [ "$#" -lt 4 ]; then
    echo "usage: $0 <smoke-log> <real-log> <smoke-cmd> <real-cmd>" >&2
    exit 2
fi
SMOKE_LOG="$1"; REAL_LOG="$2"; SMOKE_CMD="$3"; REAL_CMD="$4"

mkdir -p "$(dirname "$SMOKE_LOG")" "$(dirname "$REAL_LOG")"

echo "=== smoke ==="
echo "  $SMOKE_CMD"
set +e
eval "$SMOKE_CMD" > "$SMOKE_LOG" 2>&1
rc=$?
set -e

fail() { echo; echo "REFUSING TO LAUNCH: $1"; echo "  see $SMOKE_LOG"; exit 1; }

[ "$rc" -eq 0 ] || fail "the smoke exited $rc"

# Julia's soft-scope warning, in either the boxed or the plain form.
if grep -qi "soft scope" "$SMOKE_LOG"; then
    echo; grep -i -B1 -A1 "soft scope" "$SMOKE_LOG" | head -20
    fail "the smoke log contains a soft-scope warning -- treat it as an error"
fi
# Any Julia warning. Boxed warnings start with the box-drawing character.
if grep -qE "^┌ Warning:|^Warning:|WARNING:" "$SMOKE_LOG"; then
    echo; grep -nE "^┌ Warning:|^Warning:|WARNING:" "$SMOKE_LOG" | head -10
    fail "the smoke log contains a Julia warning"
fi
if grep -qE "^ERROR|LoadError|UndefVarError|MethodError|BoundsError" "$SMOKE_LOG"; then
    echo; grep -nE "^ERROR|LoadError|UndefVarError|MethodError|BoundsError" "$SMOKE_LOG" | head -10
    fail "the smoke log contains an error"
fi

echo "  smoke clean: exit 0, no warnings, no soft scope, no errors"
echo
echo "=== launching the long run ==="
echo "  $REAL_CMD"
nohup bash -c "$REAL_CMD" > "$REAL_LOG" 2>&1 &
echo "  pid $!  log $REAL_LOG"
