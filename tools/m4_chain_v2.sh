#!/usr/bin/env bash
# m4_chain_v2.sh -- the remaining Stage 1 pairs for columns 5 and 6, in order.
#
# Written 22 Nov 2026 after a session-bound waiter died silently and left the
# machine idle for 22 hours: pair 1 of 4 completed at 09:40 on 21 Nov and
# nothing followed. This script is launched detached (setsid + nohup) so it
# outlives the session that started it, and it appends to its own log after
# every pair, so progress is visible from the filesystem alone.
set -u
cd "$(dirname "$0")/.." || exit 1
LOG="logs/m4_chain_v2.log"
say() { echo "- $(date '+%Y-%m-%d %H:%M') $*" >> "$LOG"; }

say "chain start, pairs remaining: c5 e3e4, c6 e1e2, c6 e3e4"
for spec in "5 3 4" "6 1 2" "6 3 4"; do
    set -- $spec
    col=$1; shift
    tag="c${col}_e$(echo "$@" | tr -d ' ')"
    say "pair $tag start"
    if bash tools/m4_pair.sh "$col" "$@" > "logs/m4_pair_${tag}.log" 2>&1; then
        say "pair $tag DONE"
    else
        say "pair $tag FAILED (exit $?); chain stops here"
        exit 1
    fi
done
say "CHAIN COMPLETE: Stage 1 finished for 2Pv2 and 2Sv2"
