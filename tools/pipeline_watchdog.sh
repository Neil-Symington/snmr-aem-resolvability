#!/usr/bin/env bash
# pipeline_watchdog.sh -- relaunch tools/pipeline_v2c.sh if its heartbeat is stale.
#
# Run from cron every 15 minutes. The pipeline touches its heartbeat every 5,
# so anything older than 15 minutes means it died. Relaunching is safe and
# cheap: every completed step has a marker and is skipped, and a lock stops
# two copies running at once.
#
# Does nothing once logs/chain_state/done_ALL exists.
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
STATE=${WD_STATE:-logs/chain_state_uniform}
HB="$STATE/heartbeat"
LOG=logs/pipeline_watchdog.log
STALE=900          # seconds

[ -f "$STATE/done_ALL" ] && exit 0
# WD_STATE and WD_SCRIPT let one watchdog serve whichever pipeline is live;
# they default to the uniform gate, which is the one running.
say() { echo "- $(date '+%Y-%m-%d %H:%M') $*" >> "$LOG"; }

if [ ! -f "$HB" ]; then
    say "no heartbeat file; starting the pipeline"
elif [ $(( $(date +%s) - $(stat -c %Y "$HB") )) -gt "$STALE" ]; then
    say "heartbeat is $(( ($(date +%s) - $(stat -c %Y "$HB")) / 60 )) min old; relaunching"
else
    exit 0
fi
setsid nohup bash "$REPO/${WD_SCRIPT:-tools/pipeline_uniform.sh}" > /dev/null 2>&1 < /dev/null &
say "relaunched, pid $!"
