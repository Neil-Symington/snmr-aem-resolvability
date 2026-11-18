#!/usr/bin/env bash
# resume_release_smoke.sh -- SIGCONT the paused release-gate smoke once the two
# true-kernel controls have finished. Paused rather than killed on 19 Nov so no
# smoke progress is lost; the controls and the 3S/2S validated gate come first.
set -u
cd "$(dirname "$0")/.." || exit 1
PID="${1:?usage: $0 <pid>}"
echo "$(date '+%H:%M') waiting for the controls before resuming pid $PID"
while ! grep -q "CONTROLS_AND_REPORTS_DONE" logs/ctl_now.log 2>/dev/null; do
    kill -0 "$PID" 2>/dev/null || { echo "$(date '+%H:%M') pid $PID is gone; nothing to resume"; exit 0; }
    sleep 60
done
kill -CONT "$PID" && echo "$(date '+%H:%M') resumed pid $PID"
