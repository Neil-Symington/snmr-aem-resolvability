#!/usr/bin/env bash
# run_task.sh -- one unit of a batched field run. Idempotent: a task whose done
# marker exists is skipped, so a resubmitted farm redoes only what is missing.
#
#   bash fieldrun/hpc/run_task.sh <config> aem   <site>          AEM Occam, then AEM trans-D
#   bash fieldrun/hpc/run_task.sh <config> snmr  <site> <e> <r>  one SNMR ladder: kernel draw e, replicate r
#   bash fieldrun/hpc/run_task.sh <config> post  <site>          combine, detect, profile, gate, reports
#   bash fieldrun/hpc/run_task.sh <config> index                 the viewer's index page
#
# Run from the repository root (the Julia project). TASKLOG is the batch's
# bookkeeping folder: <TASKLOG>/tasks/<id>.log, done/<id>, failed/<id>.
#
# WHY "post" IS STEPWISE. A site with no detected interval has no profile and
# no gate, and both stages then exit 1 having done nothing wrong. So within
# post, combine and detect are required; profile, gate and the matplotlib
# report are attempted, their failures recorded, and the viewer page is written
# from whatever exists. Each step has its own marker, so a rerun repeats only
# the steps that did not finish.
set -u
CFG=$1; KIND=$2; SITE=${3:-}; E=${4:-}; R=${5:-}
: "${TASKLOG:?TASKLOG must name the batch bookkeeping folder}"
mkdir -p "$TASKLOG"/tasks "$TASKLOG"/done "$TASKLOG"/failed

case "$KIND" in
    aem|post) ID="${KIND}_${SITE}" ;;
    snmr)     ID="snmr_${SITE}_e${E}r${R}" ;;
    index)    ID="index" ;;
    *) echo "run_task.sh: unknown task kind '$KIND'"; exit 2 ;;
esac
LOG="$TASKLOG/tasks/$ID.log"
if [ -f "$TASKLOG/done/$ID" ] && [ "$KIND" != index ]; then
    echo "  skip   $ID (done)"; exit 0
fi
rm -f "$TASKLOG/failed/$ID"

J() { julia --project=. "fieldrun/stages/$1" --config "$CFG" "${@:2}"; }

# step NAME CMD... : run one step of a task, once, appending to the task log
step() {
    local name=$1; shift
    local mark="$TASKLOG/done/${ID}.${name}"
    [ -f "$mark" ] && { echo "[$(date +%T)] $name already done" >> "$LOG"; return 0; }
    echo "[$(date +%T)] $name: $*" >> "$LOG"
    if "$@" >> "$LOG" 2>&1; then
        : > "$mark"; echo "[$(date +%T)] $name ok" >> "$LOG"; return 0
    else
        echo "[$(date +%T)] $name FAILED" >> "$LOG"; return 1
    fi
}

t0=$(date +%s)
echo "  start  $ID  $(date +%T)"
echo "===== $ID on $(hostname) $(date)  commit $(git log --oneline -1 2>/dev/null)" >> "$LOG"
rc=0
case "$KIND" in
    aem)
        step occam  J run_aem_occam.jl  --site "$SITE" &&
        step transd J run_aem_transd.jl --site "$SITE" || rc=1 ;;
    snmr)
        step unit J run_smr_transd.jl --site "$SITE" --ens "$E" --rep "$R" || rc=1 ;;
    post)
        # required
        step combine J run_smr_transd.jl --site "$SITE" --combine &&
        step detect  J detect_aquifer.jl --site "$SITE" || rc=1
        if [ $rc -eq 0 ]; then
            soft=""
            step profile J run_aem_profile.jl --site "$SITE" &&
            step gate    J gate_water_quality.jl --site "$SITE" || soft="$soft profile/gate"
            step report  J report_smr.jl --site "$SITE" || soft="$soft report"
            step viewer  J report_html.jl --site "$SITE" || rc=1
            [ -n "$soft" ] && echo "  note   $ID: did not complete:$soft (see $LOG)" &&
                echo "incomplete:$soft" > "$TASKLOG/failed/$ID.soft"
        fi ;;
    index)
        J report_html.jl --index >> "$LOG" 2>&1 || rc=1 ;;
esac

dt=$(( $(date +%s) - t0 ))
if [ $rc -eq 0 ]; then
    : > "$TASKLOG/done/$ID"
    echo "  ok     $ID  $(date +%T)  (${dt} s)"
else
    echo "$(date)" > "$TASKLOG/failed/$ID"
    echo "  FAILED $ID  $(date +%T)  (${dt} s) -- see $LOG"
fi
exit $rc
