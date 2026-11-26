#!/usr/bin/env bash
# farm.sh -- run a task list, SLOTS at a time, on one machine (a Gadi node or
# the Mac). One line per task: the arguments run_task.sh takes after <config>.
#
#   TASKLOG=<folder> bash fieldrun/hpc/farm.sh <config> <tasklist> <slots>
#
# Always exits 0 once the list is worked through: a failed site must not hold
# back the other sites' next phase (the PBS dependencies are afterany). Each
# downstream task checks its own inputs, and the failed/ markers say what to
# resubmit.
set -u
CFG=$1; LIST=$2; SLOTS=$3
: "${TASKLOG:?TASKLOG must name the batch bookkeeping folder}"
cd "$(dirname "$0")/../.."                     # the repository root
export TASKLOG
n=$(grep -c . "$LIST" || true)
echo "farm: $n tasks from $(basename "$LIST"), $SLOTS at a time, $(date)"

# One BLAS thread per process: the parallelism is the chains and the slots.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1

# Load the packages once before the slots start, so N processes do not all
# try to precompile at the same moment on a fresh node or a changed checkout.
julia --project=. -e 'using HiQGA, SMRPInversion, SNMRForward, JLD2' ||
    echo "farm: WARNING package load failed; tasks will report it"

grep . "$LIST" | xargs -P "$SLOTS" -L 1 bash fieldrun/hpc/run_task.sh "$CFG" || true

nd=0; nf=0
while read -r kind site e r; do
    case $kind in snmr) id="snmr_${site}_e${e}r${r}" ;; index) id=index ;; *) id="${kind}_${site}" ;; esac
    [ -f "$TASKLOG/done/$id" ] && nd=$((nd + 1))
    [ -f "$TASKLOG/failed/$id" ] && nf=$((nf + 1))
done < <(grep . "$LIST")
echo "farm: done $nd, failed $nf, of $n  $(date)"
exit 0
