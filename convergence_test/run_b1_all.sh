#!/usr/bin/env bash
# B1 baseline: production settings, four independent ensembles, both columns.
# This doubles as the Part A data, because the production chain files were
# never on this machine -- see the Part A report.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"

# wait for any earlier Stage 1 work to drain before taking the whole machine
while pgrep -f "tools/aem_HYP.jl" > /dev/null 2>&1; do sleep 20; done
echo "machine clear at $(date -Is)"

for TAG in 2P 3P; do
    echo "### B1 $TAG starting $(date -Is)"
    bash "$HERE/run_stage1_setting.sh" b1 "$TAG" 200001 4
    echo "### B1 $TAG finished $(date -Is)"
done
echo "B1_ALL_DONE $(date -Is)"
