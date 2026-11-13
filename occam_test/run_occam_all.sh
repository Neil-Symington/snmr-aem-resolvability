#!/usr/bin/env bash
# Steps 1 and 2 of the Occam brief, on every column: AEM Occam to chi2 <= N,
# then SNMR Occam of water content on that conductivity. Writes
# data/aem_occam<TAG>.jld2 and data/smr_occam<TAG>.jld2 per column.
set -u
cd "$(dirname "$0")/.." || exit 1
mkdir -p occam_test/logs
for TAG in 2P 2S 3P 3S; do
    L="occam_test/logs/occam_${TAG}.log"
    t0=$(date +%s)
    julia --project=. tools/quicklook.jl "$TAG" > "$L" 2>&1
    rc=$?
    echo "OCCAM_DONE tag=$TAG rc=$rc seconds=$(( $(date +%s) - t0 ))" | tee -a "$L"
done
echo "OCCAM_ALL_DONE"
