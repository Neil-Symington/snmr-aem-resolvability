#!/usr/bin/env bash
# Bundle the SMALL diagnostic files from the 2P production run so they can be
# analysed on another machine. Reads only; writes one tarball to /tmp.
#
#   bash convergence_test/collect_production.sh
#
# The misfit files are plain text and hold chain, iter, the five acceptance
# rates, node count, misfit and temperature -- everything Part A needs except
# the model samples. The cond file holds Stage 2's four ensembles separately.
# The big points_/models_ binaries are deliberately NOT included.
set -u
cd "$(dirname "$0")/.." || exit 1
OUT=/tmp/production_2P_diagnostics.tar.gz
LIST=/tmp/_collect_list.$$
: > "$LIST"

for d in data/aem2P data/aem2PT data/aem2; do
    [ -d "$d" ] && find "$d" -maxdepth 1 -name "misfits_*" -o -maxdepth 1 -name "aem_setup.jld2" \
        -o -maxdepth 1 -name "aem_posterior.jld2" 2>/dev/null | grep "^$d" >> "$LIST"
done
[ -d data/smr_runs ] && find data/smr_runs -maxdepth 1 -name "misfits_*" >> "$LIST" 2>/dev/null
find data -maxdepth 1 -name "cond_C2_*.jld2" >> "$LIST" 2>/dev/null
[ -d convergence_test/results ] && find convergence_test/results -type f >> "$LIST" 2>/dev/null

sort -u "$LIST" -o "$LIST"
n=$(wc -l < "$LIST")
if [ "$n" -eq 0 ]; then
    echo "found nothing to collect under data/ -- is this the machine with the production run?"
    rm -f "$LIST"; exit 1
fi
echo "collecting $n files:"
while read -r f; do printf "  %10s  %s\n" "$(du -h "$f" | cut -f1)" "$f"; done < "$LIST"
tar czf "$OUT" -T "$LIST"
rm -f "$LIST"
echo
echo "wrote $OUT  ($(du -h "$OUT" | cut -f1))"
echo "send that file back."
