#!/usr/bin/env bash
# m4_pair.sh -- run several Stage 1 ensembles of ONE column CONCURRENTLY.
#
#   bash tools/m4_pair.sh [--smoke] <column> <ens> [<ens> ...]
#
# tools/m4_run.sh runs ensembles one at a time: four workers per ensemble on six
# physical cores is what leaves two cores for the gate work. This script is the
# concurrency experiment Neil asked for on 17 Nov -- run two at once and keep it
# only if each ensemble takes less than about 1.7x its single-run time, so that
# total throughput improves.
#
# Stage 0 must already have been done for the column by tools/m4_run.sh; this
# script refuses to run without it, and never writes a Stage 0 artefact of its
# own. Seeds, labels and directory names are IDENTICAL to m4_run.sh's, so an
# ensemble run here is indistinguishable from one run there:
#
#     SEED  = 20261117 + 100*column + ensemble
#     LABEL = _m4c<column>e<ensemble>
#     DIR   = data/aem<TAG><LABEL>
#
# Manifest entries are appended by this one process AFTER every ensemble has
# finished, so there is never a concurrent writer on manifest.toml.
set -u -o pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

SMOKE=0
if [ "${1:-}" = "--smoke" ]; then SMOKE=1; shift; fi
[ "$#" -ge 2 ] || { echo "usage: $0 [--smoke] <column> <ens> [<ens> ...]" >&2; exit 2; }

COL="$1"; shift
ENS="$*"

NSAMPLES=200001; NCHAINS=4; ROOTOUT=data/m4; KIND=""
if [ "$SMOKE" = 1 ]; then NSAMPLES=2001; ROOTOUT=data/m4/smoke; KIND=smoke; fi
MANIFEST="$ROOTOUT/manifest.toml"
STATUS="$REPO/M4_STATUS.md"
[ "$SMOKE" = 1 ] && STATUS="$REPO/$ROOTOUT/STATUS_pair_smoke.md"

status() { echo "- $(date '+%Y-%m-%d %H:%M') $*" | tee -a "$STATUS"; }
die()    { status "PAIR FAILED: $*"; exit 1; }

# kept identical to tools/m4_run.sh; duplicated rather than sourced because that
# script executes a run when sourced.
column_tag() { case "$1" in 1) echo 3P;; 2) echo 2P;; 3) echo 3S;; 4) echo 2S;; *) echo "";; esac; }

check_stage1() {
    local d="$1" f ok=0
    for f in aem_setup.jld2 aem_posterior.jld2 misfits_hypsynth_aem_s.bin \
             models_hypsynth_aem_s.bin points_hypsynth_aem_s.bin; do
        [ -s "$d/$f" ] || { echo "  missing or empty: $d/$f"; ok=1; }
    done
    return $ok
}

TAG=$(column_tag "$COL")
[ -n "$TAG" ] || die "column $COL is not wired"
CDIR="$ROOTOUT/$COL"
S0="$CDIR/stage0_m4_${TAG}${KIND:+_$KIND}.jld2"
[ -s "$S0" ] || die "column $COL: Stage 0 artefact $S0 missing -- run tools/m4_run.sh first"
[ -f "$MANIFEST" ] || die "column $COL: $MANIFEST missing -- run tools/m4_run.sh first"

status "column $COL ($TAG) ${KIND} PAIR start, ensembles: $ENS"
PAIR_T0=$(date +%s)

PIDS=""; LAUNCHED=""
for E in $ENS; do
    SEED=$((20261117 + 100 * COL + E))
    LABEL="_m4${KIND}c${COL}e${E}"
    DIR="data/aem${TAG}${LABEL}"
    LOG="$CDIR/stage1_e${E}.log"
    if [ "$SMOKE" = 0 ] && check_stage1 "$DIR" > /dev/null 2>&1; then
        status "column $COL e$E already complete in $DIR, skipping"
        continue
    fi
    status "column $COL e$E start (concurrent), seed $SEED -> $DIR"
    julia --project=. tools/aem_HYP.jl "$NSAMPLES" "$TAG" "$NCHAINS" "$LABEL" "$SEED" \
        > "$LOG" 2>&1 &
    PIDS="$PIDS $!"
    LAUNCHED="$LAUNCHED $E"
    # per-ensemble start time, for the throughput comparison
    eval "T0_$E=$(date +%s)"
done

[ -n "$LAUNCHED" ] || { status "column $COL: nothing to launch"; exit 0; }

RC=0
for p in $PIDS; do wait "$p" || RC=1; done
PAIR_MINS=$(( ($(date +%s) - PAIR_T0) / 60 ))

# ---- one writer, after every ensemble has finished -----------------------
for E in $LAUNCHED; do
    LABEL="_m4${KIND}c${COL}e${E}"
    DIR="data/aem${TAG}${LABEL}"
    SEED=$((20261117 + 100 * COL + E))
    LOG="$CDIR/stage1_e${E}.log"
    check_stage1 "$DIR" >> "$LOG" || die "column $COL e$E: empty or missing output in $DIR"
    ln -sfn "$REPO/$DIR" "$CDIR/aem_e${E}"
    {
        echo
        echo "[column_$COL.stage1_e$E]"
        echo "seed = $SEED"
        echo "dir = \"$DIR\""
        echo "minutes = $PAIR_MINS   # wall clock of the whole concurrent group"
        echo "concurrency = $(set -- $LAUNCHED; echo $#)"
        echo "group = \"$(echo $LAUNCHED | tr ' ' ',')\""
    } >> "$MANIFEST"
done
[ "$RC" -eq 0 ] || die "at least one ensemble in the group exited non-zero; see $CDIR/stage1_e*.log"
status "column $COL ${KIND} PAIR done: ensembles$LAUNCHED in $PAIR_MINS min wall clock at concurrency $(set -- $LAUNCHED; echo $#)"
