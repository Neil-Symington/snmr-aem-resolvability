#!/usr/bin/env bash
# m4_run.sh -- M4 inversions: Stage 0 and Stage 1 for the numbered columns.
#
#   bash tools/m4_run.sh [--smoke] <column> [<column> ...]
#
# Columns (DECISIONS.md, 17 Nov):
#   1 = 3P  test one, fresh      configs/synth_3P.toml
#   2 = 2P  test two, fresh      configs/synth_2P.toml
#   3 = 3S  test one, brackish   configs/synth_3S.toml
#   4 = 2S  test two, brackish   configs/synth_2S.toml
#   5 = 2Pv2 test two v2, fresh    configs/synth_2Pv2.toml
#   6 = 2Sv2 test two v2, brackish configs/synth_2Sv2.toml
#
# Columns 5 and 6 wired 21 Nov. They are test two rebuilt with TEST ONE'S
# [variability] -- jitter logit, corr_length 1.0, log_sd_phi 0.35,
# log_sd_clay 0.45, log_sd_sigma_w 0.25 -- so the two grounds differ in
# geometry and confinement and no longer also in roughness. Same units, depths,
# water and window as v1. Their windows are v1's, 17.5-33.5 m.
#
# Columns 3 and 4 wired 17 Nov on the second machine, their truth columns having
# been built from paper_3P_brackish.toml and confined_aquifer_brackish.toml.
# Their windows and packages are those of the fresh column of the same pair,
# because the pair differs only in the aquifer's pore water.
#
# Calls the existing stage scripts with arguments and environment only; no
# stage code is changed. Stage 1 runs four independent ensembles per column,
# one at a time (four physical cores, four workers per ensemble), at production
# settings: 200,001 samples, 4 chains, Tmax 2.5, layer limit 40, one chain at
# T = 1, burn-in 0.5. Seed for ensemble e of column c: 20261117 + 100c + e.
#
# Stage 2 is NOT run here. That code change has since been made: as of
# m4-stage2-pool-20261117, stages/stage2_snmr.jl takes --pool with a list of
# Stage 1 folders and uses the median of the POOLED post-burn-in T=1 draws.
# Pooling one folder reproduces that folder's sigma50 exactly.
#
# Outputs, all new:
#   data/m4/<c>/site_m4_<TAG>.toml        Stage 0 site config (extends synth_<TAG>)
#   data/m4/<c>/stage0_m4_<TAG>.jld2      Stage 0 sealed artefact
#   data/aem<TAG>_m4c<c>e<N>/             Stage 1 ensemble N (aem_HYP.jl fixes the
#                                         aem<TAG><LABEL> prefix)
#   data/m4/<c>/aem_e<N> -> the above     symlink, for the agreed layout
#   data/m4/manifest.toml                 seeds, hashes, digests, timings
#   M4_STATUS.md                          one line per event
# --smoke writes under data/m4/smoke/ and data/aem<TAG>_m4smokec<c>e1/ instead,
# with one ensemble of 2,001 samples, and still exercises every save.
set -u -o pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO" || exit 1
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

SMOKE=0
if [ "${1:-}" = "--smoke" ]; then SMOKE=1; shift; fi
[ "$#" -ge 1 ] || { echo "usage: $0 [--smoke] <column> [<column> ...]" >&2; exit 2; }

NSAMPLES=200001; NENS=4; NCHAINS=4
ROOTOUT=data/m4; KIND=""
if [ "$SMOKE" = 1 ]; then NSAMPLES=2001; NENS=1; ROOTOUT=data/m4/smoke; KIND=smoke; fi
MANIFEST="$ROOTOUT/manifest.toml"
STATUS="$REPO/M4_STATUS.md"
[ "$SMOKE" = 1 ] && STATUS="$REPO/$ROOTOUT/STATUS_smoke.md"
mkdir -p "$ROOTOUT"

status() { echo "- $(date '+%Y-%m-%d %H:%M') $*" | tee -a "$STATUS"; }
die()    { status "FAILED: $*"; exit 1; }

column_tag() { case "$1" in 1) echo 3P;; 2) echo 2P;; 3) echo 3S;; 4) echo 2S;;
                            5) echo 2Pv2;; 6) echo 2Sv2;; *) echo "";; esac; }
column_win() { case "$1" in 1) echo "7.5 23.5 6.0 40.0";; 2) echo "17.5 33.5 4.0 34.0";;
                                3) echo "7.5 23.5 6.0 40.0";; 4) echo "17.5 33.5 4.0 34.0";;
                                5) echo "17.5 33.5 4.0 34.0";; 6) echo "17.5 33.5 4.0 34.0";; esac; }

# every file in a Stage 1 directory must exist and be non-empty
check_stage1() {
    local d="$1" f ok=0
    for f in aem_setup.jld2 aem_posterior.jld2 misfits_hypsynth_aem_s.bin \
             models_hypsynth_aem_s.bin points_hypsynth_aem_s.bin; do
        [ -s "$d/$f" ] || { echo "  missing or empty: $d/$f"; ok=1; }
    done
    return $ok
}

if [ ! -f "$MANIFEST" ]; then
    {
        echo "# M4 inversion manifest. Written by tools/m4_run.sh; appended per step."
        echo "[run]"
        echo "commit       = \"$(git rev-parse HEAD)\""
        echo "runner       = \"tools/m4_run.sh\""
        echo "runner_sha256 = \"$(sha256sum tools/m4_run.sh | cut -d' ' -f1)\""
        echo "digest_sha256 = \"$(sha256sum tools/m4_digest.jl | cut -d' ' -f1)\""
        echo "host         = \"$(hostname)\""
        echo "physical_cores = $(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l)"
        echo "concurrency  = 1"
        echo "julia        = \"$(julia --version)\""
        echo "nsamples = $NSAMPLES"
        echo "nchains  = $NCHAINS"
        echo "tmax = 2.5"
        echo "nmax_k = 40"
        echo "nchainsatone = 1"
        echo "burninfrac = 0.5"
        echo "env = \"JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1\""
    } > "$MANIFEST"
fi

for COL in "$@"; do
    TAG=$(column_tag "$COL")
    [ -n "$TAG" ] || die "column $COL is not wired in this runner"
    TRUTH="data/truth_murray${TAG}.jld2"
    [ -s "$TRUTH" ] || die "column $COL: $TRUTH missing"
    CDIR="$ROOTOUT/$COL"; mkdir -p "$CDIR"
    status "column $COL ($TAG) ${KIND} start"

    # ── identity of the truth ──────────────────────────────────────────────
    {
        echo
        echo "[column_$COL]"
        echo "tag = \"$TAG\""
        julia --project=. tools/m4_digest.jl "$TRUTH" $(column_win "$COL") || exit 1
    } >> "$MANIFEST" || die "column $COL: digest"

    # ── Stage 0 into the new directory ─────────────────────────────────────
    SITE="$CDIR/site_m4_${TAG}.toml"
    cat > "$SITE" <<EOF
# Stage 0 site config for M4 column $COL, written by tools/m4_run.sh.
# Everything comes from configs/synth_${TAG}.toml except the name, the output
# dir and snmr_ensemble. M4's Stage 2 has not run, so no ensemble is named:
# the parent's path is an earlier run's file (and absent for 3P).
extends = "../../../$( [ "$SMOKE" = 1 ] && echo "../" )configs/synth_${TAG}.toml"
name  = "m4_${TAG}${KIND:+_$KIND}"

[data]
snmr_ensemble = ""

[output]
dir = "$CDIR"
EOF
    S0LOG="$CDIR/stage0.log"
    julia --project=. stages/stage0_inputs.jl --config "$SITE" > "$S0LOG" 2>&1 \
        || die "column $COL: Stage 0 exited non-zero, see $S0LOG"
    S0="$CDIR/stage0_m4_${TAG}${KIND:+_$KIND}.jld2"
    [ -s "$S0" ] || die "column $COL: Stage 0 artefact missing or empty ($S0)"
    echo "stage0 = \"$S0\"" >> "$MANIFEST"
    status "column $COL Stage 0 done: $S0"

    # M4_STAGE0_ONLY=1 stops after Stage 0, so the ensembles can be run
    # concurrently by tools/m4_pair.sh instead. Added 17 Nov for the
    # concurrency experiment; the default path is unchanged.
    if [ "${M4_STAGE0_ONLY:-0}" = 1 ]; then
        status "column $COL ($TAG) Stage 0 only, as asked; ensembles not started here"
        continue
    fi

    # ── Stage 1, four independent ensembles, one at a time ────────────────
    for E in $(seq 1 "$NENS"); do
        SEED=$((20261117 + 100 * COL + E))
        LABEL="_m4${KIND}c${COL}e${E}"
        DIR="data/aem${TAG}${LABEL}"
        LOG="$CDIR/stage1_e${E}.log"
        if [ "$SMOKE" = 0 ] && check_stage1 "$DIR" > /dev/null 2>&1; then
            status "column $COL e$E already complete in $DIR, skipping"
            continue
        fi
        status "column $COL e$E Stage 1 start, seed $SEED -> $DIR"
        t0=$(date +%s)
        julia --project=. tools/aem_HYP.jl "$NSAMPLES" "$TAG" "$NCHAINS" "$LABEL" "$SEED" > "$LOG" 2>&1
        rc=$?
        mins=$(( ($(date +%s) - t0) / 60 ))
        [ "$rc" -eq 0 ] || die "column $COL e$E: aem_HYP.jl exited $rc after $mins min, see $LOG"
        check_stage1 "$DIR" >> "$LOG" || die "column $COL e$E: empty or missing output in $DIR"
        ln -sfn "$REPO/$DIR" "$CDIR/aem_e${E}"
        {
            echo
            echo "[column_$COL.stage1_e$E]"
            echo "seed = $SEED"
            echo "dir = \"$DIR\""
            echo "minutes = $mins"
        } >> "$MANIFEST"
        status "column $COL e$E Stage 1 done in $mins min, seed $SEED"
    done
    status "column $COL ($TAG) Stage 1 complete. Stage 2 not run: it cannot pool four ensembles."
done
status "m4_run.sh ${KIND} finished: columns $*"
