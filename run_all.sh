#!/usr/bin/env bash
# run_all.sh -- the whole synthetic workflow, from a scenario config to a
# verdict table. Edit the block below and run it.
#
#   ./run_all.sh                 # smoke: one column end to end, ~1 h
#   ./run_all.sh --force         # redo stages whose output already exists
#   ./run_all.sh --dry-run       # print what it would run, run nothing
#
# WHY THIS EXISTS. run_site.jl drives Stages 0 and 3-6 but deliberately does
# NOT run Stages 1 and 2: each is a trans-dimensional MCMC of hours, and the
# driver printed the commands rather than starting a job that outlives your
# terminal. That left seven commands to copy by hand, in order, with three
# traps in them, and the field half of the repository had fieldrun/run_batch.sh
# while the synthetic half had nothing. This is the synthetic equivalent.
#
# IT IS NOT A NEW SOURCE OF TRUTH. Every scientific choice still lives in the
# site config under configs/ and the scenario config under synthetics/config/.
# What is set here is only what used to be argv: which configs, which tag, how
# many samples. If a number matters to a result, it belongs in the config, not
# in this file.
#
# SMOKE IS NOT A SMALLER ANSWER, IT IS A DIFFERENT ONE. At 2001 samples the
# hierarchical noise scale runs far above its production value and the
# posterior is not the same shape. Use smoke to prove the paths run. Never
# quote a number from it.
set -u
cd "$(dirname "$0")"

# ─── EDIT HERE ───────────────────────────────────────────────────────────────

MODE="smoke"          # "smoke" (minutes to an hour) or "production" (a day)

# The scenario, built by synthetics/scripts/build_truth.jl.
TRUTH_CONFIG="synthetics/config/confined_aquifer.toml"

# The site config run_site.jl drives Stages 0 and 3-6 with. smoke_2P.toml
# differs from synth_2P.toml only in its stage1_draws tag and a reduced gate
# grid; the grid is column-specific, so a different column needs a different
# one (see docs/columns.md).
SITE_CONFIG="configs/smoke_2P.toml"

TAG="2P"              # names every intermediate: data/truth_murray<TAG>.jld2
LABEL="_smoke"        # Stage 1 writes data/aem<TAG><LABEL>/; "" for none
TIER="typical"        # optimistic | typical | sceptical | field
CONDS="C2"            # C2 is the corrected log prior on the Stage 1 median

# Drop a name to skip that stage. Order is fixed and enforced below.
STAGES="truth stage1 draws smrdata stage2 site"

#                     samples   AEM chains   AEM draws
SMOKE_NSAMPLES=2001;  SMOKE_CHAINS=2;  SMOKE_DRAWS=50
PROD_NSAMPLES=200001; PROD_CHAINS=4;   PROD_DRAWS=600

# ─── STOP EDITING ────────────────────────────────────────────────────────────

FORCE=0; DRY=0
while [ $# -gt 0 ]; do
    case "$1" in
        --force)   FORCE=1; shift ;;
        --dry-run) DRY=1; shift ;;
        *) echo "unknown argument: $1"; echo "usage: $0 [--force] [--dry-run]"; exit 2 ;;
    esac
done

case "$MODE" in
    smoke)      NSAMPLES=$SMOKE_NSAMPLES; CHAINS=$SMOKE_CHAINS; DRAWS=$SMOKE_DRAWS ;;
    production) NSAMPLES=$PROD_NSAMPLES;  CHAINS=$PROD_CHAINS;  DRAWS=$PROD_DRAWS ;;
    *) echo "MODE must be smoke or production, not '$MODE'"; exit 2 ;;
esac

for f in "$TRUTH_CONFIG" "$SITE_CONFIG"; do
    [ -f "$f" ] || { echo "no config at $f"; exit 1; }
done

# The Julia version is the commonest setup failure and the error it causes
# names something else entirely, so check it before doing an hour of work.
# install.jl enforces the same range; this is the same check, early.
julia -e 'VERSION >= v"1.12" && VERSION < v"1.13" ||
          (println("ERROR: this project needs Julia 1.12.x; you are on $VERSION.\n" *
                   "  juliaup add 1.12 && juliaup default 1.12"); exit(1))' || exit 1

# Read the scenario name from the config with Julia, the same way run_batch.sh
# does, so the shell and the TOML cannot disagree about what gets built.
TRUTH_NAME=$(julia -e 'import TOML; print(TOML.parsefile(ARGS[1])["name"])' "$TRUTH_CONFIG") \
    || { echo "could not read [name] from $TRUTH_CONFIG"; exit 1; }
TRUTH_SRC="synthetics/output/${TRUTH_NAME}_truth.jld2"
TRUTH_DST="data/truth_murray${TAG}.jld2"

LOGDIR="logs/run_all_${TAG}${LABEL}"
mkdir -p "$LOGDIR" data

# What each stage writes, so a re-run can tell done from not-started.
out_of() {
    case "$1" in
        truth)   echo "$TRUTH_DST" ;;
        stage1)  echo "data/aem${TAG}${LABEL}/aem_posterior.jld2" ;;
        draws)   echo "data/aem_draws_${TAG}.jld2" ;;
        smrdata) echo "data/smr_data_murray${TAG}T.jld2" ;;
        stage2)  echo "data/cond_${CONDS%%,*}_${TIER}${TAG}T.jld2" ;;
        site)    echo "" ;;      # run_site.jl is itself resumable, via --from
    esac
}

# Rough cost, printed so a stage that looks hung can be told from one that is.
cost_of() {
    case "$1$MODE" in
        truthsmoke|truthproduction)   echo "seconds" ;;
        stage1smoke)                  echo "~10 min" ;;
        stage1production)             echo "~3.5 h" ;;
        drawssmoke|drawsproduction)   echo "minutes" ;;
        smrdatasmoke|smrdataproduction) echo "minutes" ;;
        stage2smoke)                  echo "~10 min" ;;
        stage2production)             echo "hours" ;;
        sitesmoke)                    echo "~45 min, nearly all Stage 4" ;;
        siteproduction)               echo "~1.5 h, nearly all Stage 4" ;;
    esac
}

declare -a TIMED=()
FAILED=""

run_stage() {
    local name="$1"; shift
    local out; out=$(out_of "$name")
    local log="$LOGDIR/$name.log"

    if [ -n "$out" ] && [ -e "$out" ] && [ "$FORCE" -eq 0 ]; then
        echo "  SKIP  $name -- $out exists (--force to redo)"
        TIMED+=("$name skipped")
        return 0
    fi
    echo "  RUN   $name  ($(cost_of "$name"))  -> $log"
    if [ "$DRY" -eq 1 ]; then
        echo "        $*"
        TIMED+=("$name dry-run")
        return 0
    fi
    local t0; t0=$(date +%s)
    if "$@" > "$log" 2>&1; then
        local dt=$(( $(date +%s) - t0 ))
        printf "        done in %dm %02ds\n" $((dt/60)) $((dt%60))
        TIMED+=("$name $((dt/60))m$((dt%60))s")
        if [ -n "$out" ] && [ ! -e "$out" ]; then
            # Exiting 0 without writing the output is the failure mode that let
            # a whole field run print "complete" having produced nothing.
            echo "  !!!   $name exited 0 but did not write $out"
            FAILED="$FAILED $name"
            return 1
        fi
    else
        echo "  !!!   $name FAILED -- see $log"
        tail -5 "$log" | sed 's/^/        /'
        FAILED="$FAILED $name"
        return 1
    fi
}

wants() { case " $STAGES " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

for s in $STAGES; do
    case "$s" in truth|stage1|draws|smrdata|stage2|site) ;;
        *) echo "unknown stage '$s' in STAGES"; exit 2 ;;
    esac
done

echo "run_all.sh  mode=$MODE  tag=$TAG  label=${LABEL:-none}  tier=$TIER  conds=$CONDS"
echo "  scenario $TRUTH_CONFIG -> $TRUTH_NAME"
echo "  site     $SITE_CONFIG"
echo "  samples  $NSAMPLES   AEM chains $CHAINS   AEM draws $DRAWS"
echo "  logs     $LOGDIR"
# Stage 2's chain and ensemble counts are NOT settable: stage2_snmr.jl holds
# both at 4 as constants and reads no flag for either, so nothing here pretends
# to set them. --nchains, which README and run_site.jl both still print, is
# silently ignored by that stage.
echo "  Stage 2 runs 4 chains x 4 ensembles regardless; that stage takes no flag for it"
[ "$MODE" = smoke ] && echo "  SMOKE: proves the paths run. Do not quote a number from it."
echo

if wants truth; then
    run_stage truth bash -c "julia --project=. synthetics/scripts/build_truth.jl '$TRUTH_CONFIG' && \
                             cp '$TRUTH_SRC' '$TRUTH_DST'" || true
fi
# Stage 1 takes POSITIONAL arguments, not flags: parse(Int, "--nsamples") is
# what you get for passing them, and a reviewer did.
wants stage1  && [ -z "$FAILED" ] && \
    run_stage stage1 julia --project=. tools/aem_HYP.jl "$NSAMPLES" "$TAG" "$CHAINS" "$LABEL" || true
wants draws   && [ -z "$FAILED" ] && \
    run_stage draws julia --project=. tools/aem_draws.jl "$TAG" "$DRAWS" "$LABEL" || true
wants smrdata && [ -z "$FAILED" ] && \
    run_stage smrdata julia --project=. tools/smr_data.jl transd --tag "$TAG" --label "$LABEL" || true
# The T suffix is not a typo: smr_data.jl's trans-D branch writes it so the
# Occam and trans-D kernels cannot collide, and Stage 2 must be told the same.
wants stage2  && [ -z "$FAILED" ] && \
    run_stage stage2 julia --project=. stages/stage2_snmr.jl \
        --tier "$TIER" --conds "$CONDS" --tag "${TAG}T" --label "$LABEL" \
        --nsamples "$NSAMPLES" --ndraws "$DRAWS" || true

# The trap that cost a reviewer an hour: a site config can name a different
# Stage 1 tag than the one just built (synth_2P.toml names 2M, synth_3P.toml
# names 3M, because the paper's Stage 1 ran under the older tag and sigma_bulk
# is identical). run_site.jl warns early but carries on, and dies in Stage 5 an
# hour later. Copy it here instead, and say so.
if wants site && [ -z "$FAILED" ]; then
    WANT_DRAWS=$(julia -e 'import TOML; print(get(TOML.parsefile(ARGS[1])["data"], "stage1_draws", ""))' "$SITE_CONFIG")
    if [ -n "$WANT_DRAWS" ] && [ ! -e "$WANT_DRAWS" ] && [ -e "data/aem_draws_${TAG}.jld2" ]; then
        echo "  NOTE  $SITE_CONFIG wants $WANT_DRAWS; copying data/aem_draws_${TAG}.jld2 to it."
        echo "        Same sigma_bulk under a different tag. Stage 5 would die an hour in without this."
        [ "$DRY" -eq 1 ] || cp "data/aem_draws_${TAG}.jld2" "$WANT_DRAWS"
    fi
    run_stage site julia --project=. run_site.jl "$SITE_CONFIG" || true
fi

echo
printf '%s\n' "stage timings:" "${TIMED[@]/#/  }"
if [ -n "$FAILED" ]; then
    echo "INCOMPLETE -- failed:$FAILED"
    exit 1
fi
if [ "$DRY" -eq 1 ]; then
    echo "dry run only, nothing was executed"
    exit 0
fi
if wants site; then
    V=$(julia -e 'import TOML; c = TOML.parsefile(ARGS[1])
                  print(c["output"]["dir"], "/verdict_", c["name"], ".md")' "$SITE_CONFIG")
    if [ -e "$V" ]; then
        echo "complete -- verdict table at $V"
    else
        # "the loop ended" is not "it worked". Saying complete over a run that
        # produced no output is the defect this line exists to prevent.
        echo "INCOMPLETE -- no verdict table at $V; check $LOGDIR/site.log"
        exit 1
    fi
else
    echo "complete -- stages run: $STAGES"
fi
