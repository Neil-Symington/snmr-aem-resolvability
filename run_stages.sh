#!/usr/bin/env bash
# run_stages.sh -- the engine behind run_smoke.sh and run_production.sh.
#
# NOT AN ENTRY POINT. Run ./run_smoke.sh or ./run_production.sh; each sets its
# own parameters and then sources this file. There are two of them, rather than
# one script with a mode switch, because a mode switch set the sample counts and
# left the label and the site config behind: five stages then skipped in silence
# and a smoke answer came out wearing a production label.
#
# What it does: runs the stages named in STAGES, in a fixed order, one log per
# stage under logs/, a timing for each, and a non-zero exit naming the stage on
# any failure OR on a stage that exits 0 without writing its output.
#
# Flags, handled here for both callers:
#   --dry-run       print the commands, run nothing
#   --force         redo every stage whatever is on disk
#   --reuse-stale   keep an output built with other settings (says so loudly)
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"      # the repository root: the Julia project

for v in RUN_KIND TAG TRUTH_CONFIG SITE_CONFIG LABEL TIER CONDS STAGES \
         NSAMPLES CHAINS DRAWS; do
    eval ": \${$v?run_stages.sh is sourced by run_smoke.sh or run_production.sh, \
which set this. Run one of those, not this file.}" || exit 2
done
case "$RUN_KIND" in smoke|production) ;;
    *) echo "RUN_KIND must be smoke or production, not '$RUN_KIND'"; exit 2 ;;
esac

FORCE=0; DRY=0; REUSE_STALE=0
while [ $# -gt 0 ]; do
    case "$1" in
        --force)        FORCE=1; shift ;;
        --dry-run)      DRY=1; shift ;;
        --reuse-stale)  REUSE_STALE=1; shift ;;
        *) echo "unknown argument: $1"
           echo "usage: $0 [--dry-run] [--force] [--reuse-stale]"; exit 2 ;;
    esac
done

for f in "$TRUTH_CONFIG" "$SITE_CONFIG"; do
    [ -f "$f" ] || { echo "no config at $f"; exit 1; }
done

# Each script is internally consistent as shipped, so this only fires when one
# has been edited into a contradiction. It is cheap and it catches the mistake
# that produced a smoke verdict table labelled as production.
case "$RUN_KIND" in
  production)
    wrong=""
    case "$LABEL" in *smoke*) wrong="$wrong
    LABEL is \"$LABEL\"  ->  set LABEL=\"\" or LABEL=\"_prod\"" ;; esac
    case "$(basename "$SITE_CONFIG")" in smoke*) wrong="$wrong
    SITE_CONFIG is $SITE_CONFIG  ->  set SITE_CONFIG=\"configs/synth_${TAG}.toml\"
                 (a smoke config also carries the REDUCED gate grid, 4+3 nodes
                  against production's 16+15, so Stage 4 answers a smaller
                  problem)" ;; esac
    if [ -n "$wrong" ]; then
        cat <<EOF
REFUSING: this is run_production.sh but smoke settings are in the block above.
$wrong

  Those outputs share filenames with the smoke ones, so the run would have
  reused smoke files and produced a smoke answer under a production label.
  Fix the lines above, then --dry-run to see every stage say RUN.
EOF
        exit 2
    fi ;;
  smoke)
    case "$(basename "$SITE_CONFIG")" in
      smoke*) ;;
      *) cat <<EOF
NOTE: $SITE_CONFIG is not a smoke config, so it carries the PRODUCTION gate grid
      of 16+15 nodes. Stage 4 will then take hours rather than about 45 minutes:
      its cost is an Occam solve plus a sweep plus a bisection PER NODE, and the
      sample count does nothing to it. Only 2P ships a smoke config; copy the
      [gate] block out of configs/smoke_2P.toml for another column.

EOF
      ;;
    esac ;;
esac

# The Julia version is the commonest setup failure and the error it causes names
# something else entirely, so check it before doing an hour of work. install.jl
# enforces the same range; this is the same check, early.
julia -e 'VERSION >= v"1.12" && VERSION < v"1.13" ||
          (println("ERROR: this project needs Julia 1.12.x; you are on $VERSION.\n" *
                   "  juliaup add 1.12 && juliaup default 1.12"); exit(1))' || exit 1

# Read the scenario name from the config with Julia, so the shell and the TOML
# cannot disagree about what gets built.
TRUTH_NAME=$(julia -e 'import TOML; print(TOML.parsefile(ARGS[1])["name"])' "$TRUTH_CONFIG") \
    || { echo "could not read [name] from $TRUTH_CONFIG"; exit 1; }
TRUTH_SRC="synthetics/output/${TRUTH_NAME}_truth.jld2"
TRUTH_DST="data/truth_murray${TAG}.jld2"

LOGDIR="logs/run_${RUN_KIND}_${TAG}${LABEL}"
mkdir -p "$LOGDIR" data

out_of() {
    case "$1" in
        truth)     echo "$TRUTH_DST" ;;
        stage1)    echo "data/aem${TAG}${LABEL}/aem_posterior.jld2" ;;
        draws)     echo "data/aem_draws_${TAG}.jld2" ;;
        smrdata)   echo "data/smr_data_murray${TAG}T.jld2" ;;
        stage2)    echo "data/cond_${CONDS%%,*}_${TIER}${TAG}T.jld2" ;;
        quicklook) echo "" ;;    # cheap enough to always redo, and it is a check
        site)      echo "" ;;    # run_site.jl is itself resumable, via --from
    esac
}

# Rough cost, printed so a stage that looks hung can be told from one that is.
cost_of() {
    case "$1.$RUN_KIND" in
        truth.*)              echo "seconds" ;;
        quicklook.*)          echo "~2 min, no MCMC" ;;
        draws.*|smrdata.*)    echo "minutes" ;;
        stage1.smoke)         echo "~10 min" ;;
        stage1.production)    echo "~3.5 h" ;;
        stage2.smoke)         echo "~10 min" ;;
        stage2.production)    echo "hours" ;;
        site.smoke)           echo "~45 min, nearly all Stage 4" ;;
        site.production)      echo "~1.5 h, nearly all Stage 4" ;;
    esac
}

# ── PROVENANCE, which is what makes a skip safe ─────────────────────────────
# aem_draws_<TAG>.jld2, smr_data_murray<TAG>T.jld2 and cond_*_<TIER><TAG>T.jld2
# carry neither the run kind nor the label, because the tools that write them
# take a label and do not put it in the filename. So a smoke file and a
# production file are the same path, and existence alone cannot tell them apart.
# Skipping on existence reused smoke outputs in a production run, silently --
# the same shape as every other defect this project has recorded: something that
# looks like it governs a stage and does not. This ledger records the settings
# each output was built under, and a stage is skipped ONLY on an exact match.
PROV="data/.run_provenance"

sig_of() {
    case "$1" in
        truth)   printf 'truth|%s' "$TRUTH_CONFIG" ;;
        stage1)  printf 'kind=%s|label=%s|nsamples=%s|chains=%s' \
                        "$RUN_KIND" "$LABEL" "$NSAMPLES" "$CHAINS" ;;
        draws)   printf 'kind=%s|label=%s|ndraws=%s' "$RUN_KIND" "$LABEL" "$DRAWS" ;;
        smrdata) printf 'kind=%s|label=%s' "$RUN_KIND" "$LABEL" ;;
        stage2)  printf 'kind=%s|label=%s|tier=%s|conds=%s|nsamples=%s|ndraws=%s' \
                        "$RUN_KIND" "$LABEL" "$TIER" "$CONDS" "$NSAMPLES" "$DRAWS" ;;
        *)       printf '' ;;
    esac
}
prov_recorded() {
    [ -f "$PROV" ] || return 1
    awk -F'\t' -v s="$1" -v t="$TAG" '$1==s && $2==t {print $3; found=1; exit}
                                      END{exit !found}' "$PROV"
}
prov_record() {
    mkdir -p "$(dirname "$PROV")"
    if [ -f "$PROV" ]; then
        awk -F'\t' -v s="$1" -v t="$TAG" '!($1==s && $2==t)' "$PROV" > "$PROV.tmp"
    else
        : > "$PROV.tmp"
    fi
    printf '%s\t%s\t%s\n' "$1" "$TAG" "$(sig_of "$1")" >> "$PROV.tmp"
    mv "$PROV.tmp" "$PROV"
}

# ONE skip decision, so the engine and the stages cannot disagree about it.
# They did: the engine printed "Rebuilding" and stage2_snmr.jl then refused,
# because that stage has its own input-signature guard and will not overwrite an
# output built from different inputs without --force. The run died at Stage 2
# AFTER Stage 1 had spent 3.5 h. Prints one of: absent, match, unrecorded, changed.
skip_state() {
    local out want have
    out=$(out_of "$1")
    if [ -z "$out" ] || [ ! -e "$out" ]; then echo absent; return; fi
    want=$(sig_of "$1")
    have=$(prov_recorded "$1") || have="__none__"
    if   [ "$have" = "$want" ];    then echo match
    elif [ "$have" = "__none__" ]; then echo unrecorded
    else                                echo changed
    fi
}

# Is the engine about to rebuild over an existing output? A stage with its own
# guard must be told, or it will refuse what the engine has already decided.
rebuilding_over_existing() {
    [ "$(skip_state "$1")" = absent ] && return 1
    [ "$FORCE" -eq 1 ] && return 0
    [ "$REUSE_STALE" -eq 1 ] && return 1
    case "$(skip_state "$1")" in match) return 1 ;; *) return 0 ;; esac
}

declare -a TIMED=()
FAILED=""

run_stage() {
    local name="$1"; shift
    local out; out=$(out_of "$name")
    local log="$LOGDIR/$name.log"

    if [ -n "$out" ] && [ -e "$out" ] && [ "$FORCE" -eq 0 ]; then
        local want have state
        want=$(sig_of "$name")
        have=$(prov_recorded "$name") || have="__none__"
        state=$(skip_state "$name")
        if [ "$state" = match ]; then
            echo "  SKIP  $name -- $out, built with these exact settings"
            TIMED+=("$name skipped")
            return 0
        fi
        # Existence is not provenance. Rebuild by default and say why: a wrong
        # answer in the right filename costs more than the time to redo it.
        if [ "$state" = unrecorded ]; then
            echo "  STALE $name -- $out exists but no run recorded building it,"
            echo "        so what produced it is unknown. Rebuilding."
        else
            echo "  STALE $name -- the settings changed since $out was built:"
            echo "          built with: $have"
            echo "          asked for : $want"
            echo "        Rebuilding."
        fi
        if [ "$REUSE_STALE" -eq 1 ]; then
            echo "        --reuse-stale given, so keeping it. Its provenance is NOT this run's."
            TIMED+=("$name reused-stale")
            return 0
        fi
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
        [ -n "$(sig_of "$name")" ] && prov_record "$name"
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
    case "$s" in truth|quicklook|stage1|draws|smrdata|stage2|site) ;;
        *) echo "unknown stage '$s' in STAGES"; exit 2 ;;
    esac
done

echo "$RUN_KIND run  tag=$TAG  label=${LABEL:-none}  tier=$TIER  conds=$CONDS"
echo "  scenario $TRUTH_CONFIG -> $TRUTH_NAME"
echo "  site     $SITE_CONFIG"
echo "  samples  $NSAMPLES   AEM chains $CHAINS   AEM draws $DRAWS"
echo "  logs     $LOGDIR"
# Stage 2's chain and ensemble counts are NOT settable: stage2_snmr.jl holds
# both at 4 as constants and reads no flag for either, so nothing here pretends
# to set them.
echo "  Stage 2 runs 4 chains x 4 ensembles regardless; that stage takes no flag for it"
[ "$RUN_KIND" = smoke ] && \
    echo "  SMOKE: proves the paths run. Do not quote a number from it."
echo

if wants truth; then
    run_stage truth bash -c "julia --project=. synthetics/scripts/build_truth.jl '$TRUTH_CONFIG' && \
                             cp '$TRUTH_SRC' '$TRUTH_DST'" || true
fi
# The cheap check. Deterministic Occam on both halves, so it needs no workers
# and no sample count: it reports whether each inversion converged and fit, and
# refuses to report anything that would look like an inference result.
wants quicklook && [ -z "$FAILED" ] && \
    run_stage quicklook julia --project=. tools/quicklook.jl "$TAG" || true

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
# stage2_snmr.jl keeps its OWN input-signature guard, which is right: it stops
# a hand-run stage from quietly overwriting an answer built from other inputs.
# But it means the engine has to tell it when the engine has already decided.
# Without this the engine said "Rebuilding", the stage refused, and the run
# ended at Stage 2 having spent 3.5 h on Stage 1. Worse, when the existing file
# carried no signature at all the stage SKIPPED and wrote nothing, the engine
# saw the old file still there and recorded THIS run's provenance against it --
# laundering a smoke output into one that looks like production.
if wants stage2 && [ -z "$FAILED" ]; then
    S2_FORCE=""
    rebuilding_over_existing stage2 && S2_FORCE="--force"
    run_stage stage2 julia --project=. stages/stage2_snmr.jl \
        --tier "$TIER" --conds "$CONDS" --tag "${TAG}T" --label "$LABEL" \
        --nsamples "$NSAMPLES" --ndraws "$DRAWS" $S2_FORCE || true
fi

# The trap that cost a reviewer an hour: a site config can name a different
# Stage 1 tag than the one just built (synth_2P.toml names 2M, synth_3P.toml
# names 3M, because the paper's Stage 1 ran under the older tag and sigma_bulk
# is identical). run_site.jl warns early but carries on, and dies in Stage 5 an
# hour later. Copy it here instead, and say so.
if wants site && [ -z "$FAILED" ]; then
    WANT_DRAWS=$(julia -e 'import TOML; print(get(TOML.parsefile(ARGS[1])["data"], "stage1_draws", ""))' "$SITE_CONFIG")
    SRC_DRAWS="data/aem_draws_${TAG}.jld2"
    if [ -n "$WANT_DRAWS" ] && [ -e "$SRC_DRAWS" ]; then
        if [ ! -e "$WANT_DRAWS" ]; then
            echo "  NOTE  $SITE_CONFIG wants $WANT_DRAWS; copying $SRC_DRAWS to it."
            echo "        Same sigma_bulk under a different tag. Stage 5 would die an hour in without this."
            [ "$DRY" -eq 1 ] || cp "$SRC_DRAWS" "$WANT_DRAWS"
        elif ! cmp -s "$SRC_DRAWS" "$WANT_DRAWS"; then
            # The copy is only made when the target is ABSENT, so on a second
            # run the target is last run's file and Stage 5 reads that instead
            # of this run's draws, silently. Not overwritten here because the
            # target may legitimately be the paper's own draws under its older
            # tag, which this run has no business replacing. So: say so.
            cat <<EOF
  WARNING  Stage 5 reads $WANT_DRAWS, which EXISTS and DIFFERS from this run's
           $SRC_DRAWS. It is NOT what this run produced -- most likely the
           previous run's draws, or the paper's own under its older tag.
           Either  rm $WANT_DRAWS   to have this run's copied in,
           or set [data].stage1_draws = "$SRC_DRAWS" in $SITE_CONFIG.
EOF
        fi
    fi
    run_stage site julia --project=. run_site.jl "$SITE_CONFIG" || true
fi

echo
# Guarded by length, NOT written as "${TIMED[@]}" bare: bash before 4.4 treats
# the expansion of an EMPTY array under `set -u` as an unbound variable and
# dies. macOS ships bash 3.2, so on a Mac an empty STAGES would have ended in
# an unbound-variable error instead of "complete".
if [ ${#TIMED[@]} -gt 0 ]; then
    printf '%s\n' "stage timings:" "${TIMED[@]/#/  }"
fi
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
        [ "$RUN_KIND" = smoke ] && \
            echo "  and it is a SMOKE table. Do not quote a number from it."
    else
        # "the loop ended" is not "it worked". Saying complete over a run that
        # produced no output is the defect this line exists to prevent.
        echo "INCOMPLETE -- no verdict table at $V; check $LOGDIR/site.log"
        exit 1
    fi
else
    echo "complete -- stages run: $STAGES"
fi
