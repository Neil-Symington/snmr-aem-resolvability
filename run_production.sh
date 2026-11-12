#!/usr/bin/env bash
# run_production.sh -- the whole synthetic workflow at the committed parameters.
#
#   ./run_production.sh --dry-run    # print what it would run, run nothing
#   ./run_production.sh              # run it, about a day
#
# BUDGET ROUGHLY A DAY. Stage 1 is ~3.5 h at 200,001 samples over 4 chains,
# Stage 2 is hours, and Stage 4 inside the site stage is ~1.5 h. Start it under
# tmux or nohup so it survives closing the terminal:
#
#   nohup ./run_production.sh > production.out 2>&1 &
#   tail -f production.out
#
# IT WILL NOT REUSE SMOKE OUTPUTS. Three of the six outputs are named with
# neither the run kind nor the label -- aem_draws_<TAG>.jld2,
# smr_data_murray<TAG>T.jld2 and cond_*_<TIER><TAG>T.jld2 -- so a smoke file and
# a production file share a path and existence alone cannot tell them apart. The
# engine keeps a ledger of the settings each output was built under and skips a
# stage only on an exact match; anything else is reported STALE and rebuilt.
# That is why this is a separate script rather than a switch in a shared one.
#
# Run ./run_smoke.sh first on a column you have not run before. It costs an hour
# and tells you the paths work. Cheaper still, for "do the inversions even fit":
#
#   julia --project=. tools/quicklook.jl <TAG>

# ─── EDIT HERE ───────────────────────────────────────────────────────────────

TAG="2P"                                              # names every intermediate
TRUTH_CONFIG="synthetics/config/confined_aquifer.toml"
SITE_CONFIG="configs/synth_2P.toml"                   # FULL gate grid, 16+15
LABEL=""                                              # no suffix for production
TIER="typical"                                        # optimistic|typical|sceptical|field
CONDS="C2"                                            # the corrected log prior

# Drop a name to skip that stage. Order is fixed and enforced by the engine.
# `--from 3` on run_site.jl is the way to redo only the gate and downstream.
STAGES="truth stage1 draws smrdata stage2 site"

NSAMPLES=200001      # trans-D samples per chain
CHAINS=4             # AEM chains. Stage 2 runs 4 regardless; see below.
DRAWS=600            # AEM posterior draws carried to Stage 2

# ─── STOP EDITING ────────────────────────────────────────────────────────────

RUN_KIND="production"
. "$(dirname "$0")/run_stages.sh"
