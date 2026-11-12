#!/usr/bin/env bash
# run_smoke.sh -- the whole synthetic workflow at smoke settings.
#
#   ./run_smoke.sh --dry-run     # print what it would run, run nothing
#   ./run_smoke.sh               # run it, about 1 h
#   ./run_smoke.sh --force       # redo every stage, whatever is on disk
#
# WHAT SMOKE IS FOR. It exercises every stage and every file path and proves
# that the chain runs end to end. It is NOT a small version of the real answer:
# at 2001 samples the hierarchical noise scale runs far above its production
# value and the posterior is a different shape, and the gate grid is reduced
# from 16+15 nodes to 4+3. Numbers from it move by large factors between runs.
# DO NOT QUOTE ANY NUMBER FROM A SMOKE RUN.
#
# If all you want to know is whether the inversions work on a new column, the
# two-minute check is cheaper and says so honestly:
#
#   STAGES="truth quicklook"   below, or  julia --project=. tools/quicklook.jl <TAG>
#
# For the real thing use ./run_production.sh, which is the same workflow with
# the committed parameters. The two scripts exist separately so that no single
# switch has to be remembered: changing one line in a combined script used to
# leave the label and the site config behind and produce a smoke answer wearing
# a production label.

# ─── EDIT HERE ───────────────────────────────────────────────────────────────

TAG="2P"                                              # names every intermediate
TRUTH_CONFIG="synthetics/config/confined_aquifer.toml"
SITE_CONFIG="configs/smoke_2P.toml"                   # reduced gate grid, 4+3
LABEL="_smoke"                                        # keeps outputs separate
TIER="typical"                                        # optimistic|typical|sceptical|field
CONDS="C2"                                            # the corrected log prior

# Drop a name to skip that stage. Order is fixed and enforced by the engine.
STAGES="truth stage1 draws smrdata stage2 site"

NSAMPLES=2001        # trans-D samples per chain
CHAINS=2             # AEM chains. Stage 2 runs 4 regardless; see below.
DRAWS=50             # AEM posterior draws carried to Stage 2

# ─── STOP EDITING ────────────────────────────────────────────────────────────

RUN_KIND="smoke"
. "$(dirname "$0")/run_stages.sh"
