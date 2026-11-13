# Decisions — synthetic study production re-run

14 Nov 2026. Branch `occam-workflow-test`. Each entry is a decision taken, not
a finding. Findings live in the reports named against them.

## Columns

**Test one = `synthetics/config/paper_3P.toml`.** Unchanged from the
manuscript. Described there as the cleaner, fresher case. Truth
`data/truth_murray3P.jld2`, hash `d9de86c430ab0165`. Confirmed the same
physical column as the manuscript's `HYP_SYNTH/truth_murray3P.jld2` (hash
`1cb14adbbb777aea`): all fields agree to 5.6e-17, the newer file carrying ten
extra keys. Package conductance 6–40 m = 0.9142 S in both.

**Test two = `synthetics/config/confined_aquifer.toml`**, the release column.
Truth `data/truth_murray2P.jld2`, hash `33bf02a46442711e`. **The older 2P truth
file is retired** (`HYP_SYNTH/truth_murray2P.jld2`, hash `6421d57fb6293968`).

**New numbers will not match the manuscript's, by construction.** The retired
column gives a package truth of 7.99 S over 4–34 m; the release column gives
5.60 S. The two are different columns, not two serialisations of one: the
retired one is 1.1–1.5x more conductive in nearly every unit at similar clay
fractions. A later reader must not read the difference as a regression.

**Saline versions of both are built from the same configs**, changing only the
target aquifer's pore-water conductivity to 0.30 S/m (about 1,920 mg/L at
`mg_per_sm = 6400`), from 0.08 on test two and 0.03 on test one. Configs
`confined_aquifer_saline.toml` and `paper_3P_saline.toml`.

## Sampler

Four independent ensembles per inversion. Layer limit 40. Seeds from config and
logged (M1.4). The chains are unseeded in the current code and this is what M1.4
changes.

**B2:** run by the reviewer on the release 2P only; layer limit 10; window
conductance medians still differ 1.28x; corrected ESS 22. Second column stopped.
See `b2_2P_diagnose.txt` and `b2_2P_multichain_ess.txt` in the reviewer's
bundle. **B3 not run.**

## Uncertainty

**Sampled bands on conductance are not used as uncertainty.** Coverage of truth
by the sampled band and by the gate range is reported for every column.

The decision path uses the gate's accepted range. On both columns checked the
true window conductance lies inside it: test one 0.2063 S in 0.0296–0.3447 S,
test two 0.6856 S in 0.0319–7.079 S. (Those gate ranges are themselves under
revision in M1.1; see `gate_check/results/edge_profile.md`.)

## Gate threshold — UNDER REVIEW

The accepted band is `chi2 <= N + 3.84` (`stage4_gate_v2.jl:260`,
`config.toml:385`), an **absolute** criterion. The 3.84 is the one-parameter 95%
chi-squared increment, specified in `PLAN_CHAINED_HYPOTHESES.md:237` as a
**relative** one, `delta-chi2 <= 3.84`. Anchoring it to N instead was a
deliberate guard against an unreliable chiref, the reason given in
`ek_realdata/REPORT_v1_firstfive.md:559-563`: the unconstrained solve was
observed to stall above the constrained sweep, so it is not reliably the
minimum.

Current status: **the rule is unchanged and still drives the verdict.** Both
bands are now computed and stored per column, with truth inside or outside for
each (`gate_check/results/both_rules.md`). Under the absolute rule truth is
inside on all four columns; under the relative rule it is outside on 3P and the
band is empty on 2P and 2S, the latter being direct evidence that the refits are
not genuine minimisations. Whether the relative rule becomes usable depends on
M1.1 item 3.

## Petrophysics

**The salinity step uses porosity from the prior**, U(0.25, 0.40), with the SNMR
water content as a **lower limit** on `theta_total` (M1.6, config switch, on by
default, reported both ways). This is a design decision, stated here and in the
paper. SNMR measures free water, which is a lower bound on total water; using it
as an equality would bias pore-water conductivity high by about 1.44x on these
columns.

**The SNMR kernel is built on the pooled AEM median** (condition C2). Neil's
standing choice. The per-draw kernel is parked.

## Status of earlier work

| what | state |
|---|---|
| B1, test two (release 2P) and test one (3P), four ensembles each | **superseded but retained** |
| Occam comparison, 2P and 2S | superseded but retained; out of scope from here |
| conditions table (`conditions_table/`) | superseded: clay terms were wrong by 5–16x |
| convergence diagnostics harness (`convergence_test/`) | retained, feeds M2 |
| gate edge profile (`gate_check/`) | current, feeds M1.1 |
