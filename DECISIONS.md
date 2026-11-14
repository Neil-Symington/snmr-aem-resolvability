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

## Gate band definition

**The band is min chi2 subject to the forced window conductance and the stated
parameter bounds. Not chi2 at the MAP.** The gate asks whether any admissible
model fits the data; smoothness is a device for stabilising the solve, not a
criterion for admissibility. Min chi2 errs toward a wider band and a looser
bound, which is the conservative direction.

In practice: at each node the lowest chi2 is kept over all starts, both sweep
directions and lambda2 in {0.01, 0.1, 1, 10}. chi2 at the MAP for lambda2 = 1 is
recorded alongside so the difference is visible; the verdict uses min chi2. The
parameter bounds define admissibility and are written into the output, and a
node whose min chi2 is reached with parameters on their bounds is flagged.

The same standard applies to the unconstrained solve: multi-start, lowest chi2.
The pipeline then asserts `chiref <= min over nodes of the constrained chi2`
within tolerance and fails loudly otherwise. That is the consistency check that
would have caught the stall described in `gate_check/results/refit_diagnosis.md`.

**Parameter bounds are not a config value.** `src/occam_core.jl` holds the clamp
on log10 sigma in mutable `Ref`s, set by `setbounds!(lo, hi)` and read by
`getbounds()`; the default (-3, 1) corresponds to `sigma_bounds = [1e-3, 10]`.
The per-cell clamp is what limits how low the window conductance can be forced,
so extending the search below 1e-3 S/m (M1.1 item 2) is a `setbounds!` call and
must be recorded per run, not assumed from the config.

## Reference misfit (chiref)

**chiref = min(unconstrained multi-start, min over constrained nodes).** Both
numbers and their gap are reported per column. The old definition — the
unconstrained solve alone — is not reliable.

Measured on release 2P (truth `33bf02a46442711e`, window 17.5-33.5 m): the
unconstrained multi-start gave **28.4587** from three starts agreeing to three
decimals, while a constrained node reached **27.1749** at C = 2.8452 S. A
constraint cannot improve a fit, so the unconstrained solve is not the minimum,
and multi-starting it does not cure that. Same pathology as
`ek_realdata/REPORT_v1_firstfive.md:559-563`.

The consistency assertion stays as a **reported diagnostic with its gap, not a
hard failure**, and the lower of the two values is used — which is what the field
tables already quoted in practice.

This bears directly on the acceptance rule, which remains under review: **a
trustworthy chiref is what makes the relative rule usable at all.** Until chiref
can be relied on, `chi2 <= chiref + 3.84` is a threshold on a number that may be
too high, giving a band that is too wide — a looser bound, so conservative, but
wrong.

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
