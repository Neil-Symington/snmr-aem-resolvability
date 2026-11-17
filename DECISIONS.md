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

**Brackish versions of both are built from the same configs**, changing only
the target aquifer's pore-water conductivity to 0.4375 S/m (about 2,800 mg/L at
`mg_per_sm = 6400`), from 0.08 on test two and 0.03 on test one. Configs
`confined_aquifer_brackish.toml` and `paper_3P_brackish.toml`. These columns
are called "brackish", not "saline", in reports, tables and figures.

**Second pair of columns set to about 2,800 mg/L, the upper end of the brackish
class, by Neil, 17 Nov. Previously about 1,920 mg/L.** The configs were
`confined_aquifer_saline.toml` and `paper_3P_saline.toml` at 0.30 S/m until
then; they were renamed because no code reads them by name. Stage code is
identical to tag `m4-inversions-20261117`.

**Truth files are identified by content digest from 17 Nov**, not by
`hash(read(path))`: two machines building the same column write different
bytes. `tools/m4_digest.jl` gives the SHA-256 of the raw Float64 values of z,
σ_bulk, θ_SMR, θ_total, clay and φ, and the window and package conductances.

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

## 3P gate edges, converged (M1.1)

Located by warm-started bisection, `gate_check/edge_bisect.jl`, both edges
meeting the convergence test on two successive refinements (edge moving less
than 0.02 dex and chi2 changing less than 0.5). 48 solves.

| rule | threshold | edge | state |
|---|---|---|---|
| absolute, N + 3.84 | 45.840 | **0.4594 S** | converged (final move 0.0008 dex) |
| relative, chiref + 3.84 | 28.929 | **0.3667 S** | converged (final move 0.0002 dex) |
| old single-pass gate | 45.840 | 0.3447 S | reported as OPEN, was solver reach |

Both lower edges are OPEN: the bound limit at sigma 1e-4, not closed by data.

Clay-free salinity ceiling at those edges, mg/L TDS. **The SEDIMENTARY figures
below are at the accepted bracket node, not a converged edge** -- bisection on
that variant is running; the WIDE figures are converged.

| edge | widest | central | tightest |
|---|---|---|---|
| absolute 0.4594 S, prior theta >= 0.25 | **2560** | 1242 | 726 |
| absolute 0.4594 S, SNMR theta >= 0.265 | 2291 | 1195 | 726 |
| relative 0.3667 S, prior theta >= 0.25 | 2043 | 991 | 580 |
| relative 0.3667 S, SNMR theta >= 0.265 | 1829 | 954 | 580 |

**The widest-corner ceiling is below 3,000 mg/L under both rules**, so test one
rules out the saline class at the most permissive corner of the porosity and
exponent box. That survives the edge moving out from the old 0.3447 S to
0.4594 S, which raised the widest-corner ceiling from 1921 to 2560 mg/L.

Caveat from the whole-profile report: models at the top of the band carry ten
parameters pinned at their clamps, so the upper edge is partly a statement about
the bounds. See `gate_check/results/whole_profile_3P.md`.

## Fresh and brackish pairs

**The fresh and brackish columns of each pair differ only in the aquifer's pore
water.** Nothing else in either config changes: same geometry, same porosity,
same clay fractions, same clay conductivity, same loop and same noise.

**Builder warning, noted and accepted (Neil, 17 Nov):** at 0.4375 S/m the
aquifer is now more conductive than the aquitard beneath it. That is left as it
is. It is a deliberately hard case for the gate -- the conductive target sits
above a less conductive confining layer, which is the opposite of the usual
contrast -- and changing it would weaken the control.

## 3P gate, converged under both variants (17 Nov)

Config `synthetics/config/paper_3P.toml`; window 7.5-23.5 m.

| variant | rule | edge | state |
|---|---|---|---|
| WIDE | absolute | 0.4594 S | converged |
| WIDE | relative | 0.3667 S | converged |
| **SEDIMENTARY** | **absolute** | **0.4320 S** | **converged**, final move 0.0002 dex |
| SEDIMENTARY | relative | 0.3605 S | converged, final move 0.0009 dex |

Containment holds: 0.4320 <= 0.4594 and 0.3605 <= 0.3667. **Truth at 0.2063 S is
inside the SEDIMENTARY band (0.0506-0.4320 S).** The lower edge is set by the
floor, minimum reachable 0.016 S.

Clay-free ceiling at the SEDIMENTARY absolute edge: **2,407 mg/L** at the widest
porosity-exponent corner, 1,168 central, 683 tightest; **2,155 mg/L** at the
widest corner with the SNMR porosity limit. Below 3,000 mg/L, so test one rules
out the saline class at the most permissive corner under the variant that drives
the verdict.

**The stated limits moved the ceiling by only 6%** -- 2,560 to 2,407 mg/L -- while
narrowing the band substantially at the bottom. The porosity and exponent ranges
dominate the ceiling: a factor of 3.5 across their box against 6% from the
geological bounds. **To be said plainly in the M5 report.**

## Which variant and which rule drive the verdict (17 Nov)

**SEDIMENTARY (0.001-2 S/m) drives the verdict. WIDE (1e-4 to 10 S/m) is
reported for comparison.** Both bands, ceilings and pinned counts appear side by
side; the data-only result is never dropped.

**The ABSOLUTE rule (chi2 <= N + 3.84) drives the verdict. The relative rule
(chi2 <= chiref + 3.84) is reported alongside.**

Reason, measured rather than assumed: on 2P the relative threshold rests on a
single bound-pinned node. The constrained minimum is 25.3414 **at C = 0.0051 S,
the lowest node on the grid, whose model has 8 cells on bounds**, and it sits
3.1172 below the unconstrained solve's 28.4587. So the relative threshold on
that column is anchored to the least trustworthy model in the sweep. On 3P both
rules give the same class result, so nothing is lost by preferring the absolute
one.

**The lower edge under SEDIMENTARY is set by the floor, not by the data.** A
floor of 0.001 S/m over a 16 m window cannot produce a window conductance below
**0.016 S**, so nodes below that return no admissible solve at all. It is
reported as "set by the floor, minimum reachable conductance 0.016 S", never as
a data limit and never as open.

### 2P, in plain words

**chi2 varies by about 0.5 across window conductances from 0.005 to 3 S.**
Measured: 25.34, 25.42, 25.53, 25.83, 25.65, 25.44, 25.53 at C = 0.0051, 0.0160,
0.0506, 0.1600, 0.5060, 1.6000 and 3.1924 S. The band on test two is not wide
because the edges are uncertain; it is wide because the likelihood is flat over
nearly three orders of magnitude. That is the finding, and it is a statement
about what the AEM data can resolve in a conductive column, not about the gate.

## Gate bounds: two variants (approved 16 Nov)

**Broad physical limits for a sedimentary system; no lithology assumed; SNMR
supplies a floor inside the saturated interval.**

| variant | log10 sigma bounds | sigma | role |
|---|---|---|---|
| WIDE | (-4, 1) | 1e-4 to 10 S/m | kept for comparison |
| **SEDIMENTARY** | (-3, 0.30103) | **0.001 to 2 S/m** | **drives the verdict** |

Applied to every cell alike. In the SEDIMENTARY variant only, the floor inside
the SNMR-defined window is the larger of 0.001 S/m and

    sigma_w(freshest plausible water) * theta_low ^ m_max

with theta_low the lowest p05 of window water content across the ensembles.

**At the provisional 50 mg/L the SNMR floor does not bind on either column.**
sigma_w = 0.00781 S/m gives a floor of 0.000627 on 3P (theta_low 0.265) and
0.000440 on 2P (theta_low 0.2199), both below the sedimentary 0.001. The
mechanism is built and inert at this value. It would begin to bind above
**80 mg/L on 3P and 114 mg/L on 2P**.

**Unit-mean check, all four columns: PASS.** All 26 unit means lie inside
0.001-2 S/m, the extremes being 2P/2S unweathered Palaeozoic at 0.00656 S/m and
3P/3S Blanchetown Clay at 0.944 S/m.

**Per-unit ranges**: config hook built, left off. For users with a driller's
log; not used in the paper's runs.

### Cells below the sensitivity depth: EXCLUDED FROM COUNTS, not fixed

Asked to choose between fixing them to the reference model and excluding them
from pinned counts, I chose **excluding them**, for four reasons.

1. Fixing cells removes degrees of freedom the data do not constrain, which
   lowers chi2 artificially and can only **narrow** the band -- the direction
   that risks excluding the truth.
2. Below the sensitivity depth the reference model is itself arbitrary, so
   fixing to it asserts a value the data cannot support, and it still enters the
   forward response however weakly.
3. Excluding is purely a **reporting** change. It achieves what was asked --
   not penalising a node for pinning where there is no sensitivity -- without
   touching the solve.
4. It is auditable: counts above and below the sensitivity depth are reported
   separately, so a reader can recombine them.

The sensitivity depth is computed from the noise-weighted Jacobian, as the
column norm of W*J per cell normalised by its maximum, with the depth taken
where that falls below a stated fraction.

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

## The gate's acceptance evidence: true-model check, 17 Nov

`gate_check/truecheck.jl`, both release columns, SEDIMENTARY bounds, absolute
rule (threshold 45.84 on 42 data). Full write-up in
`gate_check/results/TRUECHECK.md`.

**Min chi2 at the true window conductance, seeded from the generating model, is
26.205 on 3P (C_true = 0.2063 S) and 27.964 on 2P (C_true = 0.6856 S).** Both are
accepted. The cost of forcing the truth's conductance over the best available
unconstrained fit is **+0.92 (3P)** and **+2.11 (2P)**, inside the 3.84 the
relative rule allows — so the truth is accepted under the **relative** rule as
well as the absolute one, on both columns. Both solves met the constraint to
~1e-16 with **no cell on a bound**, so neither is a clamp artefact.

**This is the gate's acceptance evidence.** It answers "would this gate accept the
truth?" by measurement rather than assumption: on both release columns, yes, under
either rule, near the unconstrained optimum, with nothing pinned.

**The generating model's own chi2 is 35.449 on BOTH columns, and that is a
property of the noise realisation, not of either column.** The noise enters as
sigma times a standard normal draw and `makenoisydata!` reseeds from its own
`rseed`, so the generating model's chi2 is the sum of 42 squared draws from the
one shared realisation. It is not a test of a column and must not be reported as
one.

**"chi2 of the discretised truth" is not a well-defined quantity.** The same 3P
truth gives chi2 35.4 under the `piecewise_interp` mapping Stage 0 generates with,
146.1 under a thickness-weighted arithmetic average and 104.6 under a geometric
one, while the three agree on the window conductance to 0.1%. The spread is thirty
times the 3.84 allowance. Any statement about the truth's misfit must name the
mapping. This is a reason for defining the gate on forced conductance with
minimisation over the profile, which no choice of mapping can move. **One sentence
to this effect goes in the report** (M5).

**On 2P the grid cannot represent C_true better than 8.4% high** (0.7435 S for a
true 0.6856 S), because the window 17.5-33.5 m straddles the 18.0 m clay-cap
boundary. On 3P, whose window sits inside one unit, the same mapping is 0.3% low.
This is a floor on recovery, in the direction of over-stating conductance; it is
small against 2P's ~3x over-estimate and does not explain it. **The 8.4% is stated
beside the 2P results**, not in a methods aside.

**Trap recorded:** `occam`'s `fixval` under `ctype = :cond` is conductance in S
(`src/occam_core.jl:172`), not conductivity. Passing `C/thick` forced one
sixteenth of the intended conductance; on 3P that was below the floor-limited
minimum of 0.0160 S, so every window cell pinned to the floor and the solver
returned `settled = true` with a plausible chi2 at the wrong conductance.
`truecheck.jl` now verifies the achieved conductance from the returned model.
`profile_minchi2.jl` never had the bug: it requires `conres < 1e-8` before a
solve enters the band, which is also why `Inf` at a profile node means **no solve
met the constraint**, not a rejection.

## 2P coarse upper edge, 17 Nov

`gate_check/results/EDGE_2P.md`. SEDIMENTARY, absolute rule. **The upper edge lies
between 6.3697 S (accepted, min chi2 36.577 pooled) and 7.1469 S (rejected, 53.923)
— bracketed to 0.0500 dex**, better than the 0.1 dex asked. With the 9-node grid,
2P's absolute band is **0.0506–~6.7 S, a factor of 126**: its window conductance is
effectively unconstrained by the data over two orders of magnitude, which is the
same fact as the known flatness of its likelihood.

**Relative rule, from the existing coarse nodes only:** 0.0506–5.0596 S, upper edge
bracketed by 5.0596 accepted / 6.3697 rejected at `chiref + 3.84 = 29.70`, a width
of 0.0999 dex.

The node at 8.0190 S remains `Inf` after cold starts from `mref` and `m0` and
lambda2 in {0.01, 0.1, 1, 10}, so the earlier `Inf` there was not a warm-start
artefact. It does not matter: 7.1469 S below it solves cleanly and is rejected, so
the edge is bracketed without it and no further solves were spent.

**Two lines of that run's output must not be quoted:** `truth … absolute OUTSIDE`
and `relative band NaN – NaN` are artefacts of a three-node grid beginning at
6.3697 S, above the truth at 0.6856 S. The truth's standing is the acceptance
evidence above.

**A single profile run gives an upper bound on min chi2, not min chi2.** At
6.3697 S this run reached 37.888 where the 9-node run reached 36.577, with a
*smaller* lambda2 set — the difference is which neighbour the warm chain arrived
from. Quote the lowest solve ever obtained at a node. No verdict changes.

## Pinned cells by depth, both variants, 17 Nov

`gate_check/results/PINNED_TWO_VARIANT.md`. From stored models only. Absolute
rule. Deep counts reported separately, never folded in.

| column | variant | absolute band | largest clamp-free C | truth |
|---|---|---|---|---|
| 3P | WIDE | …–0.4509 S | 0.3192 S | 0.2063 S, below |
| 3P | SEDIMENTARY | 0.0506–0.4320 S | **0.1795 S** | 0.2063 S, **above** |
| 2P | WIDE | 0.0051–6.3697 S | 6.3697 S, all of it | 0.6856 S, below |
| 2P | SEDIMENTARY | 0.0506–~6.7 S | 5.0596 S | 0.6856 S, below |

"Clamp-free C" is the largest accepted forced conductance with no **shallow** cell
on a bound. **One case of four puts the truth above its clamp-free limit: 3P under
SEDIMENTARY, 0.2063 S against 0.1795 S.** Near the truth the pinning is one cell
(47–50 m) at 0.2014 S and three (27–31 m and 47–50 m) at 0.2260 S. The truth is
well inside the band itself and inside the relative rule; only the
no-cell-on-a-bound subset stops short. **That is direct evidence for not adopting
"no cell on a bound" as a criterion, and against the release path that still
does.**

**Where the pinning happens.** On 3P it is never inside the forced window. It sits
at 0–9 m above the window and 22–36 m plus 47–53 m below it — the compensation the
forced conductance demands, reaching the floor. At 3P's converged edge of 0.4320 S,
13 of 65 cells are on a bound and all 13 are shallow, which is why that edge is
reported as an edge under stated limits, not a data-only bound.

**Defect found in the sensitivity depth, not fixed.** `sensitivity_depth.jl` takes
the deepest cell reaching `frac` of the maximum sensitivity. Both columns carry a
secondary bump near 180–240 m (0.0093 on 3P, 0.0119 on 2P). At the default
`frac = 0.01` the bump clears the threshold on 2P and not on 3P, so 2P's split
depth reads 207.9 m with 9 cells excluded against 3P's 84.6 m with 31 — the wrong
way round, since 2P's conductive clay cap must reduce the depth of investigation,
not quadruple it. A 200 m feature at 1% of the peak is numerical.

It moves 2P's counts a long way: at C = 0.0506 S the same 17 cells on a bound
split 10 shallow / 7 deep at frac 0.01 and 1 shallow / 16 deep at frac 0.02. On 3P
the two fracs give identical tables. **The report therefore quotes frac = 0.02 and
prints frac = 0.01 in full alongside**; the definition was not changed, `frac`
being a documented parameter. A robust definition (an envelope, or a crossing on
the decaying limb) is wanted before the sensitivity depth carries anything
load-bearing. It does not affect any verdict here: only 2P's WIDE clamp-free limit
moves (6.3697 vs 5.0596 S), and its truth at 0.6856 S is far inside either.

## BLOCKER ON THE FREEZE: the release gate still uses "no cell on a bound"

Found 17 Nov while adding the conductance assertion. **Not changed — this is a
decision, not a fix, and the band it sets is the paper's band.**

`stages/stage4_gate_v2.jl:198` sets

    good[i] = nf.nbound == 0 && nf.conres < 1e-8 && nf.settled

and `good` is not merely reported. At line 270 the band is computed from
`chis[.!good] .= NaN`, and in the two-pass mode the coarse sweep's
`reach = findall(SC.good)` sets the **bracket of the fine grid**, so a node with
any cell on a bound is labelled "unreachable" and the fine search never visits it.

That contradicts the band definition recorded above — min chi2 subject to the
forced window conductance **and the stated parameter bounds** — under which a
solution with cells on a bound is a legitimate constrained minimum. It is also the
criterion that was explicitly set aside ("do not adopt 'no cell on a bound' as a
criterion"), and M1.2's whole point was to make per-cell stated limits part of the
definition rather than a disqualification.

Measured on the existing SEDIMENTARY coarse grids, absolute rule, threshold 45.84:

| column | band as defined (bounds legitimate) | band with `nbound == 0` | nodes lost |
|---|---|---|---|
| 3P | 0.0506–0.3582 S (6 of 9 nodes) | **0.1426–0.1795 S** (2 nodes) | 4 |
| 2P | 0.0506–6.3697 S (7 of 9 nodes) | 0.1600–5.0596 S (4 nodes) | 3 |

The nodes lost on 3P are 0.0506 S (chi2 25.305, nbound 3), 0.2260 S (25.471,
nbound 1), 0.2845 S (26.500, nbound 4) and 0.3582 S (28.930, nbound 15) — every
one of them far below the 45.84 threshold. The criterion is not excluding bad fits;
it is excluding good fits that happen to touch a bound.

**The consequence that blocks the freeze: the truth falls outside 3P's band under
the release path's criterion.** C_true = 0.2063 S against a `nbound == 0` band of
0.1426–0.1795 S. Under the definition as recorded, the band is 0.0506–0.4320 S
(converged edge) and the truth is comfortably inside, which is the acceptance
evidence above. So the two differ on the one question the gate exists to answer,
and 3P's band narrows from a factor of 7.1 to a factor of 1.26.

The minimal change consistent with everything already recorded is

    good[i] = nf.conres < 1e-8 && nf.settled

with `nbound` reported per node, as the profile driver already does. That is a
change to the release band and needs Neil's decision, so it is **not** applied.
Freezing before it is settled would freeze a gate that rejects the truth on 3P.

## All columns share one noise realisation

Stage 0 reports `seed 123` on every column; `makenoisydata!` reseeds from its own
`rseed`, so **3P, 2P, 3S and 2S all carry the same noise draw**. Two consequences,
both standing:

- The generating model's chi2 is the same constant (35.449) on every column,
  being the sum of 42 squared draws from that one realisation. **It is not to be
  reported per column.**
- A paired fresh-vs-brackish comparison is helped by this: the noise is common, so
  a 3P-to-3S difference is the water change and not a different draw. It also
  means every result rests on one realisation, and no spread over noise draws can
  be read from any of it.

## Columns 3 and 4, brackish pair

Built and running on this machine from tag `m4-inversions-20261117b`, Stage 0 and
Stage 1 only, four ensembles each, seeds `20261117 + 100*column + ensemble`.
`tools/m4_run.sh` wired for `3 -> 3S` and `4 -> 2S`, each taking the window and
package of the fresh column of its pair. True window conductances **1.1254 S
(3S)** and **1.3715 S (2S)**, against 0.2063 S (3P) and 0.6856 S (2P) — factors
of 5.46 and 2.00. Content digests `6c31fff4e0ae47c5` and `a0444dc066bfd97c`.

In both columns the forced window lies in the aquifer, not across the
aquifer/aquitard boundary (3S: 7.5-23.5 m inside Loxton-Parilla 6.0-40.0 m;
2S: 15.5 of 16 m inside Calivil 18.0-34.0 m), so the window conductance responds
to the aquifer's own pore water — the quantity the brackish edit changed.

**The second 2S warning is the case the paper is about.** Calivil Fm is a resistor
enclosed by conductors (clay cap above, clay aquitard below), and airborne EM
resolves such a unit poorly, so its salinity can be set by the regularisation
rather than the data. The water edit does not touch that geometry. This is the
column's purpose, not a defect in it.

**The accepted builder warning concerns pore water, not bulk conductivity.** On
unit means the pore-water inversion with depth is 0.4375 -> 0.4348 S/m (0.6%) on
3S and 0.4173 -> 0.3202 S/m (23%) on 2S, while the bulk contrast is preserved:
the aquitard stays 11.3x (3S) and 6.7x (2S) more conductive than the aquifer. The
2S builder issues a **second** warning that was not part of the accepted one --
Calivil is a resistor enclosed by conductors, "the hard case" -- and that geometry
is unchanged by the water edit.

## Status of earlier work

| what | state |
|---|---|
| B1, test two (release 2P) and test one (3P), four ensembles each | **superseded but retained** |
| Occam comparison, 2P and 2S | superseded but retained; out of scope from here |
| conditions table (`conditions_table/`) | superseded: clay terms were wrong by 5–16x |
| convergence diagnostics harness (`convergence_test/`) | retained, feeds M2 |
| gate edge profile (`gate_check/`) | current, feeds M1.1 |
