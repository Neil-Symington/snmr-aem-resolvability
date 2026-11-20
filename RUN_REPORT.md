# RUN_REPORT — synthetic study, production re-run

Started 18 Nov 2026. This file records what was run, what it produced, and what
was superseded or lost. It is the provenance record for the paper's numbers.

Branch `m1-production-rerun`. Tags: `m4-inversions-20261117b` (Stage 0 and 1
inputs and runner), `m4-stage2-pool-20261117` (pooled Stage 2 kernel),
`snmr-data-truth-20261118` (SNMR generating kernel on the true conductivity).

## Superseded and lost results

### The 13 Nov SNMR results are superseded

`tools/smr_data.jl` built the **generating** kernel on an AEM inversion's
conductivity rather than on the truth, so the synthetic SNMR measurements every
Stage 2 run before 18 Nov inverted had been forward-modelled through an inversion
result. Measured against the true `σ_bulk` at the same depths, the conductivity
used departed by a median 55% (3P), 54% (2P), 45% (3S) and 52% (2S), and was
bit-identical to `aem_occam<TAG>.jld2`.

Corrected data differ from the old data by a median 0.5–5% in amplitude but by
**up to 24σ (3P), 52σ (2P), 71σ (3S) and 46σ (2S)** in units of the measurement
noise, which is what an inversion sees. The 13 Nov results are therefore
superseded on the numbers, not only on provenance.

Superseded files were **moved, not deleted**, to `data/superseded_20261118/`:
six `smr_data_murray*.jld2` and five Stage 2 outputs built on them
(`cond_C2_occam{2P,2S,3P,3S}.jld2`, `cond_C2_typicalFIXT.jld2`). Detail in
`SNMR_DATA_FAULT_20261118.md`.

### Run A of Stage 2 on 3S and 2S was LOST, 18 Nov

The first full-length pooled Stage 2 runs on 3S and 2S (both completed 07:38,
18 Nov) **no longer exist**. They were overwritten at 10:51:58 and 10:52:04 by
re-runs, because Stage 2's output name was a function of the column alone and
`tools/stage2_run.sh` passed `--force` to get past the provenance guard that had
been refusing exactly this. There is no archived copy anywhere in the tree.

The loss matters because the runs are not bit-reproducible — the seeding is inert
(see `DECISIONS.md`, 18 Nov) — so a re-run is a different realisation, not a
reproduction. The only surviving record of run A is this table:

| 3S, pooled kernel | run A, 07:38 (LOST) | run B, 10:52 (on disk) |
|---|---|---|
| window water p50 | 0.2710 v/v | 0.2672 v/v |
| range of the four ensemble medians | 0.2653–0.2741 | 0.2611–0.2692 |
| truth 0.2704 against that range | inside | above all four |
| pooled p05–p95 | 0.2631–0.2763 | 0.2602–0.2733 |
| truth inside the pooled interval | yes | yes |
| R-hat / ESS bulk, window average | 1.507 / 9 | 1.720 / 7 |
| s p50 per ensemble | 0.84 / 1.77 / 1.99 / 1.82 | 1.71 / 2.35 / 1.54 / 1.91 |
| chi2 per datum vs true noise | 2.96 / 1.54 / 1.81 / 1.95 | 1.92 / 3.01 / 1.71 / 2.26 |

2S's run A was lost the same way. Its reported figures were: s p50 per ensemble
1.15 / 1.42 / 1.27 / 1.29, chi2 per datum 2.12 / 0.99 / 0.97 / 2.23, window water
pooled p50 0.2474 (p05–p95 0.2358–0.2618), R-hat 1.545, ESS bulk 8.

Taken together, run A and run B are eight independent ensembles whose medians
span 0.2611–0.2741 and bracket the truth at 0.2704 — which is why Stage 2 went to
**eight ensembles per column**, and why **every run now writes its own folder**
and `--force` is never passed over a full-length result.

### A day's worth of Stage 2 was unrunnable, 17–18 Nov

`stages/stage2_snmr.jl` could not load at all between commits `a39a3e9` and the
18 Nov repair: a docstring added for `pooled_sigma` was inserted between
`inputsig`'s docstring and `inputsig`, leaving one docstring documenting another.
It parses cleanly and fails only on evaluation. Closed by
`test/test_stage_loads.jl`.

### The reviewer's Stage 1 folders arrived without a manifest

3P and 2P's four Stage 1 ensembles each arrived on 18 Nov as two zip archives,
`aem3P_m4c1e.zip` and `aem2P_m4c2e.zip`, 57 MB apiece. All three arrival checks
pass: eight of eight directories with all five expected files non-empty at 14 MB
each; every seed exactly `20261117 + 100*column + ensemble`, so 20261218-21 and
20261318-21, at 200,001 samples and 4 chains; and the truth window conductances
0.2063 S (3P) and 0.6856 S (2P), which are the values the gate's acceptance
targets reference.

**No manifest came with either archive.** The reviewer's machine writes
`data/m4/manifest.toml` recording its host, commit, Julia version, runner
checksum and the truth content digest it computed, and that file was not
included. The consequence is specific: the content digests **cannot be compared
across machines**. The checks above were run against the truth files on this
machine, so they confirm that these Stage 1 outputs are consistent with the
local truths and carry the specified seeds; they do not confirm that the
reviewer's machine built them from a bit-identical truth column. For the
record, the local digests are `853f0cd9c471432b...` (3P) and
`4c5811b2e3dc1643...` (2P), over z, sigma_bulk, theta_SMR, theta_total, clay
and phi.

Nothing in the run depends on closing that gap, but it is the one provenance
link in the chain that is asserted rather than verified.


### The lower band edge was the upper edge, 19 Nov

`gate_check/edge_bisect.jl` reported the lower band edge wrongly from the day it
was written until commit 54512a1 (tag `edge-lower-reflect-20261119`).
`EdgeFinder.bracket` locates an UPPER crossing only — highest accepted node, and
the lowest rejected node above it. The tool computed the correct lower bracket,
printed it, and then called the module on the UNREFLECTED solved set, so the
"lower edge" search re-refined the upper edge and returned it.

On 3S the two sides came out 1.7957 S and 1.7963 S (absolute rule), and 1.6641 S
and 1.6654 S (relative rule): the same edge twice, against a printed lower
bracket of −1.50 dex (chi2 28.687) to −1.75 dex (chi2 48.209).

**Blast radius: 3S only.** Every other tag short-circuits to `state = "open"`
before the call, having no rejected node below the accepted range — 3P,
3P_sed2, 2S and the three smokes. 2P never went through this tool.
`stages/stage4_gate_v2.jl:433-441` reflects inline and is CORRECT, so the
release path was never affected; it is left untouched with the freeze pending.

The reflection now lives in `EdgeFinder` behind `side = :lower`.
`test/test_edgefinder.jl` case 5 puts a V-shaped chi2 with both crossings to it:
against the old code the lower search returns −1.45629, the upper crossing;
the fix returns −1.74371 against a truth of −1.74436.

### 3S upper edge, absolute rule: NOT CONVERGED, 19 Nov

Reported from the solved nodes, not re-run (Neil's ruling, 19 Nov). Window
thickness 16.0000 m, threshold chi2 ≤ 45.840.

| | dex | chi2 | C |
|---|---|---|---|
| highest accepted | −0.950195 | 45.414 | **1.7944 S** |
| lowest rejected | −0.949707 | 45.860 | **1.7964 S** |

The bracket is 0.0005 dex wide, so the edge is pinned to about 0.1%. The
NOT CONVERGED verdict is from the third convergence criterion — a bracket node's
chi2 must not move by more than 0.5 when re-solved warm from the newest midpoint
— and not from bracket width.

The relative-rule upper edge CONVERGED at **1.6641 S** (state `located`).

Five warnings were raised during that search, where re-solving a node warm from
a nearer start RAISED its chi2, which continuation should make impossible:

```
node -1.0:     27.645   -> 28.574
node -0.75:  1427.356   -> 2122.094
node -1.0:     27.645   -> 40.697
node -1.0:     27.645   -> 35.472
node -0.9375:  62.516   -> 64.603
```

Three more came from the lower-edge run of 20 Nov, on the same column:

```
absolute rule
  node -1.75:    48.209   -> 48.993
relative rule
  node -1.625:   37.786   -> 38.559
  node -1.5:     27.922   -> 28.440
```

Eight in total on 3S. Node −1.0 is the accepted bracket node of the upper edge
and moved three times, once by 13.05.
The stored minimum is never overwritten by a worse solve, so the reported edge
is unaffected; what the warnings show is that the Occam solve at a forced
conductance is not a deterministic function of its start. Not chased.

## What is current

| stage | 3P | 2P | 3S | 2S |
|---|---|---|---|---|
| Stage 0 + Stage 1, four ensembles | on the reviewer's machine | on the reviewer's machine | **done 18 Nov 05:42** | **done 18 Nov 05:42** |
| SNMR data, true-conductivity kernel | done | done | done | done |
| Stage 1, four ensembles, from the reviewer | **arrived and verified 18 Nov** | **arrived and verified 18 Nov** | n/a | n/a |
| Stage 2, eight ensembles, pooled kernel | **done, 21 min** | **done, 21 min** | **done, 22 min** | **done, 22 min** |
| Stage 2 control, true-conductivity kernel | queued after the freeze | queued after the freeze | **done** | **done** |
| Stage 0 and Stage 3 through the release path | **done** | **done** | not yet | not yet |
| Stage 4 gate, acceptance windows | running | running | n/a | n/a |

Stage 1 for 3S and 2S: 200,001 samples, 4 chains, Tmax 2.5, layer limit 40, one
chain at T = 1, burn-in 0.5; seeds `20261117 + 100*column + ensemble`; eight of
eight directories complete and verified; manifest `data/m4/manifest.toml`.
Concurrency two ensembles at a time, measured at 1.16x the single-run time per
ensemble against a 1.7x margin.

All columns share one AEM noise realisation and one SNMR noise realisation
(`seed = noiseseed = 20261012`), which is intended: it makes the fresh/brackish
comparison paired. No spread over noise draws can be read from these runs.
