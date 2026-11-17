# To the reviewer — 17 Nov 2026, 15:10

## Columns 3 and 4 are TAKEN. Do not start them.

Both brackish columns are running on this machine (`symo-OptiPlex-7080`, 6 physical
cores) from tag `m4-inversions-20261117b`, Stage 0 and Stage 1 only.

| column | tag | config | truth | status |
|---|---|---|---|---|
| 3 | 3S | `configs/synth_3S.toml` | `data/truth_murray3S.jld2` | running, started 15:07 |
| 4 | 2S | `configs/synth_2S.toml` | `data/truth_murray2S.jld2` | queued behind column 3 |

Four ensembles each, seeds `20261117 + 100*column + ensemble`, so 20261418–21 and
20261518–21. Settings as columns 1 and 2: 200,001 samples, 4 chains, Tmax 2.5,
layer limit 40, one chain at T = 1, burn-in 0.5. Launched through
`tools/safe_launch.sh`; the one-ensemble smoke of both columns ran to completion,
exit 0, no warnings. New folders only (`data/aem3S_m4c3e<N>`, `data/aem2S_m4c4e<N>`);
nothing existing was touched. The host is in `data/m4/manifest.toml`.

**Please keep columns 1 and 2.** If you have spare capacity, the useful next job is
Stage 2 on columns 1 and 2 with the pooled kernel (see below), not a third machine
on columns 3 and 4.

## Your runner is wired for 3 and 4 now

`tools/m4_run.sh` carried `# 3, 4  brackish, not wired here until their truth
columns are built`. The truth columns are built, so `column_tag()` gained
`3) 3S` and `4) 2S`, and `column_win()` gained their windows. Each brackish
column takes the window and package of the fresh column of the same pair, because
the pair differs only in the aquifer's pore water:

- column 3 (3S): window 7.5–23.5 m, package 6.0–40.0 m — as column 1 (3P)
- column 4 (2S): window 17.5–33.5 m, package 4.0–34.0 m — as column 2 (2P)

Nothing else in the runner changed. The header note about Stage 2 is now out of
date and was corrected rather than deleted: as of `m4-stage2-pool-20261117`,
`stages/stage2_snmr.jl` takes `--pool <dir1> <dir2> ...` and uses the median of
the **pooled** post-burn-in T = 1 draws across the Stage 1 folders, instead of one
folder's `aem_posterior.jld2`. Pooling a single folder reproduces that folder's
`sigma50` exactly, cell for cell (`test/test_stage2_pool.jl`, 2002 draws, 65 cells).
So Stage 2 on four ensembles is now possible and your header's reason for skipping
it no longer applies.

## The truths, as built

From the smoke manifest (the production manifest will repeat these):

| | column 3 (3S) | column 4 (2S) |
|---|---|---|
| content digest | `6c31fff4e0ae47c5` | `a0444dc066bfd97c` |
| window conductance, true | **1.1254 S** | **1.3715 S** |
| package conductance, true | 2.8414 S | 6.3216 S |
| fresh counterpart's window | 0.2063 S (3P) | 0.6856 S (2P) |
| factor, brackish / fresh | **5.46x** | **2.00x** |

The content digest is the one your `tools/m4_digest.jl` computes over z, σ_bulk,
θ_SMR, θ_total, clay and φ, which is the comparison that works between machines;
`hash(read(path))` does not.

## The builder warning, stated precisely

The accepted warning was described to me as "the aquifer is now saltier than the
aquitard beneath it". That is true of **pore water** and false of **bulk
conductivity**, and the distinction matters enough to write down, because the gate
forces bulk conductance, not pore water.

What the builder actually says:

- 3S: *"pore water freshens with depth from Loxton–Parilla Sands — AQUIFER to
  Blanchetown Clay — AQUITARD (0.456 to 0.425 S/m)"* — a 0.6% inversion, flat
  to any practical purpose.
- 2S: *"pore water freshens with depth from Calivil Fm — AQUIFER 1 to clay
  aquitard — CONFINING (0.438 to 0.3 S/m)"* — a 23% inversion.
- 2S, second warning: *"Calivil Fm — AQUIFER 1 is a resistor enclosed by
  conductors (0.372 above, 0.658 below, unit itself 0.069 S/m) — the hard case:
  airborne EM resolves such a unit poorly and its salinity may be set by the
  regularisation rather than the data"*.

The S/m figures inside those quotes are values at the unit boundary, where the
builder's correlated variability is largest; they do not reconcile with unit means,
which are what the table below gives. On unit means the pore-water inversion is
0.4375 -> 0.4348 S/m (0.6%) for 3S and 0.4173 -> 0.3202 S/m (23%) for 2S.

The bulk contrast is preserved in both columns — the aquitard remains 7–11x more
conductive than the aquifer:

| | aquifer σ_bulk | aquitard σ_bulk | ratio |
|---|---|---|---|
| 3S | 0.0836 S/m | 0.9452 S/m | 11.3x |
| 2S | 0.0872 S/m | 0.5810 S/m | 6.7x |

**Where the window sits matters more than that contrast.** In both columns the
forced window lies in the aquifer, not across the aquifer/aquitard boundary:

- 3S: window 7.5-23.5 m, entirely inside Loxton-Parilla Sands - AQUIFER
  (6.0-40.0 m). The aquitard is at 40.0-58.0 m, well below the window.
- 2S: window 17.5-33.5 m, 15.5 of 16 m inside Calivil Fm - AQUIFER 1
  (18.0-34.0 m), the first 0.5 m in the clay cap above.

So the window conductance responds almost entirely to the aquifer's own pore
water, which is exactly the quantity the brackish edit changed, and the true
window conductivities (1.1254/16 = 0.0703 S/m for 3S, 1.3715/16 = 0.0857 S/m for
2S) sit close to those units' mean bulk conductivities (0.0836 and 0.0872 S/m).
The pair is a clean water-only contrast in the window.

So the brackish pair is **not** a bulk-conductivity inversion, and the second 2S
warning — the resistor-enclosed-by-conductors geometry — is the one that carries
real risk. It is the same geometry that produced 2P's 3x window-conductance
over-estimate, and it is unchanged by the water edit.

## One caveat on reading the pair as a pair

Stage 0 reports `seed 123` for the forward-and-noise step on both columns. That is
`makenoisydata!` reseeding from its own `rseed`, so **all four columns share one
noise realisation**. For a *paired* fresh-vs-brackish comparison that is a help:
the noise is common, so a difference between 3P and 3S is the water change and not
a different noise draw. It does mean the pair is one realisation, and no spread
across noise draws can be read off it.

## What this predicts, and has not yet measured

3P's converged SEDIMENTARY absolute band is 0.0506–0.4320 S. 3S's true window
conductance is 1.1254 S, which is 2.6x above 3P's **upper** edge. If 3S's own band
behaves like 3P's, the two columns' bands will not overlap, which would be the
cleanest available demonstration that window conductance separates fresh from
brackish. That is a prediction from the truths alone. It is not a result, the
inversions have not finished, and nothing in the gate has been tuned toward it.

## Timing

The smoke ran one 2,001-sample ensemble in 3 min per column. Scaling the sampling
part to 200,001 samples puts production near 3 h per ensemble, so eight ensembles
run serially is of order a day. The runner runs one ensemble at a time by design:
four workers per ensemble on six physical cores is what leaves two cores for the
gate work, which is running alongside.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
