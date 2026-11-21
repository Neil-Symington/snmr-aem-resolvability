# Test two v2 (2Pv2 / 2Sv2) -- unit table

Built 21 Nov 2026 from `synthetics/config/confined_aquifer_v2.toml` and
`confined_aquifer_brackish_v2.toml`. Same units, depths, water and window
(17.5-33.5 m) as v1; the only change is `[variability]`, which is now test
one`s verbatim: jitter logit, corr_length 1.0, log_sd_phi 0.35,
log_sd_clay 0.45, log_sd_sigma_w 0.25, seed 20261012.

The two columns differ ONLY in the aquifer`s pore water. Verified cell by
cell: `sigma_w` and `sigma_bulk` differ in 84 of 800 cells, 15.875-36.625 m
(the aquifer plus its gradational contacts); porosity, clay, saturation and
both water contents are bit-identical.

| top (m) | base (m) | unit (internal) | general name | lithology role | porosity | clay | pore water 2Pv2 | pore water 2Sv2 |
|---|---|---|---|---|---|---|---|---|
| 0.0 | 4.0 | Shepparton Fm, unsaturated | unsaturated clay | mixed | 0.33 | 0.40 | 0.05 S/m (320 mg/L) | 0.05 S/m (320 mg/L) |
| 4.0 | 18.0 | Shepparton Fm, clay cap | clay aquitard | mixed | 0.35 | 0.45 | 0.05 S/m (320 mg/L) | 0.05 S/m (320 mg/L) |
| 18.0 | 34.0 | Calivil Fm — AQUIFER 1 | aquifer | kaolinite | 0.30 | 0.08 | 0.08 S/m (512 mg/L) | 0.44 S/m (2800 mg/L)  **changed** |
| 34.0 | 46.0 | clay aquitard — CONFINING | clay aquitard | mixed | 0.42 | 0.65 | 0.30 S/m (1920 mg/L) | 0.30 S/m (1920 mg/L) |
| 46.0 | 62.0 | Renmark Gp — AQUIFER 2 | lower aquifer | kaolinite | 0.28 | 0.10 | 0.60 S/m (3840 mg/L) | 0.60 S/m (3840 mg/L) |
| 62.0 | 78.0 | weathered Palaeozoic | weathered bedrock | kaolinite | 0.12 | 0.30 | 1.00 S/m (6400 mg/L) | 1.00 S/m (6400 mg/L) |
| 78.0 | 85.0 | saprock | saprock | kaolinite | 0.06 | 0.12 | 1.20 S/m (7680 mg/L) | 1.20 S/m (7680 mg/L) |
| 85.0 | 199.9 | unweathered Palaeozoic | unweathered bedrock | kaolinite | 0.02 | 0.03 | 1.50 S/m (9600 mg/L) | 1.50 S/m (9600 mg/L) |

## True window values, 17.5-33.5 m (L = 16.000 m)

| | 2Pv2 | 2Sv2 |
|---|---|---|
| window conductance | 0.7449 S | 1.3962 S |
| water, total | 4.6956 m (0.2935 v/v) | 4.6956 m (0.2935 v/v) |
| water, SNMR | 3.7527 m (0.2345 v/v) | 3.7527 m (0.2345 v/v) |
| pore water | 561 mg/L | 2,652 mg/L |
| clay fraction | 0.1355 | 0.1355 |
| clay share of conductivity | 76.1% | 40.7% |
| conductance 0-17.5 m | 6.1147 S | 6.1477 S |
| conductance 33.5-60 m | 7.9235 S | 8.0124 S |

## Names

ANYTHING THAT LEAVES THE REPO FOR THE PAPER -- figures, report pages, tables --
uses the general name only. The formation names in the `unit (internal)` column
are for tracing a result back to its config and must not appear in figure
output. Ruled 21 Nov 2026.

The deepest unit, 85 m and below, had no general name in the ruling's list,
which covers seven of the eight units. It is called `unweathered bedrock` here,
beneath `saprock`, by the same logic as `weathered bedrock` above it. Confirm
or replace it.


## For the ceiling figures: apply the SNMR lower limit on porosity

Added 22 Nov 2026, on Neil's instruction.

The salinity ceilings as published draw porosity from `[priors].porosity_range`
in `config.toml`, which is `[0.25, 0.40]`. `stages/stage5_quality.jl:98` reads
that range directly and nothing in the file refers to the SNMR water content,
so the ceilings are computed over porosities the SNMR data already exclude.

**The ceiling figures need the lower limit applied.** Porosity then runs from

    max(0.25, pooled SNMR p05 over the window)   to   0.40

The upper end does not change, so the least cautious ceiling is unchanged; the
central and most cautious ones fall. Values of record, on the shared window
convention in `src/window.jl`:

| column | pooled SNMR p05 | lower limit | most cautious ceiling, before -> after |
|---|---|---|---|
| 3P | 0.2760 | 0.2760 | 2,406 -> **1,994** mg/L |
| 3S | 0.2601 | 0.2601 | 9,998 -> **9,273** mg/L |
| 2S | 0.2332 | 0.25 (the config floor wins) | unchanged |

2S is left out of the ceiling table entirely: both its band edges are open, so
its "upper" figure is the top of the solved grid and a bound, not a ceiling.

2Pv2 and 2Sv2 follow the same rule once their Stage 2 runs land; their numbers
will be added here.

The window mean itself must use `src/window.jl`: theta is piecewise-constant
over 60 LAYERS whose bottoms are `zsmr`, the first running 0 to 1.0 m. Reading
`zsmr` as 60 cell boundaries instead moves the 3P window mean from 0.2811 to
0.2786.



## tools/naive_probability.jl -- the standard-practice contrast

Added 22 Nov 2026. **Discussion only. Not a workflow stage** -- nothing reads
its output and no stage depends on it.

It shows what the usual approach -- invert the AEM, assume rock properties,
divide Archie out -- would have concluded from the conductance alone, so the
paper can say what the joint method adds. It is read-only: it assembles
existing Stage 1 draws and runs no inversion.

    julia --project=. tools/naive_probability.jl --cols 3P,3S,2P,2S
    julia --project=. tools/naive_probability.jl --cols 2Pv2,2Sv2 --sampling both

`--sampling` is `loguniform` (the configured prior, the default), `uniform`,
or `both`; `--scales` multiplies the clay upper limit, default 0.5, 1, 2.

Results of record, log-uniform clay prior, pooled over four runs:

| column | fresh | brackish | saline | median |
|---|---|---|---|---|
| 3P | 100.0% | 0.0% | 0.0% | 0 mg/L |
| 3S | 38.8% | 39.1% | 22.1% | 2,033 mg/L |
| 2P | 36.6% | 14.0% | 49.3% | 2,938 mg/L |
| 2S | 30.8% | 11.5% | 57.6% | 3,818 mg/L |

Two things the discussion draws from it. **2P's saline share runs from 38.7%
to 63.1% across four identical runs**, which is the Stage 1 convergence failure
showing through (window conductance R-hat 1.195, ESS 23); a single AEM run
would have landed anywhere in that range. And the fresh and brackish partners
are barely separated -- 30% apart in median salinity where the true pore waters
differ by 14x and 5x.

**If a figure is wanted**, the natural one is a stacked bar per column with the
three class shares, the per-run spread as an overlay on 2P, and the uniform
result beside it as the contrast. The numbers above are all it needs; the
script prints them for any column that has Stage 1 draws on disk.


## Files

- `data/truth_murray2Pv2.jld2`
- `data/truth_murray2Sv2.jld2`
- `configs/synth_2Pv2.toml`, `configs/synth_2Sv2.toml` (site configs)
- `figures/fig_truth_all.jl` draws all six columns
