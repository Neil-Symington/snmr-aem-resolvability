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

## Files

- `data/truth_murray2Pv2.jld2`
- `data/truth_murray2Sv2.jld2`
- `configs/synth_2Pv2.toml`, `configs/synth_2Sv2.toml` (site configs)
- `figures/fig_truth_all.jl` draws all six columns
