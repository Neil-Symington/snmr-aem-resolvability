# Which columns are which — 14 Nov 2026

Hashes are Julia's 64-bit `hash(read(path))`, the same function
`stages/stage0_inputs.jl:219` uses for `truth_file_sha`, so they can be
compared directly with what any stage0 artefact recorded. SHA-256 prefixes are
given too where useful.

## 1. B1 ran 3P for the unconfined case, and 3P is the paper's test one

**Evidence that it was 3P, not 1P.** No `data/aem1P_*` directory exists. The
running processes are `tools/aem_HYP.jl 200001 3P 4 _b1e3` and `_b1e4`, and
`convergence_test/logs/driver_b1_all.log` records only two legs, `tag=2P` and
`tag=3P`. The `1P` tag existed for about three minutes on 13 Nov as a
2,001-sample validation run, which was deleted along with
`data/truth_murray1P.jld2` before the 3P leg started.

**Evidence that 3P is the paper's column.** The reviewer's two candidate
package conductances separate them cleanly over 6–40 m:

| file | julia hash | bytes | package 6–40 m |
|---|---|---|---|
| `data/truth_murray3P.jld2` (B1 used this) | `d9de86c430ab0165` | 94,934 | **0.9142 S** |
| `synthetics/output/paper_3P_truth.jld2` | `d9de86c430ab0165` | 94,934 | **0.9142 S** |
| `synthetics/outputs/HYP_SYNTH/truth_murray3P.jld2` (the paper's) | `1cb14adbbb777aea` | 73,622 | **0.9142 S** |
| `synthetics/output/unconfined_aquifer_truth.jld2` | `bbf82033e1d06991` | 94,366 | 0.5806 S |

`outputs/synth_3P/stage0_synth_3P.jld2` records `truth_file_sha =
1cb14adbbb777aea`, which is the HYP_SYNTH file. It and the file B1 used are
**the same physical column**: `σ_bulk`, `z`, `θ_SMR`, `clay`, `σw` and `φ` agree
to at most 5.6e-17. The newer file differs only by carrying ten extra keys
(`config`, `p_clay`, `σ_clay`, `σ_water`, `petro_model` and so on) that the
older build did not write.

**So B1's test one is correct.** The 0.5806 S column is
`unconfined_aquifer.toml` and B1 never touched it.

## 2. 2P is NOT confirmed — the release repo builds a different column from the paper

This is the answer to the second question and it is not the one I expected.

| | the paper's 2P | the release's 2P |
|---|---|---|
| site config | `<parent>/configs/synth_2P.toml` | `release/configs/synth_2P.toml` (sha256 `06515d16…`) |
| truth path it declares | `synthetics/outputs/HYP_SYNTH/truth_murray2P.jld2` | `data/truth_murray2P.jld2` |
| truth julia hash | `6421d57fb6293968` | `33bf02a46442711e` |
| built from | (pre-dates the release scenario files) | `synthetics/config/confined_aquifer.toml` (sha256 `efa2744e…`) via `run_production.sh:30` |
| **package 4–34 m** | **7.9868 S** | **5.6022 S** |
| **window 17.5–33.5 m** | **0.8727 S** | **0.6856 S** |

`outputs/synth_2P/stage0_synth_2P.jld2` records `truth_file_sha =
6421d57fb6293968`, confirming the paper's run used the HYP_SYNTH file.

They are materially different columns, not two serialisations of one. The
paper's is 1.1–1.5x more conductive in almost every unit at similar clay
fractions, so the difference is in the water and petrophysics, not the geometry:

| unit | sigma_bulk release | paper | factor |
|---|---|---|---|
| Shepparton Fm, unsaturated | 0.3589 | 0.5348 | 1.49 |
| Shepparton Fm, clay cap | 0.3498 | 0.5094 | 1.46 |
| Calivil Fm — AQUIFER 1 | 0.0440 | 0.0534 | 1.21 |
| clay aquitard — CONFINING | 0.5773 | 0.6478 | 1.12 |
| Renmark Gp — AQUIFER 2 | 0.0921 | 0.1129 | 1.23 |

### This retracts my conductance reconciliation

On 13 Nov I reported that truth over 4–34 m is 5.60 S, that the stated 8.0 S
therefore implied integration to about 38–40 m, and that the paper's package
figure was a depth-range labelling error. **That was wrong.** The depth range is
right; the truth file is different. The paper's own column gives **7.9868 S over
exactly 4–34 m**, which is the stated 8.0 S.

Read against its own truth, each run behaves the same way:

| run | truth | posterior | high by |
|---|---|---|---|
| the paper | 7.99 S | 9.1 S | 14% |
| reviewer's 13 Nov production | 5.60 S | 6.65 S | 19% |

So there was never a discrepancy to explain — only two columns. The window
finding does survive: 2.32 S against a 0.8727 S paper truth is 2.66x, and
2.07 S against a 0.6856 S release truth is 3.02x, so the AEM still overstates
the aquifer window by about three in both.

### What this contaminates

The reviewer's 13 Nov run used `run_production.sh`, so it is on the **release**
column. Everything I have run is too. That makes the following internally
consistent with each other but **not the paper's test two**:

- B1 2P, all four ensembles
- the Occam 2P comparison
- 2S, which I built from `confined_aquifer_saline.toml`, same lineage
- the 2P and 2S rows of the conditions table

3P and 3S are unaffected.

The incoming `synth_2P_outputs_13Nov.tar.gz` will also be on the release column,
so redoing the Occam comparison against it is still valid as a
**sampler-against-Occam** test — both sides will be the same column. It just
cannot be labelled the paper's test two until the column question is settled.

### The decision

Either `synthetics/config/confined_aquifer.toml` is wrong and should be made to
reproduce the paper's 2P, or the paper's 2P is superseded and the release column
is the intended one. I cannot tell which from the repository: the release
scenario files post-date the HYP_SYNTH truths and carry no note about the
change. **This needs your call**, and if it is the former, B1 2P and the Occam
2P work need re-running on the corrected column.

## Convention from here

Every report names the config file and its hash. For a truth file the hash is
the same `filesha` the stages record, so it can be matched against any
`stage0_*.jld2`:

```julia
filesha(p) = bytes2hex(reinterpret(UInt8, [hash(read(p))]))
```

Current values:

| file | julia hash | sha256 (16) |
|---|---|---|
| `configs/synth_2P.toml` | `4399631518caef6f` | `06515d166b62b8c1` |
| `configs/synth_3P.toml` | `d2ab0d30d324298d` | `111440f92ec6adc7` |
| `synthetics/config/confined_aquifer.toml` | — | `efa2744e3e101a16` |
| `synthetics/config/paper_3P.toml` | — | `18bd80b0b237a2f6` |
| `data/truth_murray2P.jld2` | `33bf02a46442711e` | `46a1f0177fe2a75d` |
| `data/truth_murray3P.jld2` | `d9de86c430ab0165` | `d72dc63ef82d0a47` |
