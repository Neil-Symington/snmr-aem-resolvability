# Fast test: does the workflow depend on the trans-D sampler?

Branch `occam-workflow-test`, cut from `sampler-convergence-test` so it inherits
the corrected ESS estimator and the diagnostics. Runs alongside B1. Nothing
writes to `data/aem2P/`, `data/cond_*`, `outputs/` or `numbers.csv`.

## The brief's last question first: SNMR Occam already exists

**Both halves are already implemented**, in `tools/quicklook.jl`. Nothing needs
writing from scratch.

| step | code | status |
|---|---|---|
| 1. AEM Occam to chi2 = N | `src/occam_core.jl` via `tools/quicklook.jl:103-104` | **exists**, sound |
| 2. SNMR Occam of water content | `SMRPInversion.deterministic_invert` via `tools/quicklook.jl:158-163` | **exists**, sound |

Step 1 is exactly the brief's convention: `occam(...; mode = :minimise)` finds
the best misfit the forward model can reach, then `mode = :target` returns the
smoothest model still reaching chi2 <= N. Step 2 builds the kernel on the step 1
conductivity (`σ_kernel` from `smr_data.jl`, which reads `aem_occam<TAG>.jld2`)
and inverts water content with `regtype = :R2`, `λ²` searched over 1e-1 to 1e6,
`β² = 0.01`. Measured on the smoke column it reached chi2/N 1.000.

Two changes are needed, both small:

1. **quicklook.jl prints the SNMR profile but never saves it.** It writes
   `aem_occam<TAG>.jld2` for the AEM side and nothing for the water content, so
   steps 3-5 have no input. About 20 lines to save `θ_occam`, `zsmr`, chi2 and
   lambda-squared alongside it.
2. **It is a pass/fail check, not a producer.** It prints RMS against truth and
   returns a boolean. For this brief it has to emit an artefact per column.

Nothing in `occam_core.jl` needed changing. One thing worth recording: it cannot
use `build_skytem_operator`, because that lacks `calcjacobian` and the
Gauss-Newton steps need it. `build_aem_jac` from `src/aem_jac.jl` is the right
builder and quicklook already uses it.

## Steps 3-5 are the real work, and it is adapter work

Stages 3, 4, 5 and 6 all consume posterior **ensembles**: `stage3` reads a cond
file, `stage4_floor.jl:43` does `reduce(hcat, E["ens"])`, `stage5` reads a
`sigma` matrix of draws. An Occam result is one model.

The cheapest faithful route, and the one the brief's wording asks for ("apply
the existing rules", "run the existing gate"), is **not** to modify the stages.
It is to write ensemble-shaped files holding a single draw:

- `cond_C2_occam<TAG>.jld2` with `ens = [reshape(θ_occam, nz, 1)]`
- an AEM draws file with `sigma` as an `nz x 1` matrix

then run stages 3, 4, 5 and 6 unmodified. The percentiles they compute collapse
onto the single model, which is the honest representation of a deterministic
inversion: no posterior, so p05 = p50 = p95.

Known risks, to be checked rather than assumed:

- Stage 5 subsamples 600 draws; with one draw it may error rather than degrade.
- Stage 4's gate is a constrained optimisation (`occam_core.jl` `mode = :map`,
  `ctype = :cond`), which is already deterministic. It consumes Stage 1 draws
  only to set its starting model, which is exactly what the brief's last
  question asks about.
- Anything that divides by a band width will divide by zero on a single draw.

Estimate: half a day to a day, most of it finding which stages object to a
one-draw ensemble. The Occam solves themselves are minutes per column.

## Which four columns, and a correction to my own work

The brief says all four synthetic columns. The repository defines six site
configs, so the four are taken to be the two tests and their saline controls:

| tag | column | truth config | status |
|---|---|---|---|
| 2P | test two, confined Calivil/Renmark | `confined_aquifer.toml` | existed |
| 3P | **test one**, unconfined Loxton-Parilla | `paper_3P.toml` | built |
| 2S | test two saline control | `confined_aquifer_saline.toml` | **written, built** |
| 3S | test one saline control | `paper_3P_saline.toml` | **written, built** |

`2W` and `3W` are a separate lineage ("bare lineage, 88.6 m loop") and are left
out. Say the word and they are two more runs.

**A correction.** In the Part A work I built test one from
`synthetics/config/unconfined_aquifer.toml` and tagged it `1P`. That was wrong.
`configs/synth_3P.toml` labels itself "test one, unconfined Loxton-Parilla
sands, tier P", so **test one is 3P**. The two columns share a geometry but are
not the same column: `p_clay` is 3.0 in `paper_3P.toml` against 1.5 in
`unconfined_aquifer.toml`, the variability differs (3P carries
`log_sd_sigma_w = 0.25` and logit jitter, which the other lacks), and
`unconfined_aquifer.toml` still uses the retired `clay_mineral` dropdown rather
than a free `sigma_clay`. Measured, sigma_bulk differs by up to 0.79 S/m and
theta_SMR by 0.18 v/v.

B1's second leg was repointed from 1P to 3P before it started, so no compute was
wasted, and `data/truth_murray1P.jld2` is deleted. The Part A report needs its
test-one rows relabelled; its 2P findings are unaffected.

### The saline controls had no truth builder

There is no saline truth config in the repository and no code that derives one;
`parent` and `role` in `configs/synth_2S.toml` are metadata the stages carry
through, nothing more. The old truths in
`synthetics/outputs/HYP_SYNTH/truth_murray2S.jld2` are from the previous
workflow and in an older format (74 kB against 99 kB today), so they are not
usable directly, and they do not store their config.

The recipe was recovered by differencing the old 2P and 2S water conductivity
per unit:

```
Shepparton Fm, clay cap        0.049 -> 0.061   x1.25
Calivil Fm — AQUIFER 1         0.094 -> 0.297   x3.15
clay aquitard — CONFINING      0.329 -> 0.346   x1.05
everything else                        unchanged x1.00
```

So the saline control is a water-only change **in the target aquifer only**; the
1.25x and 1.05x on its neighbours are transition-zone smoothing bleeding across
the boundaries, not deliberate. The old files carry `σW_SALINE = 0.3`, and
0.297 S/m is that value after smoothing. The new configs therefore set the
target aquifer's `sigma_w` to **0.30 S/m**, about 1,920 mg/L at
`MGL_PER_SM = 6400`, from 0.08 (about 520 mg/L) in 2P and 0.03 in 3P.

That it straddles the 1,500 mg/L verdict threshold is the point of the control.

## The three extra questions

**Edge blur.** Occam regularises with `:R2`, so its interval edges will be
smoother than the truth's step and smoother than the sampler's, which can place
a discontinuity. Measured as the depth interval over which water content crosses
from 10% to 90% of its in-interval mean, at the top and base of each interval,
for truth, sampler and Occam.

**Sensitivity to the target misfit.** Steps 1-5 repeated on test two with the
target loosened 10% (chi2 <= 1.1 N). One extra line in the driver, because
`occam()` already takes `ndata` as the target.

**Does the gate depend on its starting model?** Run the Stage 4 gate on test two
twice, once started from the Occam conductivity and once from the sampler's
median, and compare admissible range, edge status, bound and floor. The gate is
a constrained optimisation and ought to be start-independent; if it is not, that
is a finding in its own right and it bears on the production numbers, not only
on this test.

## Reporting

One table per column, sampler against Occam against truth: interval top and
base, interval mean water content, window and package conductance of the
starting model with the depth range stated at a unit boundary, admissible
conductance range and edge status, bound, floor, verdict at 1,500 and 3,000
mg/L, run time. Factors throughout, conductivity not resistivity, no tuning
toward truth.

Package depth ranges end on unit boundaries: test two 4-34 m (clay cap plus
AQUIFER 1, which is where truth is 5.60 S) and test one 6-40 m (the aquifer,
which is unconfined so there is no cap).
