# The synthetic SNMR measurements are generated from an inversion result

**Reported, not fixed.** 18 Nov 2026. Stop condition (3b) of Neil's brief is met.

## What was asked

> a. The synthetic SNMR MEASUREMENTS are generated from the TRUE column: true
> water content, kernel computed on the TRUE conductivity, plus noise. Never
> from an inversion result.
> b. The INVERSION kernel is computed on the pooled Stage 1 median.
>
> If it only sets the inversion kernel, use the pooled median and proceed. If it
> generates the measurements from that conductivity, stop.

## What `tools/smr_data.jl` does

It generates the measurements from that conductivity. The water content is the
truth; **the conductivity is not**.

```julia
T = load(joinpath(OUT, "truth_murray$(TAG).jld2"))
w = T["θ_smrgrid"]                       # TRUE water content          <- correct
af = USE_TRANSD ? joinpath(OUT, "aem$(TAG)$(LABEL)", "aem_posterior.jld2") :
                  joinpath(OUT, TAG == "2" ? "aem_occam.jld2" : "aem_occam$(TAG).jld2")
σ50 = (USE_TRANSD ? A["σ50"] : A["σ_occam"])[keep]    # INVERSION sigma  <- the fault
ct = SNMRForward.ConductivityModel(σ50, tA)
F  = SNMRForward.MRSForward_square(R, zsmr, QGRID, ϕ, 0.0, Be, ct)
V0 = abs.(SNMRForward.forward(F, w))                  # the measurements
```

`V0` plus noise is written as `typical_V` and `poor_V`, which is exactly what
Stage 2 inverts. So the measurements are forward-modelled through a kernel built
on an **AEM inversion's** conductivity. The truth file carries `σ_bulk` and the
script never reads it.

The `--label` conductivity therefore does **not** only set an inversion kernel.
Under `transd` it selects which Stage 1 ensemble's p50 **generates the data**.

## Measured, not inferred

Every existing SNMR data file stores the conductivity it used as `σ_kernel`.
Compared with the true `σ_bulk` of the same column at the same depths:

| file | median \|diff\| from truth | max | verdict |
|---|---|---|---|---|
| `smr_data_murray3P.jld2` | 54.9% | 553% | not the truth |
| `smr_data_murray2P.jld2` | 54.3% | 419% | not the truth |
| `smr_data_murray3S.jld2` | 45.0% | 447% | not the truth |
| `smr_data_murray2S.jld2` | 51.6% | 402% | not the truth |

And compared bit-for-bit with the Occam inversions, `σ_kernel` is **identical**
(max absolute difference < 1e-12) to `σ_occam` from `aem_occam<TAG>.jld2` for all
four columns.

## The 13 Nov production run on 2P: yes, it was built that way

`data/smr_data_murray2P.jld2`, dated 13 Nov, carries the suffix `2P` and not
`2PT`. In `smr_data.jl`, `SUF = USE_TRANSD ? TAG*"T" : TAG`, so the absent `T`
means the **Occam branch** ran, and its `σ_kernel` is bit-identical to
`aem_occam2P.jld2`. The 13 Nov 2P SNMR measurements were generated through a
kernel built on **2P's Occam inversion conductivity**, not on 2P's true
conductivity. The same holds for 3P, 3S and 2S.

This matters more than a factor on a kernel. The SMR kernel is
conductivity-sensitive, and 2P is the column whose aquifer is a resistor
enclosed by conductors — precisely where an inversion's conductivity departs
most from the truth. The 54% median departure above is on the column the paper's
hard case rests on.

## Whether it is an inverse crime as well

Not automatically, and not under `--pool`. Stage 2 rebuilds its own kernel from
the conductivity it is given (`stage2_snmr.jl:178`,
`ConductivityModel(σ_aem, diff(z_aem))`) and never reads the `kernel` matrix
stored in the data file. So:

- Stage 2 **with** `--pool`: inversion kernel is the pooled Stage 1 median, data
  were generated from one Occam profile — a mismatch, but not a crime.
- Stage 2 **without** `--pool` on the same single profile that generated the
  data: generation and inversion operators coincide, which is an inverse crime
  and would flatter the SNMR result.

The 13 Nov runs predate `--pool`.

## Noise seeds, as asked

One seed, shared by every column. `smr_data.jl` has `SEED = 20261012` and
`NOISESEED = SEED` unless `--noiseseed` is passed. All four files record
`seed = 20261012` and `noiseseed = 20261012`.

- The per-moment noise **shape** (`sigma_shape`) is always drawn on `SEED`
  deliberately, being a survey property measured at East Kimberley.
- The noise **draw** (`add_noise`) uses `NOISESEED` for the typical tier and
  `NOISESEED + 7` for the poor tier.
- No file was written with `--noiseseed`, so **3P, 2P, 3S and 2S share one SNMR
  noise realisation**, as they already share one AEM realisation.

## A second fault, mine, found by trying to run Stage 2: it could not run at all

The timing probe failed immediately with

    ERROR: LoadError: cannot document the following expression:

`stage2_snmr.jl` had **two docstrings back to back**. The `pooled_sigma`
docstring I added in `a39a3e9` ("Stage 2 can pool the kernel across several
Stage 1 folders") was inserted between `inputsig`'s docstring (lines 182-193) and
`inputsig` itself (line 230), so the first docstring was documenting the second.

The file **parses** cleanly — `Meta.parseall` returns OK — and fails only at
load, which is why no syntax check caught it. It went unnoticed because
`tools/m4_run.sh` deliberately skips Stage 2 and `test/test_stage2_pool.jl`
exercises `pooled_sigma` without loading the stage. **Stage 2 has therefore been
unrunnable since `a39a3e9`**, which is also the commit whose test I reported as
passing — the test was sound, the stage around it was not.

Fixed by moving `inputsig`'s docstring back next to `inputsig`, with the reason
recorded in place. This one is repair of my own regression in already-agreed
work, not a change of scope; Stage 2 now loads and completes.

## Stage 2 timing, as asked

Measured, then scaled. A 2,001-sample run of all four ensembles took **95 s**
wall clock, of which the sampling was **36 s** and the rest Julia start-up and
compilation. Sampling is linear in samples, so 100,001 samples scales to about
**30 minutes**, plus a minute of start-up and a minute or two for `--pool` to
re-assemble roughly 8,000 post-burn-in draws from four 14 MB chain sets.

**Call it 30-35 minutes per column**, four workers, one column at a time. Two
columns concurrently would fit the box comfortably — it is a far lighter job than
Stage 1, which took 160 minutes an ensemble.

## What I have not done

No fix, no re-generation, nothing run on the strength of it. Stage 2 is not
started on 3S or 2S. The one thing I did run is a 2,001-sample **timing probe**
on the existing Occam path, purely to measure the per-sample rate; its output
`data/cond_C2_typical3S.jld2` is **scientifically void** — built on the very data
this report questions — and should be deleted. I tried to delete it in the same
command and the action was refused as destructive, so it is still there.

## The shape of a fix, for your decision only

`smr_data.jl` would read `σ_bulk` from the truth it already loads, interpolate it
to the AEM layering it already builds, and use that for the generating kernel —
keeping the inversion kernel exactly where it is, as a separate concern fed by
`--pool`. That makes the generating conductivity independent of any inversion and
removes the `--label` dependence from data generation entirely. It also means
every SNMR result in the project to date was inverted against data that no longer
exist in that form.
