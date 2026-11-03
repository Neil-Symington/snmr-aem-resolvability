# Can an AEM survey tell fresh groundwater from clay?

A seven-stage workflow that takes an airborne electromagnetic (AEM) sounding and a
surface nuclear magnetic resonance (SNMR) sounding and returns a **verdict** on a
pre-registered hypothesis about an aquifer: is there an interval of a stated
thickness, holding a stated amount of water, whose pore water is fresher than a
stated limit?

The answer is often "the data cannot tell", and the point of the workflow is to say
so defensibly rather than to produce a number. Two properties make that possible:

* **The conversion from conductance to salinity is clay-free**, so the salinity it
  reports is an *upper* bound. The data can therefore certify that water is fresh and
  can never certify that it is saline. A saline verdict comes only from declared
  ceilings on clay, porosity and the cementation exponent, and is labelled as such.
* **No stage after Stage 0 may read a truth value.** On the synthetic columns the
  firewall is a file boundary: Stage 0 seals a whitelisted artefact and every later
  stage reads only that. `test/test_stage0_roundtrip.jl` checks that the seal
  reproduces the forward operator exactly.

## Install

Julia 1.10 or newer.

```julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Three dependencies are geophysics packages rather than general ones:

| package | what it does | where it comes from |
|---|---|---|
| `HiQGA` | trans-dimensional MCMC, the SkyTEM forward, the joint operator | registered |
| `SNMRForward` | the SNMR kernel | registered |
| `SMRPInversion` | the SNMR inversion | **not registered, and not `main`** — pinned in `[sources]` to branch `add-occam-deterministic-inversion` |

`SMRPInversion`'s `main` branch does **not** carry the hierarchical heteroscedastic
noise model or the Occam bootstrap that every result here depends on. `Project.toml`
pins the branch and `Pkg.instantiate()` will fetch it; if you are on Julia 1.10,
`[sources]` is ignored and you must add it by hand.

`PyPlot` needs a working matplotlib. Only the figure scripts use it, so the pipeline
runs without one.

## Quick start

```bash
julia --project=. run_site.jl configs/synth_2P.toml --dry-run   # validate, run nothing
julia --project=. run_site.jl configs/synth_2P.toml --from 3    # the fast half
```

`--from 3` is the usual entry point. Stages 0 and 3–6 take minutes. **Stages 1 and 2
take hours** — each is a trans-dimensional MCMC of 100,000–200,000 samples over four
chains — so they run from `tools/` rather than from the driver, and `run_site.jl`
prints the command instead of starting a job that will outlive your terminal.

Nothing is shipped under `data/`. Rebuild the synthetic inputs first:

```bash
julia --project=. tools/truth_murray.jl      # the truth columns
julia --project=. tools/aem_HYP.jl 200001 2M # Stage 1, hours
julia --project=. tools/aem_draws.jl 2M      # the posterior draws Stage 5 uses
```

## The stages

| stage | what it does | cost |
|---|---|---|
| 0 | read the truth or the field data; seal a whitelisted artefact | seconds |
| 1 | trans-dimensional AEM conductivity posterior | ~3.5 h |
| 2 | SNMR water-content posterior on the Stage 1 median kernel | hours |
| 3 | the pre-registered interval and the water criterion | minutes |
| 4 | the admissible conductance band, the clay-free bound, the floor | ~1 h |
| 5 | pore-water salinity under the Stage 0 priors | minutes |
| 6 | the verdict table | seconds |

## The examples

`configs/synth_3P.toml` and `configs/synth_2P.toml` are two synthetic columns — an
unconfined sand and a confined aquifer pair. Each has a saline control (`3S`, `2S`)
that shares its parent's clay, porosity and water content bit for bit and differs
*only* in pore-water conductivity, so any difference between a column and its control
can be nothing but the water. A control names its parent in `[data].parent` and the
validator refuses one that does not: comparing a control against the wrong parent has
already produced one retracted finding in this project.

`configs/synth_2PJ.toml` adds a galvanic (Schlumberger) sounding as a second data term
in the same trans-dimensional chain, through HiQGA's `JointOperator`. HiQGA ships no
DC physics, so the forward is in `src/dc_forward.jl`: the Ghosh linear-filter method
on HiQGA's own Hankel filters, verified against the two-layer analytic series to a
worst relative error of 3.3e-5, against a 65-layer subdivision of it to 1.4e-13, and
against an independent adaptive-quadrature implementation to 6.5e-5.

Release of the field data is Geoscience Australia's decision, not this repository's.

## Configuration

Everything that is a choice lives in `config.toml`, which is commented at length and
is where to read why each prior is what it is. A per-site config inherits from it with
`extends` and restates only what differs, so the whole difference between two sites
fits on one screen. The sections are `[data]`, `[priors]`, `[hypothesis]`, `[gate]`,
`[inversion]`, `[ceilings]` and `[output]`.

`docs/choosing_priors.md` has the reasoning behind the lithology-keyed cementation
exponent, the porosity ranges, the clay prior and the ceilings.

## Tests

```bash
julia --project=. test/runtests.jl
```

The forward and config tests run anywhere. The three integration tests need stage
outputs and skip, loudly, without them.

## Citing and licence

See `CITATION.cff`. MIT; see `LICENSE`.
