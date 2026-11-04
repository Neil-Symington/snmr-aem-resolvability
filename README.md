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

### Install with `install.jl`, not with `Pkg.instantiate()`

```bash
julia --project=. install.jl
```

`Pkg.instantiate()` on its own **does not work from a clean Julia depot**, and the
error it prints points at this repository rather than at the cause. `HiQGA` depends
on `PyCall`, `PyPlot` and `Conda`; its build runs `Conda.add("matplotlib")`
unconditionally and then overwrites `ENV["PYTHON"]`, so setting `PYTHON` yourself has
no effect. `Conda.jl` bootstraps the current Miniconda, which ships `python3.x` and no
`python` alias, while conda's own launcher begins `#!/usr/bin/env python` — and on
Debian or Ubuntu there is no `python` on `PATH`, only `python3`. The build dies with
exit 127 and takes the whole environment with it.

`install.jl` creates the missing alias inside conda's own bin, puts that bin on `PATH`
for the build, and finishes. It was verified from an empty depot; if your machine
already has a working conda it is a no-op and takes seconds.

`Manifest.toml` is tracked, so the versions you get are the ones used here rather than
a fresh resolve.

## Quick start

```bash
julia --project=. run_site.jl configs/synth_2P.toml --dry-run   # validate, run nothing
julia --project=. run_site.jl configs/synth_2P.toml --from 3    # the fast half
```

`--from 3` is the usual entry point. Stages 0 and 3–6 take minutes. **Stages 1 and 2
take hours** — each is a trans-dimensional MCMC of 100,000–200,000 samples over four
chains — so they run from `tools/` rather than from the driver, and `run_site.jl`
prints the command instead of starting a job that will outlive your terminal.

### Where the synthetic columns come from

Generation lives in a separate repository, **[snmr-aem-synthetics](https://github.com/Neil-Symington/snmr-aem-synthetics)** — a
config-driven generator with a column per TOML file:

```bash
git clone https://github.com/Neil-Symington/snmr-aem-synthetics
cd snmr-aem-synthetics
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/build_truth.jl config/confined_aquifer.toml
```

It is deterministic: the same config and `seed` give a bit-identical column, so a
stated parameter set reproduces a stated column exactly. The shipped configs are not
the parameter sets behind this repository's results — they differ in the porosity and
clay variances, the correlation length and the clay exponent — so building from them
demonstrates the method and produces a statistically equivalent column, not the one in
the paper.

### Reproducing the paper's numbers

Stages 1 and 2 are hours of MCMC each. Their products — the truth columns, the SNMR
data, the AEM posterior draws and the Stage 2 ensembles — are **21 MB in total** and
are published with the release. With them in `data/`, Stages 0 and 3–6 run in about an
hour and reproduce the verdict tables directly:

```bash
julia --project=. run_site.jl configs/synth_2P.toml --from 3
```

Without them you can still run `--dry-run` and the tests, but the first real stage is a
3.5-hour trans-dimensional inversion per column.

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
