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

**Julia 1.12 or newer.** `[sources]`, which pins an unregistered dependency
below, needs 1.11 — but the pinned `Manifest.toml` resolves `PrecompileTools`
1.3.4, which needs `Base.ReinferUtils` or `Base.StaticData`, and neither exists
before 1.12. On 1.11 `install.jl` fails at step 2. Install from
<https://julialang.org/downloads/> or with `juliaup`; the `Manifest.toml` was
resolved on 1.12.6.

**Budget about fifteen minutes and roughly an hour of CPU** for a first install on
a machine with no Julia depot. Most of it is Miniconda and matplotlib, pulled in by
`HiQGA`. Against a depot that already has them, a measured re-install took 37
seconds here and 46 seconds for the generator repository — under a minute, not
instant.

**The conda workaround is covered by CI through a long depot path.** The failure
needs a depot python path over 127 characters, so the CI job sets
`JULIA_DEPOT_PATH` to a deliberately long directory; a runner's default depot is
far too short to trigger it. Removing `python` from `PATH` — which an earlier
version of that job did — tests nothing, because a short path gets an absolute
shebang that never consults `PATH`.

**Do not run two installs at once on one machine.** `Conda.jl` downloads its
installer to a fixed `/tmp/installer.sh`, which is not per-depot, so concurrent
installs overwrite each other's file and one of them dies with "text file is
busy". A fresh `JULIA_DEPOT_PATH` does not isolate this.

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
pins the branch and `Pkg.instantiate()` will fetch it.

### Install with `install.jl`, not with `Pkg.instantiate()`

```bash
julia --project=. install.jl
```

`Pkg.instantiate()` on its own **does not reliably work from a clean Julia depot**,
and the error it prints points at this repository rather than at the cause.
`install.jl` expects that first attempt to fail and recovers from it. Step 1
prints the reason it stopped, indented under a `|`, and step 4 is what tells you
whether the environment is good. An earlier version of this paragraph promised a
red `Error building HiQGA` scrolling past; a reviewer found that the exception
was swallowed by a bare `catch`, so the log contained no occurrence of "rror" at
all. It does now. `HiQGA` depends
on `PyCall`, `PyPlot` and `Conda`; its build runs `Conda.add("matplotlib")`
unconditionally and then overwrites `ENV["PYTHON"]`, so setting `PYTHON` yourself has
no effect. conda chooses its launcher's shebang by **path length**, not by operating
system: if `#!<absolute path to conda's own python>` would exceed conda's
127-character limit it writes `#!/usr/bin/env python` instead, and `env` then resolves
`python` through `PATH`. Measured on one machine, everything else held fixed: a depot
python path of **125 characters gets the absolute shebang and works with no `python`
on `PATH` at all; 126 and above gets `env python` and fails** with exit 127.

`~/.julia` is about 57 characters, so a default depot never hits this — which is why
GitHub's runners install cleanly and a deep project-local or scratch depot does not.
A **system `python` does not rescue it** either: with `env python` resolving to
`/usr/bin/python3`, conda fails differently with `ModuleNotFoundError: No module named
'conda'`. Only conda's own python can run conda, so what `install.jl` does is not
"supply a python" but "put conda's python where `env` will find it".

**Putting conda's bin on `PATH` is the fix**, which is what `install.jl` does. It also
creates the `python` alias if it is missing, as a fallback for older Miniconda; on
26.9.4 the installer already ships it and that step is a no-op. It was verified from an empty depot; if your machine
already has a working conda it is a no-op and takes seconds.

`Manifest.toml` is tracked, so the versions you get are the ones used here rather than
a fresh resolve.

## Quick start

```bash
julia --project=. run_site.jl configs/smoke_2P.toml --dry-run   # validate, run nothing
```

**That is the only command that works on a fresh clone.** `--from 3` needs
Stage 1 and 2 products, which are not published, so on a clean tree it stops and
names what produces them. The second line used to be `--from 3` with no warning,
and two reviewers in a row hit it and got a path with no remedy.

To actually compute something, build a column and drive the chain yourself —
`configs/smoke_2P.toml` ships for exactly this, with a reduced gate grid:

```bash
julia --project=. synthetics/scripts/build_truth.jl synthetics/config/confined_aquifer.toml
mkdir -p data
cp synthetics/output/confined_aquifer_truth.jld2 data/truth_murray2P.jld2
julia --project=. tools/aem_HYP.jl 2001 2P 2 _smoke
julia --project=. tools/aem_draws.jl 2P 50 _smoke
julia --project=. tools/smr_data.jl transd --tag 2P --label _smoke
julia --project=. stages/stage2_snmr.jl --tier typical --conds C2 \
      --tag 2PT --label _smoke --nsamples 2001 --nchains 2
julia --project=. run_site.jl configs/smoke_2P.toml
```

That reaches a verdict table. It is **not a result** — 2001 samples and 50 draws
moved the Stage 5 TDS median by 50% between two machines — but every stage runs
and every file path is exercised. Reckon on 45 minutes to an hour and a half,
nearly all of it Stage 4.

### Running Stages 1 and 2 yourself

Neither stage has a `run_site.jl` entry point, and until a review pointed it out
nothing here named the scripts. The chain, for a tag `<T>` and an optional run
label `<L>`:

Stage 1 reads `data/truth_murray<T>.jld2`, and nothing automates getting a
generator-built column there under that name — a reviewer had to work that out
twice. Do it explicitly:

```bash
mkdir -p data
cp synthetics/output/<name>_truth.jld2 data/truth_murray<T>.jld2
```

Then:

```bash
julia --project=. tools/aem_HYP.jl 200001 <T> 4 <L>   # Stage 1, trans-D AEM
julia --project=. tools/aem_draws.jl <T> 600 <L>      # posterior draws
julia --project=. tools/smr_data.jl transd --tag <T> --label <L>
julia --project=. stages/stage2_snmr.jl --tier typical --conds C2 --tag <T>T
```

The `T` on the last line is not a typo: `smr_data.jl`'s trans-D branch writes a
`T` suffix so the Occam and trans-D kernels cannot collide. Stage 1 at 200,001
samples takes hours; 2001 samples on 2 chains finishes in minutes and exercises
every path.

**Run Stage 0 first, once.** A bare `run_site.jl` starts there; it re-seals the
inputs and takes seconds. Only then does `--from 3` work, because Stage 4 reads
Stage 0's sealed artefact — on a clean tree `--from 3` used to fail inside Stage
4 with a path and no explanation, and `run_site.jl` now refuses with the fix
instead. After that, `--from 3` is the entry point for every re-run: it skips
straight to the part worth repeating.

**Stages 1 and 2 take hours** — each is a trans-dimensional MCMC of
100,000–200,000 samples over four chains — so they run from `tools/` rather than
from the driver, and `run_site.jl` prints the command instead of starting a job
that will outlive your terminal.

**Stage 4 is not an MCMC, and `--nsamples` does nothing to it.** Its cost is an
Occam solve per conductance node, plus a λ² sweep at each node and a bisection
for the band edges. Measured, on one machine:

| `[gate]` grid | outcome |
|---|---|
| 16 coarse + 15 fine (production) | still in the pass-1 solve at 28 minutes |
| 4 coarse + 3 fine | reached the final edge bisection at 28 minutes |

So cutting the grid helps but **does not give Stage 4 a minutes-long mode** —
the cost is per node, not in the node count. There is no fast path here; budget
an hour or more, and reduce the grid only to shorten it:

```toml
[gate]
n_nodes_coarse = 4      # production: 16
n_nodes_fine   = 3      # production: 15
```

Stages 0 and 3 take seconds to minutes, and Stages 5 and 6 are quick once Stage
4 has run.

### Where the synthetic columns come from

Generation lives in **`synthetics/`**, which was the separate repository
`snmr-aem-synthetics` until the two were consolidated. There is no second
environment: one `install.jl`, one `Manifest.toml`, one `test/runtests.jl`.

```bash
julia --project=. synthetics/scripts/build_truth.jl synthetics/config/confined_aquifer.toml
julia --project=. synthetics/scripts/plot_truth.jl  synthetics/output/confined_aquifer_truth.jld2
julia --project=. synthetics/scripts/forward.jl     synthetics/output/confined_aquifer_truth.jld2
```

`synthetics/README.md` documents the config schema. The two were merged because
`src/aem_setup.jl` existed in both and diverged: the noise-seed fix landed in
one copy and not the other, so every published AEM number came from the wrong
realisation for weeks. There is now one copy, guarded by
`test/test_noise_seed.jl`, which builds the operator and compares the data rather
than matching the source text — three textual versions of that guard were each
defeated by a reviewer, sixteen times in total.

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

Reproduction across machines must be compared on **values, not file hashes**. A
truth `.jld2` embeds the Julia version string and a serialised `Dict`, so its
SHA-256 differs between machines and between Julia patch versions while every
number is identical — measured on three machines and two patch versions.

The config and noise-seed tests run anywhere. The four integration tests need stage
outputs and skip, loudly, without them.

## Citing and licence

See `CITATION.cff`. MIT; see `LICENSE`.
