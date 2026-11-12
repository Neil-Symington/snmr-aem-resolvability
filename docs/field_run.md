# Running the example field inversion

The field workflow was retired on 12 Nov 2026 and is **not in this tree**. It
lives in the tag `fieldrun-archive-20261112`, complete: six stages, the
control-file loader, the reference check, and a two-site example with real data.

These are the instructions for getting that example running again, written
because the code being archived is not a reason for the knowledge to go with it.
Everything below was measured on the machine the workflow was developed on, a
12-core box with 31 GB, during the one full-mode run that completed.

Read the whole page before starting. The example takes **1 hour 48 minutes**,
and three of the five things that can go wrong cost most of that before they
show themselves.

## 1. Restore it

```bash
git checkout fieldrun-archive-20261112 -- fieldrun
```

That brings back 40 files under `fieldrun/`, including the example data, which
is real: two SNMR soundings as Vista Clara `.mat` files, the AEM soundings and
system descriptor, a clipped nav GeoParquet, and the site and metadata tables.
Nothing else needs downloading. `fieldrun/example/SHA256SUMS` and
`PROVENANCE.md` record where each file came from.

Restoring puts the directory back in your working tree without moving the
branch, so `git status` will show it as new files. Do not commit it unless you
intend to un-retire the workflow.

## 2. Check the environment

`MAT` is still in `Project.toml` and the Manifest, so the `.mat` readers work
with no further install. If that has changed since, the restore will fail at
the first SNMR stage rather than at load.

Julia must be 1.12.x, which `install.jl` enforces as a range.

## 3. Run it

```bash
./fieldrun/run_batch.sh --config fieldrun/config/udf_example.toml
```

From the repository root, which is the Julia project. The driver reads the
control file through the same loader the stages use, so the shell and Julia
cannot disagree about which sites, how many at once, or which stages run.

The example runs `mode = "smoke"`: 2,001 samples on 2 chains with a single
kernel draw. Measured, two sites, `npar = 2`:

| stage | cost |
|---|---|
| `aem_occam` | ~4 min/site |
| `aem_transd` | ~2 min/site |
| `snmr_transd` | ~3 min/site |
| `detect` | ~3 min total |
| `profile` | **~95 min total** |
| `gate` | ~25 s total |
| | **1 h 48 min end to end** |

`mode = "full"` is the committed parameter set, 150,001 samples on 4 chains with
4 kernel draws, and adds roughly 13 min/site of SNMR trans-D on top.

## 4. Check the run, do not read it

```bash
julia --project=. fieldrun/reference/check_run.jl --config fieldrun/config/udf_example.toml
```

It compares the deterministic Occam numbers against stored values, asserts the
structural invariants downstream, prints one line per check and **exits
non-zero** on any disagreement.

This step is not optional, and it is the most important instruction on this
page. The first complete run of this workflow printed `complete` and six `ok`
lines over a run in which both SNMR sites had thrown and all three sweep stages
had skipped every site, producing no output files at all. A run that only prints
a table can only be admired.

## The five traps

**Smoke mode is not quick, and it is a different answer rather than a smaller
one.** `profile` is 95 of the 108 minutes and smoke does not shorten it at all,
because the λ² sweep and the bisection contain no Monte Carlo. Worse, at 2,001
samples the hierarchical noise scale runs to 1268× σ̂ against the full run's
7.4×, because `burninfrac = 0.5` discards only 1,000 samples. That halves the
constraint depth and flips SMR_UD21 from `AQUIFER above 43 m` to `NOT RESOLVED
above 21 m`. **No smoke posterior number should ever be quoted.**

**Worker count.** The SNMR stage's own usage line is `julia -p 4`, and
`run_batch.sh` launches plain `julia`. When those disagreed, `nworkers()` was 1,
`min(NCHAINS, nworkers())` silently ran one chain, and HiQGA rejected the
one-chain temperature ladder with `range(0.0, stop=0.3979, length=1): endpoints
differ`. The stage now calls `addprocs` itself and asserts the count, so a
mismatch fails immediately instead of running one chain in silence.

**Two sites at a time is the ceiling.** One site is an AEM master plus tempered
chains at about 6 GB, then an SNMR master plus chains at about 6.5 GB. Four
concurrent exhausted RAM and swap and per-process CPU fell to 12–20%. `npar = 2`
and `minfree_gb = 7` in `[run]` set it, and the driver waits on free memory
rather than thrashing. Do not raise `npar` on a 31 GB machine.

**`gates_pass false` is not a bad run.** The `[convergence]` gate pools R̂ and
ESS across ensembles that condition on *different* SMR kernels, one per draw, so
the pooled statistic is not measuring chain mixing. Within kernel the same runs
give split-half ≤ 0.074 and posterior-predictive χ²/N of 0.98–1.05 per ensemble;
pooled, that χ²/N reads 1.8–39. The gate fails every correct run. It is a
specification problem, not a defect in the run.

**No chain is seeded.** Neither the AEM nor the SNMR trans-D chains are seeded,
so two runs of one config differ: `sigma_draws` by up to 8.05 S/m between two
runs of the same stage. Every artefact's `seed` or `seeds` field describes
something other than the chain. At SMR_UD21 that spread spans a verdict
boundary, so a verdict from a single run is one draw from a distribution of
runs.

## What was never established

Worth knowing before you trust a result from it.

- The SNMR half of the port was **never verified** against the pre-port
  original, and cannot be verified by comparison while no chain is seeded: the
  port-versus-original difference was 9.33 S/m against 8.05 S/m between two runs
  of the port itself.
- 17 defects were found across four run attempts, and the discovery rate was
  not falling when the workflow was retired.
- It was never run from a clean clone on a second machine.

What *was* established: the Occam half is bit-identical to the pre-port
original, `max |Δm| = 0.000e+00` over all 70 cells with χ²/N agreeing to twelve
figures on both sites, asserted by `fieldrun/reference/udf_example.toml` at
1e-12 against a 1e-4 tolerance. One full-mode two-site run completed in
4 h 48 min and passed `check_run.jl`. One clean run is one data point.

## Reading further, without restoring anything

```bash
git show fieldrun-archive-20261112:fieldrun/README.md
git show fieldrun-archive-20261112:fieldrun/reference/README.md
git show fieldrun-archive-20261112:fieldrun/config/udf_example.toml
git show fieldrun-archive-20261112:fieldrun/example/PROVENANCE.md
```
