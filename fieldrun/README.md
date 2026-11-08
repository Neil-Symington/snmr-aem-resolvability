# Batch field runs

The same six stages the synthetic work uses, driven from one control file over
a list of real sites. `example/` holds a two-site extract from Geoscience
Australia's Upper Darling Floodplain surveys so the workflow can be run and
checked end to end — see `example/PROVENANCE.md`, which also states the release
position on that data.

## Status

All six field stages are now in the repository, with their support modules:

```
common/   fieldconfig.jl      the control file, read once and validated
          stage_prelude.jl    what every stage does first (see below)
          aem_lib.jl  occam_udf.jl  rockphys.jl  transform.jl
          udf_smr.jl  prerun_data.jl  prerun_lib.jl
stages/   udf_system.jl       the SkyTEM system, parsed from the shipped .des
          run_aem_occam.jl    AEM smooth inversion
          run_aem_transd.jl   AEM trans-D
          run_smr_transd.jl   SNMR trans-D on Stage-1 conductivity draws
          detect_aquifer.jl   detection
          run_aem_profile.jl  the profile likelihood
          gate_water_quality.jl  the water-quality gate
run_batch.sh                  the driver
```

**What is verified, and what is not.** Every support module and four of the six
stages LOAD from a clone, and `run_aem_occam.jl` inverts a real sounding. A full
two-site run has not been completed, and `reference/` — the expected verdicts
for UD21 and UD24, so a run can be checked rather than admired — does not exist
yet. Do not quote a number from this branch.

**Three porting blockers, cleared.** Each had been named in this file for three
rounds before the stages moved:

1. `udf_smr.jl` hardcoded the SNMR archive root and an April-2024 metadata
   spreadsheet. Both now come from `[paths] snmr_dir` and `snmr_meta`.
2. `prerun_data.jl` set **four East Kimberley paths as `const`s at load time**,
   so a clone carried four dead absolute paths that no error mentioned, because
   they are used lazily. They now resolve from an optional `[paths] ek_root` and
   raise a message naming the cause if anything reaches for them. Nothing in the
   UDF stages reads them: that file is included for its generic machinery —
   `pick_coil`, `per_q_noise`, the sigma-hat convention, the common depth grid —
   which East Kimberley shares, so the two surveys' posteriors sit on one grid.
3. `export_udf_soundings.py` needs the 4.9 GB delivery and stays a provenance
   script, not part of the batch.

**Two more found by loading the files rather than reading them**, which is worth
recording because a syntax check passed both:

- `rockphys.jl` included `../hypotheses/transform.jl`, **outside the
  repository** — the same defect that made the synthetic Stage 2 unrunnable from
  a clone for two rounds. Vendored as `common/transform.jl`.
- **`MAT` was not a dependency.** The SNMR stages read Vista Clara `.mat` files
  and the package was undeclared, so two support modules could not load at all.
  It is now in `Project.toml`, and it pulls HDF5 and about 180 Manifest lines
  with it — so any install evidence from before it was added does not cover this
  environment.

**`stage_prelude.jl`** is what every stage includes first. `udf_system.jl` parses
the SkyTEM `.des` at load time and takes its path from the control file, so a
stage that includes a support module before loading the config fails with
`UndefVarError: FieldConfig` from inside a file that names neither the cause nor
the fix. The prelude also resolves the output and log directories, so a stage
writes where the control file says rather than beside itself.

Still to come: a completed two-site run, and `reference/`.

## How it will be driven

```bash
./fieldrun/run_batch.sh --config fieldrun/config/udf_example.toml
```

Run it from the repository root: that is the Julia project. The driver reads the
control file through the same loader the stages use, so the shell and Julia
cannot disagree about which sites, how many at once, or which stages to run.

`mode = "smoke"` in `[run]` cuts the chains to 2,001 samples on 2 chains with a
single kernel draw. It is what a reader runs first: it exercises every stage and
every file path and proves nothing about a posterior. `mode = "full"` is the
committed parameter set.

**No CI job runs any fieldrun stage**, because the stages do not exist yet. An
earlier version of this file said smoke mode "is what CI runs", which was never
true — a review caught it.

Stage order matters and the driver enforces it: detection must precede the
profile likelihood, because the profile reads the detected interval. An earlier
script had these the wrong way round and the profile silently skipped with "no
detection file".

Two sites at a time is the measured ceiling on a 31 GB machine. One site is an
AEM master plus tempered chains at about 6 GB, then an SNMR master plus chains
at about 6.5 GB. Four concurrent exhausted RAM and swap and per-process CPU fell
to 12–20%, so `npar = 2` and the driver waits on free memory rather than
thrashing.

## Why one control file

The parameters now in `config/udf_example.toml` were spread across four places:
constants in `udf_realdata/udf_system.jl`, constants in
`ek_prerun/prerun_config.jl`, argv in `run_seven.sh`, and literals inside the
stage scripts. Two runs could differ in a way no log recorded. The control file
is the unit of provenance: a result is quotable only alongside the file that
produced it.

## Porting notes

Three things must change as the stage scripts move, and none are cosmetic:

1. **Absolute paths.** `udf_smr.jl` hardcodes the SNMR archive root and an
   April-2024 metadata spreadsheet; `udf_system.jl` hardcodes the `.des`. All
   three become `[paths]` entries.
2. **`prerun_data.jl` sets four East Kimberley path constants at load time.**
   The field stages include it, so a clone would carry four dead paths. They are
   only used lazily, so nothing errors — which is exactly why it would survive
   unnoticed. They need lifting into config or isolating per survey.
3. **`export_udf_soundings.py` needs the 4.9 GB delivery** and derives its site
   list from a filter expression. It stays a provenance script, not part of the
   batch.

## Relationship to the synthetic workflow

Shared: the install, `Project.toml`, `Manifest.toml`, the depth grid, the
sigma-hat noise convention and the GP prior. That sharing is the point — East
Kimberley, Upper Darling and the synthetics then sit on the same grid with the
same noise model, and can be compared without an apples-to-oranges caveat.

Not shared: the stage scripts themselves, which read field inputs rather than a
sealed synthetic truth. There is no truth to compare against here, which is why
`reference/` records observed verdicts rather than errors against a known
answer, and why bore control does the work a truth model does in the synthetics.
