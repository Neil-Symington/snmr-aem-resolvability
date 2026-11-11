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

**What is verified, and what is not.** A full two-site run through all six
stages has been completed in smoke mode and checked: `reference/check_run.jl`
passes both sites and exits 0. The Occam inversion is **bit-identical** to the
pre-port code — `max |Δm| = 0.000e+00` over all 70 cells, χ²/N to twelve
figures, both sites — which is what `reference/udf_example.toml` asserts.

The SNMR half of the port has **not** been shown faithful: that needs a
full-mode run compared against the production values, which `reference/README.md`
records as the target. Until then, do not quote an SNMR number from this branch.
Nor should a smoke number be quoted at all — smoke mode is a different answer,
not a smaller one, and it flips one of the two detection verdicts.
`reference/README.md` has the measurement.

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

Both are now done: the two-site run completes and `reference/` checks it. Still
to come: a full-mode run, which is the test of the SNMR half of the port.

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
committed parameter set — 150,001 samples on 4 chains with 4 kernel draws, which
is what the UDF production run used.

**Smoke mode is not quick.** Measured end to end on two sites with `npar = 2`:
1 h 48 min, of which `profile` is 95 minutes. Smoke only shortens the two MCMC
stages; `profile` is an Occam solve plus a λ² sweep plus a bisection per band
edge, with no Monte Carlo in it to shorten, so it costs the same either way and
dominates both. The control file carries the per-stage breakdown.

**Then check the run, do not read it.**

```bash
julia --project=. fieldrun/reference/check_run.jl --config fieldrun/config/udf_example.toml
```

It compares the deterministic Occam numbers against stored values and asserts
the structural invariants downstream, exiting non-zero on any disagreement. It
exists because the first complete run of this workflow printed "complete" and
six "ok" lines over a run in which both SNMR sites had thrown and all three
sweep stages had skipped every site; `reference/README.md` has that story and
what came out of it.

**No CI job runs any fieldrun stage.** The stages exist and run now, but a
1 h 48 min smoke run does not belong in CI, and the field inputs it needs are
not all in the repository. An earlier version of this file said smoke mode "is
what CI runs", which was never true — a review caught it.

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

## Data preparation

Two Python tools in `tools/`, neither part of the batch. Both want a python
with `pyarrow`; the project's conda env does not have one, so
`python3 -m pip install pyarrow` first.

`aseg_to_geoparquet.py` converts a GA ASEG-GDF2 delivery (the `.dfn` layout
plus the fixed-width `.dat`) to one GeoParquet file with a WKB point geometry
column, streaming it so the 4.9 GB Upper Darling delivery converts in about
three minutes. `--schema` prints the parsed layout and writes nothing, which is
the quickest way to see a delivery's gate structure. `--selftest` checks the
reader against a built-in fixture and needs no delivery files.

It refuses rather than guesses in the two places that matter. If the record
width in the `.dfn` does not match the `.dat` it stops, because shifted
fixed-width numbers still parse and the output would look right and be wrong.
If the `.dfn` carries no `EPSGCODE` on the coordinate fields it stops and prints
what it did find, because the 2017 Ord-Keep delivery states its datum only as
"(GDA94)" in free text and a GeoParquet with a wrong CRS is silently wrong.
Checked against both deliveries: 195 columns from 39 DEFN entries at a 2788
character record for Upper Darling, 286 at 4184 with a leading `DATA` marker for
Ord-Keep.

`assemble_example.py` builds `example/` from the source archives, and is
provenance rather than something a reader runs. It does not subset the ASEG
delivery; `example/aem/nav_clip.parquet` was produced outside the repository,
which is the gap the converter above closes.

## Relationship to the synthetic workflow

Shared: the install, `Project.toml`, `Manifest.toml`, the depth grid, the
sigma-hat noise convention and the GP prior. That sharing is the point — East
Kimberley, Upper Darling and the synthetics then sit on the same grid with the
same noise model, and can be compared without an apples-to-oranges caveat.

Not shared: the stage scripts themselves, which read field inputs rather than a
sealed synthetic truth. There is no truth to compare against here, which is why
`reference/` records observed verdicts rather than errors against a known
answer, and why bore control does the work a truth model does in the synthetics.
