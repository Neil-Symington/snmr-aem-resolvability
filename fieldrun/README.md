# Batch field runs

The same five stages the synthetic work uses, driven from one control file over
a list of real sites. `example/` holds a two-site extract from Geoscience
Australia's Upper Darling Floodplain surveys so the workflow can be run and
checked end to end — see `example/PROVENANCE.md`, which also states the release
position on that data.

## Status

This is a branch under construction. What is here now:

- `example/` — the two-site input set, 1.59 MB, with `SHA256SUMS`
- `config/udf_example.toml` — the control file, every stage parameter in one place
- `tools/assemble_example.py` — how `example/` was built, for provenance

- `common/fieldconfig.jl` — the control-file reader. Its failure paths are
  negative-tested by `fieldrun/test/test_fieldconfig.jl`, which the repository's
  `test/runtests.jl` runs: 76 cases, including six that assert a VALID config is
  still accepted. An earlier version of this line claimed the tests existed when
  nothing was committed — they were reviewers' scripts, on reviewers' machines.
- `stages/udf_system.jl` — the SkyTEM system, parsed from the shipped `.des`
  rather than transcribed, and cross-checked against what the control file
  declares

Still to come, in this order:

1. The remaining field stage scripts — AEM Occam, AEM trans-D, SNMR trans-D,
   detection, the profile likelihood and the water-quality gate — ported from
   `udf_realdata/` with their absolute paths lifted into `[paths]`. Three known
   blockers are listed under "Porting notes" below.
2. `run_batch.sh` — the driver, replacing `run_seven.sh`
3. `reference/` — expected verdicts and conductance quantiles for the two
   sites, so a full run can be checked rather than admired

Until 1 and 2 land there is no end-to-end run: what is here reads the control
file and builds the AEM system, and nothing else.

## How it will be driven

```bash
./run_batch.sh --config config/udf_example.toml
```

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
