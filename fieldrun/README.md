# Batch field runs

The same six stages the synthetic work uses, driven from one control file over
a list of real sites. Two AEM systems are supported: SkyTEM (the Upper Darling
Floodplain example) and the NRG XCite, a single-moment central-loop system.
`example/` holds a two-site extract from Geoscience Australia's Upper Darling
Floodplain surveys so the workflow can be run and checked end to end. See
`example/PROVENANCE.md`, which also states the release position on that data.

## Layout

```
common/   fieldconfig.jl      the control file, read once and validated
          stage_prelude.jl    what every stage does first (see below)
          smr_setup.jl        the SNMR trans-D settings, shared by the SNMR
                              stage and its report
          aem_lib.jl  occam_udf.jl  rockphys.jl  transform.jl
          udf_smr.jl  prerun_data.jl  prerun_lib.jl
prep/     prep_smr.jl         a GMR field delivery -> SNMR inputs
          prep_aem_xcite.jl   the nearest XCite sounding to each SMR site -> AEM inputs
stages/   udf_system.jl       the SkyTEM system, parsed from the shipped .des
          xcite_system.jl     the XCite system, parsed from the delivered .des
          run_aem_occam.jl    AEM smooth inversion
          run_aem_transd.jl   AEM trans-D
          run_smr_transd.jl   SNMR trans-D on AEM trans-D conductivity draws
          detect_aquifer.jl   detection
          run_aem_profile.jl  the profile likelihood
          gate_water_quality.jl  the water-quality gate
          report_smr.jl       standard SMR figures and tables (run by hand)
reference/                    the stored check of the example run
tools/                        Python data-preparation tools (not in the batch)
run_batch.sh                  the driver
```

## Running

```bash
./fieldrun/run_batch.sh --config fieldrun/config/udf_example.toml
```

Run it from the repository root, which is the Julia project. The driver reads
the control file through the same loader the stages use, so the shell and
Julia cannot disagree about which sites, how many at once, or which stages to
run. Stage order is enforced: detection must precede the profile likelihood,
because the profile reads the detected interval.

`mode = "smoke"` in `[run]` cuts the chains to 2,001 samples on 2 chains with a
single kernel draw. It exercises every stage and every file path and proves
nothing about a posterior. `mode = "full"` is the committed parameter set:
150,001 SNMR and 200,001 AEM samples, 4 chains, 4 kernel draws, as in the UDF
production run.

**Smoke mode is not quick.** Two sites with `npar = 2` took 1 h 48 min, of which
`profile` was 95 minutes. Smoke only shortens the two MCMC stages; `profile` is
an Occam solve plus a λ² sweep plus a bisection per band edge, with no Monte
Carlo in it to shorten. The control file carries the per-stage breakdown.

**Then check the run, do not read it.**

```bash
julia --project=. fieldrun/reference/check_run.jl --config fieldrun/config/udf_example.toml
```

It compares the deterministic Occam numbers against stored values and asserts
the structural invariants downstream, exiting non-zero on any disagreement. It
exists because the first complete run of this workflow printed "complete" and
six "ok" lines over a run in which both SNMR sites had thrown and all three
sweep stages had skipped every site; `reference/README.md` has the story.

**Figures and tables** for the SNMR result, per site, after the run:

```bash
julia --project=. fieldrun/stages/report_smr.jl --config <toml> [--site S]
```

It reads what `run_smr_transd.jl` left on disk and runs no MCMC. The posterior
and tables are pooled over kernel draws; the misfit and swarm plots are per
draw, and use replicate 1 only when there are several.

**Memory.** One site is an AEM master plus tempered chains at about 6 GB, then
an SNMR master plus chains at about 6.5 GB (more with a longer SNMR ladder).
Two sites at a time is the measured ceiling on a 31 GB machine; four exhausted
RAM and swap. So `npar = 2`, and the driver waits on free memory (`minfree_gb`)
rather than thrashing.

## What is verified

- **Occam is bit-identical** to the pre-port code on the example: `max |Δm| =
  0.000e+00` over all 70 cells and χ²/N to twelve figures, both sites. This is
  what `reference/udf_example.toml` asserts.
- **The SNMR half of the port has not been shown faithful** against the
  production values; `reference/README.md` records the target. Until then, do
  not quote an SNMR number from this workflow. Nor should a smoke number be
  quoted at all: smoke mode is a different answer, not a smaller one, and it
  flips one of the two detection verdicts.
- **The XCite path has no stored reference.** It has been run end to end on
  field data, but nothing in `reference/` checks it.
- **No CI job runs any fieldrun stage.** A 1 h 48 min smoke run does not belong
  in CI, and not all field inputs are in the repository. The control-file
  loader and `prerun_lib.jl` helpers are unit-tested (`test/runtests.jl`,
  `fieldrun/test/`).

## AEM systems

`[aem] system` picks the system module: `"skytem"` (the default) or `"xcite"`.
Each module exports the same API, so the AEM stages are written once. The two
systems need different `[aem]` keys, and a key belonging to the other system
is refused by name rather than ignored.

**XCite** (`stages/xcite_system.jl`, inverted with HiQGA's
`VTEM1DInversion.dBzdt`). The header of that file has the reasoning; in short:

- **Waveform:** `.des` column 5 (the Tx current), with the receiver coil and
  filters applied as two first-order low-pass corners, `[aem] lowpass_hz`.
  Column 6 already includes that path, but it has not decayed when the first
  gates open, and HiQGA's ramp convolution does not converge for gates inside a
  ramp. The corners are **fitted** so that column 5 through them reproduces
  column 6. They are an effective model of the receiver path, not the
  contractor's filter design, which is not available: a review item.
- **Gates:** the `.des` arithmetic window middles, all kept by default
  (`[aem] gates`). Some delivered windows overlap; that is kept and flagged.
- **Noise:** σᵢ = √((`nominal_unc`·dᵢ)² + `halt_noise_pv`ᵢ²), HiQGA's VTEM form.
  `[aem] use_ml = true` also fits the noise level in the trans-D (HiQGA
  `useML`), leaving the 3% + high-altitude model to set only the relative noise
  between gates. The fitted level is saved per draw.
- **Attitude:** pitch and roll are not modelled. `prep_aem_xcite.jl` refuses
  soundings tilted more than `[aem] max_tilt_deg`.
- **Prior:** `[aem] rho_bounds_log10` sets the trans-D resistivity bounds,
  because the UDF floodplain default is too conductive for other ground.

## SNMR trans-D: the ladder and its convergence check

The SNMR stage samples water content on each of `[kernel] n_draws` AEM
conductivity draws (one kernel per draw) and pools the posteriors, so AEM
uncertainty is marginalised rather than conditioned on.

**Its own ladder.** `[mcmc] nchains_snmr`, `nchainsatone_snmr` and `tmax_snmr`
set the SNMR parallel-tempering ladder separately from the AEM one. Each falls
back to the shared `nchains` / `tmax` and to one chain at T = 1, which
reproduces the UDF production run. A 6-chain ladder to Tmax 10 swapped at 0–1%
between T = 1 and the next rung in most draws on field data; 7 chains with 2 at T = 1 and
Tmax 2.5 (rung ratio about 1.2) swaps readily.

**Replicate ladders.** `[mcmc] replicates_snmr = R` runs R independent
ladders, each with its own seed, on every kernel draw. One (draw, replicate)
pair is a **unit**, with its own chain files. Units can run as separate
processes, several per node, and be combined afterwards:

```bash
julia --project=. fieldrun/stages/run_smr_transd.jl --config <toml> --site S --ens 1,2 --rep 3
julia --project=. fieldrun/stages/run_smr_transd.jl --config <toml> --site S --combine
```

Without `--ens/--rep/--combine` the stage runs every unit, then combines. With
R = 1 the file names are unchanged, so earlier runs stay readable.

**The gate** (`[convergence] rhat`, `min_ess`) is rank-normalised split-R̂ and
bulk ESS on gate scalars: θ at `[convergence] depths_m`, total water to 100 m,
log s, misfit and layer count. Which comparison it uses depends on what was run,
and is recorded as `gate_basis` in `smr_transd_<site>.jld2`:

| `gate_basis` | when | R̂ compares |
|---|---|---|
| `replicate` | R ≥ 2 | the R replicate ladders on each kernel draw, each read as its T = 1 sample in iteration order. Worst over draws; ESS summed over draws, worst scalar in each. |
| `within-kernel` | R = 1, more than one chain at T = 1 | streams by chain id within each draw |
| `across-ensemble` | otherwise | the kernel draws against each other |

The ladder must also be stationary (misfit trend) with at least one round trip
per chain. All three R̂ values are saved, whichever one gates.

**Read the chain-id R̂ as a diagnostic only.** Once swaps are frequent, each
chain id visits T = 1 in short, state-selected stretches, which inflates R̂: on
one field test it read 10–26 where replicate ladders on the same draws read
1.4–1.6. Only independent replicates test mixing properly.

**Known limitation, current sampler.** On that test (2 sites × 4 draws × 4
replicates, 150,001 samples), replicate R̂ was 1.37 and 1.56, so the gate
failed. Noise scale s, misfit and layer count are the slow quantities: s
drifts downward through the kept samples in most ladders. Scoring the first
third, two thirds and all of the kept samples, ESS grew roughly in proportion
but the worst R̂ did not fall, so longer chains alone will not pass the gate.
Water-content products reproduced closely across replicates (interval water
within about 3%).

## Preparing field inputs

**SNMR** (`prep/prep_smr.jl`): from a GMR delivery, finds each site's
`*_params.txt` and `.mat`, chooses one `.mat` per site, and writes
`snmr_metadata.csv` and `snmr/`. Both the delivery layout (`SiteInfo/`,
`Reprocessed_*/QC_FID*/`) and a flat copy are read. `[snmr] mat_preference`
decides between EMI-cancelled and uncancelled files. In the GMR deliveries
checked so far `_noemi` means **no** EMI cancellation; an earlier version of
this script read it the other way. What was chosen and what was passed
over are both printed.

**AEM, XCite** (`prep/prep_aem_xcite.jl`): the nearest record to each SMR site
that passes QC (all gates present, tilt, `max_sep_m`), one sounding per site,
streamed from the delivered ASEG `.dat` or its `.zip`. Writes the AEM
manifest, one CSV per site, and `sites.csv` (the SMR ↔ AEM pairing). A whole
block is one pass over a very large file (tens of minutes); a pre-cut line
subset takes seconds.

**Python tools** (`tools/`, not part of the batch; need `pyarrow`):

- `aseg_to_geoparquet.py` converts an ASEG-GDF2 delivery (`.dfn` + fixed-width
  `.dat`) to GeoParquet, streaming. `--schema` prints the parsed layout;
  `--selftest` needs no delivery files. It refuses rather than guesses when the
  `.dfn` record width does not match the `.dat`, or when no `EPSGCODE` is
  given.
- `assemble_example.py` builds `example/` from the source archives. It is
  provenance, not something a reader runs.

## Why one control file

The parameters now in `config/udf_example.toml` were spread across four places:
constants in `udf_realdata/udf_system.jl`, constants in
`ek_prerun/prerun_config.jl`, argv in `run_seven.sh`, and literals inside the
stage scripts. Two runs could differ in a way no log recorded. The control file
is the unit of provenance: a result is quotable only alongside the file that
produced it. Unknown keys are refused, not ignored.

**`stage_prelude.jl`** is what every stage includes first. It loads the control
file before any support module (`udf_system.jl` and `xcite_system.jl` parse
their `.des` at load time from a path in it), and resolves the output and log
directories so a stage writes where the control file says. `--site NAME`
restricts a per-site stage to one site.

## Porting history

The stages came from `udf_realdata/` and `ek_prerun/`. What had to change, and
is worth knowing if more is ported:

- Absolute paths (the SNMR archive root, an April-2024 metadata spreadsheet,
  the SkyTEM `.des`) became `[paths]` entries.
- `prerun_data.jl` set four East Kimberley paths as `const`s at load time and
  used them lazily, so a clone carried dead paths that no error mentioned. They
  now resolve from an optional `[paths] ek_root` and raise a message naming the
  cause if anything reaches for them.
- `rockphys.jl` included a file outside the repository (now vendored as
  `common/transform.jl`), and `MAT` was an undeclared dependency. A syntax check
  passed both; loading the files found them.
- `export_udf_soundings.py` needs the 4.9 GB delivery and stays a provenance
  script.

## Relationship to the synthetic workflow

Shared: the install, `Project.toml`, `Manifest.toml`, the depth grid, the
sigma-hat noise convention and the GP prior. That sharing is the point: East
Kimberley, Upper Darling and the synthetics sit on the same grid with the same
noise model, and can be compared without an apples-to-oranges caveat.

Not shared: the stage scripts themselves, which read field inputs rather than a
sealed synthetic truth. There is no truth to compare against here, which is why
`reference/` records observed verdicts rather than errors against a known
answer, and why bore control does the work a truth model does in the
synthetics.
