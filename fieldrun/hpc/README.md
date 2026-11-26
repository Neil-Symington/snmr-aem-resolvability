# fieldrun/hpc — the field workflow as Gadi task farms

`run_batch.sh` runs a whole field stack on one machine. The scripts here run the same stages as independent tasks, so a run of any size spreads over as many Gadi nodes as it needs. Every task can be resumed.

## The phases

| phase | one task per | what a task runs | size on Gadi `normal` |
|---|---|---|---|
| `aem` | site | `run_aem_occam.jl`, then `run_aem_transd.jl` | nchains + 1 processes, one core each, about 6 GB |
| `snmr` | site × kernel draw × replicate ladder | `run_smr_transd.jl --ens e --rep r` | nchains_snmr + 1 processes, two per core, about 7 GB |
| `post` | site | `--combine`, `detect_aquifer.jl`, `run_aem_profile.jl`, `gate_water_quality.jl`, `report_smr.jl`, `report_html.jl` (each with `--site`) | one process, about 4 GB |
| `index` | run | `report_html.jl --index` | one core |

Dependencies are per site. A site's SNMR ladders wait only for the shard that ran that site's AEM, and its post task waits only for the shards holding its ladders. Gadi has no PBS job arrays, so `submit.sh` cuts each phase into shards of one node each and chains them with `-W depend=afterany`.

## Use

```
bash fieldrun/hpc/submit.sh <config> --env <env.sh> --project <code> --dry-run   # look first
bash fieldrun/hpc/submit.sh <config> --env <env.sh> --project <code>             # submit
bash fieldrun/hpc/submit.sh <config> --local --slots-snmr 2                      # no PBS (a Mac)
```

- **`env.sh`** must put Julia 1.12 on the PATH and set `JULIA_DEPOT_PATH` and `REPO`.
- **Phases:** by default they follow `[run] stages`. Choose them with `--phases aem,snmr,post,index`.
- **Sizing:** override with `--slots-aem/-snmr/-post` (tasks per node) and `--t-aem/-snmr/-post` (minutes per task). `--wmax` caps a shard's walltime in hours (default 10).

### Resuming

Run the same command again and it submits only what isn't finished.

- **Bookkeeping** lives in `<logs>/batch/`:
  - `tasks/<id>.log` holds one log per task;
  - `done/<id>` marks a finished task, and `done/<id>.<step>` a finished step within it;
  - `failed/<id>` marks a failed task.
- **Rerunning a task:** delete its done marker. Also delete any output that its stage skips when the output already exists (the profile does this).

### Post: required and optional steps

- **Required:** combine and detect. If either fails, the post task fails.
- **Optional:** profile, gate and the matplotlib report. A site with no detected interval has no profile and no gate, so their failure is recorded in `failed/<id>.soft` and the viewer page is still written.

## The viewer

`report_html.jl` writes `<outputs>/viewer/<site>.html`, plus `index.html` linking every site.

- **The page itself** is Neil's `tools/make_report.jl` page from branch `m1-production-rerun`, vendored in `stages/viewer/`. Its header says what changed for field data.
- **Opening it:** each page is self-contained, so copy the `viewer/` folder back and open `index.html`.

## Sizing notes (BDU, Nov 2026)

- **SNMR packing:**
  - 12 ladders per node: 43 min for 32 ladders, 84 GB, 46% CPU.
  - 6 per node: 70 min.
  - With s integrated out, a ladder takes about 60% of the time, so 16 per node is worth trying (`--slots-snmr 16`).
- **The profile** is about 200 Occam solves per site on one core, roughly 20–25 min, and is most of `post`.
- **Unused walltime isn't charged,** so the per-task estimates err long.
