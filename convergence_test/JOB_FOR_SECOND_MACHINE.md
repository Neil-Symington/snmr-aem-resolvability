# Job for the second machine (the one that ran the 2P production job)

Two jobs. The first takes about twenty minutes and **only that machine can do
it**. The second is long and keeps it busy for a day or two.

Nothing here writes to `data/aem2P/`, `data/cond_*`, `numbers.csv` or any other
production output. New runs go to `data/aem<TAG>_<setting>e<N>/`, which cannot
collide. Job 1 reads only.

## Setup, once

```bash
cd <wherever snmr-aem-resolvability is>
git fetch origin
git checkout sampler-convergence-test
```

If that machine has local commits or edits, stash them first; the branch is
`749d730` plus the convergence-test harness, so nothing of theirs is lost.

---

## Job 1 — read the production chains (about 20 minutes, do this first)

This is the irreplaceable one. The production chain files exist **only** on that
machine, so A1, A3, A4 and A5 cannot be answered anywhere else. My machine has
nothing newer than 10 Nov.

```bash
julia --project=. convergence_test/diagnose_production.jl 2>&1 | tee /tmp/production_2P.txt
bash convergence_test/collect_production.sh
```

The first command auto-discovers `data/aem2P/` and `data/cond_C2_*2P*.jld2`. If
it says NOT FOUND, pass the paths:

```bash
julia --project=. convergence_test/diagnose_production.jl data/aem2P data/cond_C2_typical2PT.jld2
```

**Send back two things:** `/tmp/production_2P.txt` and
`/tmp/production_2P_diagnostics.tar.gz`. The tarball holds the misfit files
(plain text: chain, iteration, five acceptance rates, node count, misfit,
temperature) and the Stage 2 cond file. It deliberately leaves out the large
`points_`/`models_` binaries, so expect tens of megabytes, not gigabytes.

That gives, on the real production run:

- acceptance rates per move type, at T = 1 and in the hottest chain (A3)
- node count against iteration, and how often k actually changes (A4)
- a lower bound on the swap rate, from how often the T = 1 slot changes hands
- Q1, Q2, Q3 with median, 5th and 95th percentile (A5)
- an independent recomputation of the reviewer's R-hat and ESS (A1)

Stage 1 ran one ensemble, so its R-hat is a **segment** R-hat across four
consecutive pieces of a single chain. It detects drift and cannot detect a mode
the chain never visited. The script prints that caveat with the number every
time, because it looks like an R-hat and is not one. Stage 2 genuinely has four
ensembles, so its R-hat is the real thing.

---

## Job 2 — run B2, then B3 (about a day and a half)

My machine is doing B1, the control: production settings, four independent
ensembles, test two then test one, about twenty hours. Do **not** duplicate it.

Build the test-one column first; it has never existed on either machine:

```bash
julia --project=. synthetics/scripts/build_truth.jl synthetics/config/unconfined_aquifer.toml
cp synthetics/output/unconfined_aquifer_truth.jld2 data/truth_murray1P.jld2
```

### B2 first — maximum layers 10, nothing else changed

B2 is the cheapest setting and the most likely to work, and the brief says to
stop as soon as something passes. So it goes first, not last.

```bash
SCT_NMAX=10 bash convergence_test/run_stage1_setting.sh b2 2P 200001 4
SCT_NMAX=10 bash convergence_test/run_stage1_setting.sh b2 1P 200001 4
```

Each line is four independent ensembles, two at a time, roughly ten hours per
column and probably less than B1 because the models are smaller. Then:

```bash
julia --project=. convergence_test/diagnose_stage1.jl 2P b2
julia --project=. convergence_test/diagnose_stage1.jl 1P b2
```

**If both columns pass on Q2 and Q3, stop and send the output.** That is the
answer and nothing further needs running.

### B3 only if B2 fails — maximum layers 10 plus the seven-chain ladder

```bash
SCT_NMAX=10 bash convergence_test/run_stage1_setting.sh b3 2P 200001 7
julia --project=. convergence_test/diagnose_stage1.jl 2P b3
```

Seven chains at Tmax 2.5 is the upstream report's configuration. The ladder
becomes 1, 1.17, 1.36, 1.58, 1.84, 2.15, 2.5, which halves the gap in 1/T
between T = 1 and its only neighbour — from 0.263 to 0.142. That gap is the
single link by which the posterior chain hears anything, and log spacing in
temperature puts the coarsest rung exactly there.

Seven chains need seven workers, so the driver drops to one ensemble at a time
and this takes about twice as long per column as B2. Run test two only; test one
can wait on the result.

### What not to bother with

**B4 is skipped and needs no attempt.** There is no prior on k to tune in HiQGA:
`nmin` and `nmax` are the whole of it, enforced by refusing a birth at `nmax`
and a death at `nmin`. A prior favouring fewer layers would mean new code in
HiQGA's birth/death acceptance, which the brief rules out. B2 covers the intent
bluntly anyway, since `nmax` 40 to 10 changes the prior from uniform on {2,40}
to uniform on {2,10}.

**Do not re-run Stage 3 onward, the gate, or the saline columns.**

---

## If something goes wrong

The sampler knobs are environment variables read by `tools/aem_HYP.jl`:
`SCT_NMAX`, `SCT_TMAX`, `SCT_NCHAINSATONE`. Unset, they are the production
values, so an un-prefixed command reproduces the current sampler exactly.

Per-ensemble logs land in `convergence_test/logs/`. Each ends with a line
`ENSEMBLE_DONE e=<n> rc=<code> minutes=<m>`; a non-zero `rc` is a failure worth
sending rather than retrying.

Julia must be 1.12.x. 1.13 fails to resolve HiQGA.
