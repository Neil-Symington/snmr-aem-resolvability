# Part A — diagnosing the existing runs

Branch `sampler-convergence-test`, from `749d730` (which is `e214185` plus the
five production-run defect fixes). 13 Nov 2026.

Interim. A1–A5 need chain files that are not on this machine; the measurements
that do not need them are below, together with what is now running to supply
the rest.

## A0. The blocking fact: the production chains are not here

Part A assumes the 2P production run's output is available to re-analyse. It is
not. The newest chain binary of non-zero length anywhere under
`AEM_SMR_Resolvability_Quickrun/` is from **10 Nov**, and belongs to the retired
field workflow (`fieldrun/outputs/full_validate/`, sites SMR_UD21 and SMR_UD24).

What exists for the synthetic workflow:

| path | what | state |
|---|---|---|
| `data/aemFIX_smoke/` | Stage 1, smoke, 2001 samples | intact, 1 ensemble |
| `data/smr_runs/*FIXT_e2,e3,e4*` | Stage 2, smoke, 2001 samples | intact, 3 ensembles |
| `data/smr_runs/*FIXT_e1*` | Stage 2, smoke | **0 bytes**, truncated 11 Nov |
| `data/cond_C2_typicalFIXT.jld2` | Stage 2 pooled output, smoke | intact |
| `data/truth_murray2P.jld2` | test two truth column | intact |
| 2P production chains | — | **absent** |

The reviewer's run was on a different machine (Dell XPS 13 9360) and their
attachments (`convergence.jl`, `convergence_2P.txt`, the figures) did not come
with `BUG_REPORT.md`. So A1, A3, A4 and A5 cannot be answered from stored
output, and A2 never could be — test one had never been run at all.

**Correcting my own earlier note.** `CONVERGENCE_PROBLEM_20261113.md` says every
chain file on this machine is zero bytes. That is wrong and too broad: only the
four `e1` files were truncated. Ensembles e2, e3 and e4 survived, and they are
what the acceptance measurements below are made on.

### What I did instead, and why

The brief says to follow the intent where it conflicts with how the code works.
The intent of Part A is an independent diagnosis. That needs chains at
production settings, which is exactly what B1 is. So Part A and B1 are merged:
B1 is running now and its output answers A1–A5 as well as being the control.

One consequence has to be stated plainly. The trans-D chains are **unseeded**
(`tools/aem_HYP.jl:84-110` documents why: Julia's task-local RNG defeats worker
seeding, and HiQGA issues a fresh `remotecall_fetch` per step). So no run can
reproduce the reviewer's realisation. A1 can only ask whether the *statistics*
reproduce, not the numbers. That is the meaningful question anyway, and it is
the one I will answer.

Incidentally this is also what makes the test possible: four repeats of one
command are four genuinely independent ensembles.

## A1a. The diagnostics are independently implemented and validated

I did not write new diagnostics. `ek_prerun/prerun_lib.jl` already contains
rank-normalised split R-hat and multi-chain ESS after Vehtari et al. (2021),
built for the field pre-run. Reusing it means the synthetic and field workflows
are judged by one implementation.

It reproduces the reviewer's published benchmarks:

| test | this implementation | reviewer | theory |
|---|---|---|---|
| iid normal, 4 x 10,000 | 1.0000 / 39,351 | 1.0000 / 39,313 | 1.00 / 40,000 |
| AR(1) phi = 0.90 | 1.0019 / 2,000 | 1.0017 / 2,068 | 1.00 / 2,105 |
| AR(1) phi = 0.99 | 1.0266 / 175 | 1.0047 / 187 | 1.00 / 201 |
| one chain shifted 0.5 sigma | 1.0257 | 1.029 | > 1.01 |

Two differences, both small and both in the conservative direction:

- ESS runs about 5–13% below theory. The autocorrelation sum truncates on the
  initial **positive** sequence; the reviewer specified Geyer's initial
  **monotone** sequence. Monotone is the slightly less conservative of the two.
  Both cap at lag 2000, so both are upper bounds on the real chains.
- R-hat 1.027 against the reviewer's 1.005 on AR(1) 0.99 is realisation noise:
  with ESS near 180, R-hat itself has visible sampling variability.

If a later number disagrees with the reviewer's, it is the data, not the tool.

## A1b. The gate is unwired here — and wired everywhere else

Confirmed for `release/`: `GATE_RHAT = 1.05`, `GATE_MIN_ESS = 100`,
`GATE_DEPTHS = [5, 15, 30, 60, 100]` m are read by no code. The only use outside
their definition is `src/prerun_config.jl:176`, which interpolates
`GATE_MIN_ESS` into a **documentation string**. The reviewer's critical finding
stands exactly as written.

The part the bug report did not reach: **the same gate is fully wired in the
sibling workflows.** `ek_prerun/run_prerun_mcmc.jl`, `ek_realdata/run_smr_transd.jl`,
`udf_realdata/run_smr_transd.jl` and `ssc_prerun/run_ssc_mcmc.jl` all compute
R-hat and ESS at `GATE_DEPTHS` and branch on `GATE_RHAT`/`GATE_MIN_ESS`.
`ssc_prerun` goes further and carries a second threshold, `GATE_RHAT_PROD = 1.10`,
with a comment saying 1.05 is the literal spec gate and 1.10 is "the one to
quote".

That changes the character of the finding. This is not a check nobody got round
to writing. It is a check that exists, runs, and was dropped on the way into the
synthetic pipeline. It also means Part C needs no new science, only wiring — and
that `ssc_prerun`'s two-threshold design is a precedent for the open question of
what a failure should do.

## A3. Acceptance rates — method fixed, and a measurement of my own retracted

The `misfits_*_s.bin` files are **text**, 11 columns. From HiQGA's own format
string (`TransD_GP_MCMC.jl:1159`):

```
chain | iter | birth | death | position | property | dc | nodes | misfit | T | dcvalue
```

The critical subtlety: columns 3–7 are **windowed**, not cumulative.
`get_acceptance_stats!` recomputes them every `stat_window = 100` samples and
then zeroes the counters, and rows are written every `save_freq`. So the logged
value is a 100-iteration rate, repeated across rows until the next update.

**This invalidates a second measurement of mine.** Earlier today I read the last
row of each chain and reported it as that chain's acceptance rate, getting 0.00%
birth and death at T = 1. That is one 100-iteration window — roughly three birth
attempts — not a rate. It was noise. Averaging over window boundaries
(`(iter-1) % 100 == 0`) gives the real picture, and it is much less dramatic:

Stage 2 (SNMR), smoke, intact ensembles, mean over stat windows, %:

| ensemble | scope | birth | death | position | property | dc |
|---|---|---|---|---|---|---|
| e2 | T = 1 | 6.6 | 6.3 | 32.2 | 28.2 | 17.7 |
| e2 | hottest | 11.0 | 14.2 | 32.7 | 31.5 | 22.6 |
| e3 | T = 1 | 11.9 | 11.9 | 35.0 | 32.5 | 27.3 |
| e3 | hottest | 54.0 | 48.5 | 57.0 | 73.9 | 89.8 |
| e4 | T = 1 | 3.2 | 2.5 | 26.3 | 25.5 | 9.7 |
| e4 | hottest | 29.0 | 29.5 | 42.1 | 54.2 | 54.4 |

Stage 1 (AEM), smoke, test one, same method, %:

| ensemble | scope | birth | death | position | property | dc |
|---|---|---|---|---|---|---|
| e1 | T = 1 | 35.2 | 33.7 | 41.5 | 58.9 | 39.4 |
| e2 | T = 1 | 42.0 | 36.0 | 44.3 | 58.5 | 38.7 |

Birth and death at T = 1 are **not** zero. But the two stages are very
different: Stage 2 accepts a dimension change at T = 1 in 2.5–12% of attempts,
Stage 1 in 34–42%. In Stage 2 the hottest chain accepts 3 to 10 times more
readily than T = 1; in Stage 1 there is no such gap.

Smoke runs, 2001 samples, about 20 windows at T = 1, and the two stages are on
different columns, so these are indicative and nothing more. They are reported
because they are the first honest look at the quantity the convergence note
named as the place to look next.

## A4. Dimension mixing — the exact measure

Acceptance rates are windowed and awkward. The node count is logged every row,
so the fraction of consecutive post-burn-in T = 1 samples at which **k actually
changes** is exact and needs no bookkeeping assumptions. That is the headline
measure of dimension mixing and it is implemented in
`convergence_test/diagnose_stage1.jl`.

On the Stage 1 smoke run it is 86% and 90%, with k spanning 13–26 and 10–36 and
ensemble medians 20.5 and 21.5. Stage 1's dimension moves freely at smoke
length. Whether that survives at 200,001 samples, and what Stage 2 does, is what
B1 will say.

## A5. Decision quantities — defined, implemented, pending data

Implemented in `diagnose_stage1.jl` (Q2, Q3) and `diagnose_stage2.jl` (Q1), both
reporting median, 5th and 95th percentile per ensemble plus R-hat, bulk ESS,
tail ESS and spreads as factors.

Windows, with the one judgement called out:

| column | Q1/Q2 window | Q3 package | note |
|---|---|---|---|
| test two, 2P confined | 17.5–33.5 m | 4.0–34.0 m | window is AQUIFER 1 (18–34 m); package is clay cap + aquifer, as the brief specifies |
| test one, 1P unconfined | 7.5–23.5 m | 6.0–40.0 m | aquifer is 6–40 m and **unconfined**, so there is no clay cap; the "equivalent" package is the aquifer itself |

Test one's truth column did not exist and has been built from
`synthetics/config/unconfined_aquifer.toml` as `data/truth_murray1P.jld2`.

Stage 2 needs no re-run to be diagnosed once a cond file exists: Stage 2 already
stores the four ensembles **separately** in `ens`, a vector of four matrices. It
is only the pooling downstream that hides the disagreement.

Run on the one pooled file that does exist (smoke, FIXT, 102 draws per
ensemble), the signature is already the one the reviewer describes:

```
ensemble      p05       p50       p95
1          0.2261    0.2295    0.2314
2          0.2051    0.2193    0.2236
3          0.2011    0.2065    0.2095
4          0.2187    0.2207    0.2228

Q1 window water content   R-hat 2.559   ESS bulk 51   median spread 1.11x
```

Each ensemble is internally tight — 5th to 95th percentile spans about 0.005
v/v — while the ensembles sit 0.023 v/v apart. The disagreement between
ensembles is roughly four times the width each one reports. That is what drives
R-hat to 2.56 despite the medians differing by only a factor of 1.11, and it is
the reason pooling looks like a wide posterior rather than a failure. Smoke
settings, so the number is not the production number; the shape is the point.

## Ladder arithmetic for B3

HiQGA builds the ladder as `T = 10 .^ range(0, log10(Tmax), length = nchains -
nchainsatone + 1)` (`MCMC_Driver.jl:35-36`), confirmed against the temperature
column of a real run: 1.000, 1.357209, 1.842016, 2.500000.

| chains | temperatures | largest gap in 1/T |
|---|---|---|
| 4 (current) | 1, 1.36, 1.84, 2.5 | 0.263 |
| 7 (upstream report) | 1, 1.17, 1.36, 1.58, 1.84, 2.15, 2.5 | 0.142 |

Log spacing in T puts the **coarsest** rung between T = 1 and its only
neighbour, which is the single link by which the posterior chain receives
anything. Seven chains roughly halves it.

Note for B3: the spacing is derived from `Tmax` and `nchains`, so "tune the
spacing" means choosing those two, unless HiQGA accepts an explicit ladder —
to be checked before B3 runs.

## What is running

`convergence_test/run_b1_all.sh`: B1 at production settings, 200,001 samples,
4 chains, four independent ensembles, test two then test one. Two ensembles at a
time (each takes 4 workers; this machine has 8 threads). About 5 h per ensemble,
so about 10 h per column.

Output goes to `data/aem<TAG>_b1e<N>/`, never `data/aem<TAG>/`, so no production
posterior can be overwritten. `numbers.csv` is untouched.

Sampler knobs are now environment-overridable in `tools/aem_HYP.jl`
(`SCT_NMAX`, `SCT_TMAX`, `SCT_NCHAINSATONE`), defaulting to the production
values, so B2 and B3 need no further code change and the control is exactly the
current sampler.

## Two questions from the brief, now settled in the code

**B3: the ladder cannot be given explicitly.** `Chain()` in `MCMC_Driver.jl:25-48`
builds it from `Tmax`, `nchains` and `nchainsatone` and takes no ladder
argument. So B3 is specified as a (`nchains`, `Tmax`) pair, and "tune the
spacing so swap rates are 20-40%" becomes "choose the pair that gets there".
Seven chains at Tmax 2.5 is the upstream report's configuration and halves the
bottom rung's gap, which is where it matters.

Worth recording for `nchainsatone`, since the arithmetic is easy to get wrong:
the ladder length is `nchains - nchainsatone + 1` and the extra chains are
appended at T = 1. With 4 chains and `nchainsatone = 2` the ladder becomes
1, 1.58, 2.5 **plus** a second chain at T = 1 — a second posterior chain bought
by giving up a rung. It is the only way to get a within-run R-hat, but it
coarsens the ladder unless `nchains` rises too.

**B4 is skipped: it needs new code.** There is no prior on k to tune. `nmin` and
`nmax` are the whole of it, enforced by refusing a birth at `nmax` and a death
at `nmin` (`TransD_GP_MCMC.jl:897`, `:850`); no prior-ratio term for k appears
in the acceptance. Favouring fewer layers would mean adding one to HiQGA's
birth/death acceptance, which the brief rules out.

B2 covers part of B4's intent anyway. Taking `nmax` from 40 to 10 changes the
prior from uniform on {2,40} to uniform on {2,10}, which *is* a prior favouring
fewer layers — just the blunt version of it.

## Open item

A swap-rate measurement, which B3 is judged on. Swaps are not logged. The
temperature column recovers them by tracking which chain holds each rung, but
only at `save_freq` resolution, so it undercounts when two swaps fall inside one
interval. Worth doing for the T = 1 rung, which is the one that decides whether
the posterior chain hears anything; the undercount makes it a lower bound, which
is the safe direction.
