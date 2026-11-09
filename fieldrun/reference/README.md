# Reference values for the batch field run

A run that only prints a table can only be admired. These files let a run be
**checked**: stored values, a script that compares a finished run against them,
and a non-zero exit code when it disagrees.

```
julia --project=. fieldrun/reference/check_run.jl --config fieldrun/config/udf_example.toml
```

Run it from the repository root, after `run_batch.sh` has finished. It prints
one line per check and exits 1 on any disagreement.

## Why this is here

The first complete two-site run of this workflow ended with:

```
      SMR_UD21 snmr_transd ok
      SMR_UD24 snmr_transd ok
######## detect
      detect ok
######## udf_example complete
```

Every one of those lines was true in the only sense the runner meant: the
process exited 0. In fact both SNMR trans-D sites had thrown, the stage caught
the exceptions and returned normally, and all three sweep stages then printed
an empty table over "no trans-D file for SMR_UD21, skipping". The run produced
nothing and said it was complete.

Three things came out of that, and all three are worth knowing about if you
extend these scripts:

- a stage that catches a per-site exception must still exit non-zero,
- a sweep that produced no rows is a failed sweep, not an empty table,
- and `run_batch.sh` now exits non-zero and prints `INCOMPLETE` if any stage
  left a failure marker, because "complete" previously meant only that the
  loop had ended.

## What is checked, and what is not

| Section | What it holds | Tolerance |
|---|---|---|
| `[deterministic]` | Occam χ²/N, data count, σ(z) window geometric means | `tol`, relative (default 1e-4) |
| structural (in the script, not the file) | files present, bands ordered and positive, interval inside the grid, verdict vocabulary | exact |
| `[full]` | converged posterior values | not recorded yet |

**The Occam inversion has no Monte Carlo in it.** χ²/N, `ndata` and the window
geometric means of σ(z) reproduce to the printed digits on any machine, so they
are compared tightly and a disagreement there is a defect to go and find.

**Everything downstream of trans-D is a posterior.** In `mode = "smoke"` that
posterior is 2001 samples and is *not converged* — the stage says so itself,
printing `gates FAIL` for R̂ and ESS. Comparing its quantiles against a stored
number would be checking noise, so the script does not. It checks instead what
must hold of any valid run: that the files exist, that the conductivity band is
positive and ordered, that the target interval lies inside the inversion grid,
and that each verdict is a sentence the stages can actually produce.

So **a passing smoke check means the workflow is wired correctly end to end.**
It does not mean the numbers are quotable. Smoke values are deliberately *not*
promoted into `[full]`; that section stays absent until a full-mode run has
been done, and the script says so rather than inventing them.

The σ(z) comparison uses the **geometric** mean over each window, `10^mean(log10 σ)`,
which is what `unit_geomean` in `common/aem_lib.jl` reports. The arithmetic
mean over the same cells differs by roughly 1.6× in this ground, so using the
wrong one looks exactly like a real disagreement.

## Updating the reference

Only when a number *should* have changed — a corrected forward model, a
different noise specification, a fixed defect. Re-record it, say in the commit
message which input changed and why, and keep the mode it was recorded in
(`mode` at the top of the file) honest: the script warns when a run's mode and
the reference's mode differ.

Do not widen `tol` to make a failing check pass. The Occam numbers are
deterministic; if they move, something moved them.

## What a smoke run is worth, measured

Smoke mode is not a scaled-down answer, it is a different answer. Both of these
are real runs of the same two sites through the same code:

| | smoke, 2001×2 | full, 150001×4 (the production run) |
|---|---|---|
| noise scale `s` p50, UD21 | **1268** × σ̂ | **7.4** × σ̂ |
| UD21 constraint depth | 20.8 m | 42.7 m |
| UD21 detect verdict | `NOT RESOLVED (P < 0.9 for 1 m) above 21 m` | `AQUIFER (W >= 2 m, P >= 0.9) above 43 m` |
| UD24 detect verdict | `AQUIFER … above 39 m` | `AQUIFER … above 41 m` |
| gate verdict, both sites | `NOT CERTIFIED` | `NOT CERTIFIED` |

**The UD21 verdict flips.** `s` is a free nuisance that absorbs misfit while the
model is still poor, and `burninfrac = 0.5` on a 2001-sample chain discards only
1000 samples — nowhere near enough for it to come back down. At 1268× σ̂ the
inferred noise swamps the signal, the constraint depth halves, and a resolved
aquifer becomes an unresolved one. This is why no posterior quantile is stored
in this directory, and why a smoke run must never be quoted.

## The full-mode target

These are the production run's values, read from `udf_realdata/outputs/`. They
are recorded here as the target a full-mode run should land near, **not** as
automatic checks — see the caveat below.

| site | detect | profile band, S/m | gate |
|---|---|---|---|
| SMR_UD21 | `AQUIFER … above 43 m`, window 18.6–30.0 m (modal best-10m screen), zlim 42.7 | 0.05727–1.59200, χ²/N 0.3507 | X_cert p95 53.30 → `NOT CERTIFIED` |
| SMR_UD24 | `AQUIFER … above 41 m`, window 14.4–31.8 m (consensus interval), zlim 40.6 | 0.03749–1.05269, χ²/N 0.3100 | X_cert p95 27.61 → `NOT CERTIFIED` |

**The caveat, and it is a large one.** That production run does not pass its own
convergence gates:

```
SMR_UD21   rhat_max 1.89   ess_bulk_min 36.9   gates_pass false   stationary false
SMR_UD24   rhat_max 1.94   ess_bulk_min 89.5   gates_pass false   stationary false
```

`[convergence]` in the control file says a run failing these "is not quotable",
and by that standard these numbers are not either. Two independent full runs
would not be expected to agree closely on the posterior quantiles, which is why
they are documented rather than asserted. Treat the table as the shape a full
run should have, and treat a large departure from it as worth investigating in
*both* runs.

## Port faithfulness, as far as it has been established

The Occam inversion is bit-identical to the pre-port code: `max |Δm| = 0.000e+00`
over all 70 cells and χ²/N agreeing to twelve figures, for both sites. That is
what `[deterministic]` asserts, so the section doubles as a port check.

The SNMR half has **not** been shown faithful. A full-mode port run compared
against the table above is the test, and it has not been done.
