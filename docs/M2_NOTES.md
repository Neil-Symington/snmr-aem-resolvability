# Notes to carry into M2

Things found during M1 that belong in the pipeline, not in a one-off script.

## Julia soft scope: all loops in functions, and test that a warm start propagates

At top level, `warm = r.m` inside a nested loop binds a **new local** on each
iteration rather than updating the outer binding. It does not error. The loop
runs, the profile is produced, and every solve silently starts cold.

This bit twice in M1: once in `gate_check/profile_minchi2.jl`, where it would
have produced a profile with no warm starts at all and reported it as if it had
them, and once in `convergence_test/diagnose_stage1.jl` (`zall`, `zbnd`), where
Julia at least warned.

**Rule for M2:** every loop that carries state across iterations lives in a
function. No `global` patches in the release path.

**Test for M2:** assert that a warm start actually propagates — run two nodes
where the second's answer depends on the first's model, and check that the
second solve's iteration count (or its starting misfit) is what a warm start
would give, not what a cold one would. A test that only checks the final chi2
will pass with the chain broken, because a cold start reaches a similar answer;
it is the *cost* and the *path* that differ.

## Solve cost is dominated by start quality

Measured in M1.1: a warm solve from the adjacent node is about 45 s; a cold
solve from `mref` or `m0` is about 7 minutes — roughly nine times. Any design
that cold-starts every node is about an order of magnitude more expensive than
the warm-chained sweep the gate currently does. Budget accordingly, and do not
estimate a constrained-solve cost from a warm-chained one.

## Parameter bounds must not be a mutable global

`src/occam_core.jl` holds the clamp on log10 sigma in `Ref`s, set by
`setbounds!` and read by `getbounds()`. Two runs in one session can silently use
different bounds. **M2: thread the bounds from config as an explicit argument,
record them in the manifest, and remove the mutable global from the release
path.**

## Exit status

`run_site.jl` returned 1 on the Occam 2S column after writing every output
correctly, the same trailing-command pattern that made `run_production.sh` exit
1 after a successful run. In M4 a non-zero exit on a good column would read as a
failure. Fix the pattern wherever a script's last statement is a conditional.

## One place for the operator builder

`operator_from_stage0` is now copied verbatim in three places
(`stages/stage4_gate_v2.jl`, `gate_check/profile_minchi2.jl`,
`gate_check/cold_refit.jl`) because the stage file ends in `main()` and cannot
be included. M2 should move it, the window/overlap construction and the
constrained-solve wrapper into one module both the stage and any diagnostic can
include.
