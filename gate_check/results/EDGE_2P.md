# 2P coarse upper edge, SEDIMENTARY — 17 Nov 2026

Driver `gate_check/profile_minchi2.jl`, `--nodes -0.4,-0.35,-0.3 --mbounds
-3.0,0.30103 --lam 0.01,0.1,1,10 --starts all`, window 17.5–33.5 m (support
16.00 m, 8 of 65 cells), N = 42, truth hash `33bf02a46442`. 44 solves, 220 min.
Output `profile_2P_upper_sed.txt`, store `profile_2P_upper_sed.jld2`.

## The answer

**The absolute upper edge lies between 6.3697 S and 7.1469 S — bracketed to
0.0500 dex**, better than the 0.1 dex asked for.

| C (S) | min chi2 | nbound | winning start | absolute rule (45.84) |
|---|---|---|---|---|
| 6.3697 | 37.888 (36.577 pooled with the 9-node run) | 15 | `cold_m0@l2=0.01` | **accepted** |
| 7.1469 | 53.923 | 19 | `cold_m0@l2=0.01` | **rejected** |
| 8.0190 | `Inf` — no solve met the constraint | — | none | unresolved |

Taken with the 9-node grid, 2P's absolute band is **0.0506–~6.7 S, a factor of 126**. 2P's window conductance is effectively unconstrained by the data over
two orders of magnitude, which is the same fact as the known flatness of its
likelihood (chi2 varies about 0.5 across 0.005–3 S).

**Relative rule, from the existing 9-node grid, no new solves:** band
0.0506–5.0596 S, upper edge bracketed by 5.0596 S accepted and 6.3697 S rejected
against `chiref + 3.84 = 29.70` — a width of 0.0999 dex, which meets the 0.1 dex
target but does not beat it; no new solves were spent to tighten it, the lower
edge being the only thing below and that feeding a salinity floor that is zero.

## Why 8.0190 S is still `Inf`, and why it no longer matters

`Inf` at a node means **no solve met the constraint to `conres < 1e-8`**, not that
the data rejected it. The node is not infeasible: 0.501 S/m against a per-cell
ceiling of 2 S/m leaves room, unlike the bottom node at 0.0160 S which sits
exactly on the 0.001 S/m floor. The 9-node run had tried it with warm sweeps only
at λ² ∈ {0.01, 0.1}; this run added cold starts from `mref` and from `m0` and
λ² ∈ {1, 10}, and it still could not be solved onto the constraint. So the earlier
`Inf` was not a warm-start artefact.

It no longer matters, because 7.1469 S **below it solves cleanly and is rejected**
at chi2 53.923. The edge is bracketed without it, and no further solves were spent
on it.

## Two lines in the output that must not be quoted

The run prints `truth window C = 0.6856 S … absolute: **OUTSIDE**` and
`relative band: NaN – NaN`. Both are artefacts of a three-node grid that begins at
6.3697 S: the truth at 0.6856 S lies **below** the grid, and all three nodes sit
above the relative threshold. Neither is a statement about 2P. The truth's standing
is in `TRUECHECK.md` — min chi2 at C_true = 27.964, accepted — and the relative
band is the 9-node one above. The same artefact appeared in the earlier
`profile_2P_upper` run, whose grid ran 5.06–160 S.

## This edge is another instance of the `nbound` blocker

The accepted node at 6.3697 S has **15 of 65 cells on a bound**, and the rejected
one at 7.1469 S has 19. Under the release path's `good[i] = nf.nbound == 0 && …`
both are discarded as "unreachable", and 2P's upper edge falls back to 5.0596 S,
the highest node with nothing on a bound. So the criterion moves 2P's upper edge
by 0.1 dex as well as moving 3P's band off the truth. See the blocker section of
`DECISIONS.md`.

## One inconsistency worth recording

At 6.3697 S this run's lowest solve is 37.888, while the 9-node run reached
**36.577** at the same node with a *smaller* λ² set. Enriching the search cannot
raise a minimum, so the difference is the warm chain: the 9-node sweep arrived at
this node from neighbours at 5.0596 and 8.0190 S, this one only from 7.1469 S, and
the two warm starts land in different local minima. The band uses the lowest solve
ever obtained at a node, so the pooled value 36.577 is the one to quote. It does
not change any verdict — both are far below 45.84 — but it is a reminder that a
single profile run returns an upper bound on min chi2, not min chi2.
