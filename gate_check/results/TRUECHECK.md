# True-model check, release columns 3P and 2P — 17 Nov 2026

Script `gate_check/truecheck.jl`. Parameter bounds SEDIMENTARY
(0.001–2 S/m per cell). Absolute rule, threshold N + 3.84 = 45.84 on N = 42 data.
Truth hashes `d9de86c430ab` (3P) and `33bf02a46442` (2P). Nothing here feeds back
into the band: the band stays what `profile_minchi2.jl` computed from starts that
never see the truth.

## The answer

| | 3P | 2P |
|---|---|---|
| true window conductance C_true | 0.2063 S | 0.6856 S |
| **min chi2 at C_true**, seeded from the generating model | **26.205** | **27.964** |
| unconstrained best available (`chiref` used) | 25.283 | 25.858 |
| **cost of forcing the truth** | **+0.92** | **+2.11** |
| absolute threshold | 45.84 | 45.84 |
| verdict | **accepted** | **accepted** |

Both solves met the constraint exactly (achieved conductance 0.00% off, constraint
residual ~1e-16) with **zero cells on a bound**, so neither is a clamp artefact.

The cost of forcing the truth's window conductance is 0.92 and 2.11 in chi2,
against the 3.84 that the relative rule allows and the ~20 that the absolute rule
allows. The truth is not merely inside the band — in misfit terms it is close to
the unconstrained optimum. On 2P the solve seeded from the generating model landed
0.50 **below** the unconstrained multi-start reference of 28.459, so that reference
is start-limited; the comparison above uses the lower `chiref` the profile revised
it to, which is the honest one.

## The generating model's own misfit says nothing about the column

| | 3P | 2P |
|---|---|---|
| chi2 of the generating model | 35.449 | 35.449 |

Identical to three decimals, and that is not a coincidence. `makenoisydata!`
reseeds from its own `rseed`, so both columns carry one noise realisation, and the
noise enters as sigma times a standard normal draw. The generating model's chi2 is
then the sum of 42 squared draws from that one realisation — a constant of the
realisation, not a property of the column. Expected 42, observed 35.449.

So "does the truth fit its own data?" is a question about the noise draw, and it
is already answered by Stage 0's `chi2/N = 0.844` check. It is **not** a test of
the column, and it is not evidence about the gate. The informative number is the
one above: min chi2 at forced C_true.

## How the truth is put on the grid matters more than the 3.84 allowance

The generating model is `piecewise_interp` of the 0.25 m truth onto the 65 grid
nodes — the mapping Stage 0 itself used. Two other defensible mappings of the same
truth give:

| mapping of the same truth | 3P chi2 | 2P chi2 | 3P window C | 2P window C |
|---|---|---|---|---|
| piecewise_interp (generating) | 35.449 | 35.449 | 0.2058 S | 0.7435 S |
| arithmetic average in conductivity | 146.088 | 37.828 | 0.2060 S | 0.7680 S |
| geometric average in conductivity | 104.586 | 35.417 | 0.2040 S | 0.7218 S |

On 3P a thickness-weighted arithmetic average of the **same truth** is rejected at
chi2 146 while the node-interpolated truth is accepted at 35.4 — a spread of 111,
nearly thirty times the 3.84 allowance, with the window conductance agreeing to
0.1% throughout. On 2P the three agree to within 2.4.

Two things follow.

1. **"chi2 of the discretised truth" is not a well-defined reference** and must
   never be quoted as the truth's misfit without naming the mapping. This was hit
   here: the arithmetic average was computed first and read as a 3P discretisation
   penalty, which it is not.
2. It is a reason for the gate's design rather than against it. The gate forces a
   window conductance and minimises chi2 over the whole profile. That quantity is
   invariant to how anyone chooses to discretise a reference model, which the
   truth's own chi2 is not.

## A floor on how well 2P's window conductance can be recovered

The grid's own representation of the truth puts 2P's window conductance at
0.7435 S against a true 0.6856 S — **8.4% high**, because the window 17.5–33.5 m
straddles the 18.0 m boundary between the conductive Shepparton clay cap and the
Calivil aquifer, and node interpolation carries the cap's conductivity into the
window's first cell. On 3P, whose window lies inside one unit, the same mapping is
0.3% low.

So on 2P no method working on this grid can be expected to return C_true exactly;
8.4% high is the grid's floor, and it is in the direction of making the aquifer
look more conductive. It is small next to 2P's known ~3x window-conductance
over-estimate, and it does not explain it.

## The trap this check hit, recorded so it is not hit again

`occam`'s `fixval` under `ctype = :cond` is the **conductance in S**
(`src/occam_core.jl:172`), not the window conductivity. The first version of this
script passed `C_true/thick`, forcing one sixteenth of the intended conductance.
On 3P that fell below the floor-limited minimum of 0.0160 S, so the solve pinned
all nine window cells to the 0.001 S/m floor, stopped at 0.0160 S, and returned
`chi2 = 25.998` with `settled = true` — which read as a clean acceptance at
C_true. Nothing in the return value said the constraint had been missed by a
factor of thirteen.

`truecheck.jl` now recomputes the achieved conductance from the returned model and
discards any solve more than 1% off, rather than trusting `settled`.

`gate_check/profile_minchi2.jl` was checked and does not have this bug: it passes
`10^lt[i] * thick`, and it accepts a solve into the band only `if r.conres < 1e-8`.
The published bands are unaffected. That guard is also why the end nodes of the 2P
profile read `Inf`: an `Inf` there means **no solve met the constraint**, which is
not the same as a rejection.
