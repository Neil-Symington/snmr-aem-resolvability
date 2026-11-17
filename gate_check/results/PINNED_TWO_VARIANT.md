# Pinned cells by depth, both bound variants — 17 Nov 2026

Driver `gate_check/pinned_two_variant.jl`, from stored models only, no new
solves. Full tables in `pinned_two_variant_3P.txt` and `pinned_two_variant_2P.txt`.

Stores used: 3P WIDE `edges_3P.jld2`, 3P SEDIMENTARY `edges_3P_sed2.jld2` (the
converged gate), 2P WIDE `profile_2P_models.jld2`, 2P SEDIMENTARY
`profile_2P_sed.jld2`. Absolute rule throughout, threshold 45.84 on 42 data.

## Why by depth

A count of cells on a bound over all 65 cells mixes two different things. A
shallow cell on a bound says something about the data: the fit wanted a
conductivity outside the stated limits where the data had information. A cell at
250 m on a bound says nothing — the data cannot see it, the regularisation puts it
where the gradient points, and it lands on a limit. Only the shallow count bears
on whether a node's acceptance owes anything to where the clamps were set. Deep
counts are reported separately, never folded in, so the two can be recombined.

## The answer, at the defensible split (frac = 0.02)

| column | variant | bounds (S/m) | absolute band | largest clamp-free C | truth |
|---|---|---|---|---|---|
| 3P | WIDE | 1e-4 – 10 | …–0.4509 S | **0.3192 S** | 0.2063 S, below it |
| 3P | SEDIMENTARY | 1e-3 – 2 | 0.0506–0.4320 S | **0.1795 S** | 0.2063 S, **above it** |
| 2P | WIDE | 1e-4 – 10 | 0.0051–6.3697 S | **6.3697 S** (all of it; 5.0596 S at frac 0.01) | 0.6856 S, below it |
| 2P | SEDIMENTARY | 1e-3 – 2 | 0.0506–~6.7 S | **5.0596 S** | 0.6856 S, below it |

"Clamp-free C" is the largest accepted forced conductance with **no shallow cell
on a bound** — the part of the band that owes nothing to where the limits were set.

**One case out of four has the truth above its clamp-free limit: 3P under
SEDIMENTARY, 0.2063 S against 0.1795 S, a 15% gap.** Near the truth the pinning is
light: one shallow cell at 0.2014 S (cell 24, 47–50 m) and three at 0.2260 S
(cells 16, 17 at 27–31 m and cell 24). The truth is comfortably inside the band
itself (min chi2 at C_true = 26.205 against 45.84) and comfortably inside the
relative rule too; it is only the no-cell-on-a-bound subset that stops short.

**This is direct evidence for the instruction not to adopt "no cell on a bound" as
a criterion**, and against the release path that still does. Were it a criterion,
3P's band would end at 0.1795 S and exclude the truth. See the blocker section in
`DECISIONS.md`.

## Where the pinning happens, and why

On 3P the pinned cells are not scattered. They sit in two bands, both outside the
forced window (7.5–23.5 m):

- cells 1–6, 0–9 m, the near surface above the window;
- cells 14–19, 22–36 m, immediately below the window base, plus cells 24–25 at
  47–53 m.

That is the compensation mechanism made visible. To carry more conductance in the
window and still fit 42 gates, the model must take conductance out above and below
it, and under SEDIMENTARY limits that compensation reaches the 0.001 S/m floor. At
3P's converged edge of 0.4320 S, 13 of 65 cells are on a bound and **all 13 are
shallow** — the edge is set by the limits, which is why it is reported as an edge
under stated limits rather than as a data-only bound.

On 2P under WIDE bounds and the frac = 0.02 split, nothing shallow is pinned
anywhere in the band: a band a factor of 126 wide that owes nothing to the clamps.
This is the one headline the split depth moves — at frac = 0.01 the cells at
193–208 m count as shallow and the clamp-free limit drops to 5.0596 S. Every other
headline in the table is identical at both fracs.

## The split depth is not robust at frac = 0.01, and it changes 2P's counts

`sensitivity_depth.jl` defines the sensitivity depth as the deepest cell whose
noise-weighted sensitivity reaches `frac` of the maximum. Both columns carry a
secondary bump in normalised sensitivity near 180–240 m, peaking at 0.0093 on 3P
and 0.0119 on 2P. On 3P that bump stays under 0.01; on 2P it crosses it. So at the
default `frac = 0.01`:

| frac | 3P split | cells excluded | 2P split | cells excluded |
|---|---|---|---|---|
| 0.05 | 76.9 m | 33 | 50.2 m | 41 |
| 0.02 | 80.7 m | 32 | 53.2 m | 40 |
| **0.01** | 84.6 m | 31 | **207.9 m** | **9** |

At 0.01 the ordering is physically wrong: 2P's conductive clay cap (0.35 S/m)
must *reduce* the depth of investigation relative to 3P's resistive section, not
quadruple it. A bump at 200 m sitting at ~1% of the peak is numerical, not
information.

The consequence for the counts is large. At 2P's node C = 0.0506 S under
SEDIMENTARY, the same 17 cells on a bound split as:

- frac 0.01 (split 207.9 m): **10 shallow**, 7 deep — including cells 30–38 at
  66–102 m, and under WIDE even cells 55–56 at 193–208 m, called "shallow";
- frac 0.02 (split 53.2 m): **1 shallow** (cell 22, 42–44 m), 16 deep.

On 3P the choice makes no difference at all: frac 0.01 and 0.02 give identical
tables, identical shallow counts and the same clamp-free limits.

**So this report uses frac = 0.02 for its headline numbers**, with the frac = 0.01
tables printed alongside in full so nothing is hidden. The definition itself was
not changed — `frac` is a documented parameter and both values are reported. A
robust definition (an envelope, or a first crossing on the decaying limb) is worth
having before the sensitivity depth is used for anything load-bearing, but it is
not needed for this report, because the one column where it matters (2P) has a
clamp-free band reaching 5.0596 S under either choice and a truth at 0.6856 S far
inside it.
