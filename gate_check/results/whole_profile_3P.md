# What the rest of the profile has to become — 3P

From stored winning models in `edges_3P.jld2`. Unit boundaries are the TRUE
ones; units are named by role. Report only, as instructed: nothing here excludes
anything.

Truth, over the true boundaries:

| role | depth | mean sigma (S/m) | conductance (S) |
|---|---|---|---|
| cover/cap | 0–6 m | 0.00912 | 0.0547 |
| window | 6–40 m | 0.02689 | 0.9142 |
| below window | 40–58 m | 0.94246 | 16.9643 |
| base | 58–88 m | 0.12446 | 3.7338 |
| base | 88–200 m | 0.00773 | 0.8654 |

17 of 22 nodes accepted under the absolute rule.

## Near the truth the profile is sensible

Node 8, forced window C = 0.2845 S (truth over the forced window is 0.2063 S),
chi2 26.088, **no cells on bounds**:

| role | model conductance | truth | ratio |
|---|---|---|---|
| cover/cap | 0.0249 | 0.0547 | 0.46x |
| below window | 14.650 | 16.964 | 0.86x |
| base 58–88 | 4.779 | 3.734 | 1.28x |
| base 88–200 | 1.638 | 0.865 | 1.89x |

Everything within a factor of two, nothing pinned.

## At the top of the band the profile becomes implausible

Node 17, forced window C = 0.4509 S, chi2 42.766, **10 cells on the low bound
and 4 on the high bound**:

| role | model conductance | truth | ratio |
|---|---|---|---|
| **cover/cap** | **0.0006** | 0.0547 | **0.011x** |
| below window | 17.507 | 16.964 | 1.03x |
| base 58–88 | 2.445 | 3.734 | 0.65x |
| **base 88–200** | **28.620** | 0.865 | **33x** |

The cover collapses to **0.00010 S/m — the parameter floor — in all five of its
cells**, 91 times below truth and pinned there. Inside the window, five or six
cells also sit at the floor while others reach 0.147 S/m, so the forced
conductance is met by a spiky profile rather than a raised one. And the deep
basement carries **28.6 S against a true 0.87 S**.

The same pattern holds across nodes 14–17, the whole approach to the edge: cover
at the floor with 5 of 5 cells pinned, 9–10 cells on bounds overall.

## At the bottom of the band, the mirror image

Node 2, forced window C = 0.0090 S, chi2 26.297, **5 cells on the high bound**:
the basement 88–200 m carries **426 S against a true 0.87 S**, a factor of 490,
with five cells at the 10 S/m ceiling.

## What this means for reading the band

The upper edge at 0.4594 S is a data limit in the sense that chi2 genuinely
crosses the threshold there. But the models that reach the top of the band have
**ten parameters pinned at their clamps**, so the band's upper end is as much a
statement about where the parameter bounds were set as about what the data
allow. The same is true at the bottom.

That does not change any verdict — truth at 0.2063 S sits in the middle of the
band, in the region where no cell is pinned and every unit is within a factor of
two of truth. But a band edge should not be quoted without the count of pinned
cells beside it, and M5's tables should carry that column.

---

## Pinned-cell count per accepted node (M1.1 close-out item 2)

Unit conductances as a ratio to truth. Bounds sigma 1e-4 to 10 S/m.

| C (S) | min chi2 | pinned | rule | cover | window | below | base 58–88 | base 88–200 |
|---|---|---|---|---|---|---|---|---|
| 0.0090 | 26.297 | 5 | both | 2.59x | 1.78x | 0.89x | 1.80x | **492x** |
| **0.2845** | 26.088 | **0** | both | 0.45x | 0.82x | 0.86x | 1.28x | 1.89x |
| 0.3192 | 26.239 | 4 | both | 0.29x | 1.46x | 0.90x | 1.68x | **169x** |
| 0.3582 | 28.341 | 7 | both | **0.01x** | 0.66x | 0.90x | 1.19x | 0.70x |
| 0.3634 | 28.722 | 9 | both | 0.01x | 0.64x | 0.91x | 1.18x | 0.75x |
| 0.3687 | 29.048 | 9 | abs | 0.01x | 0.61x | 0.91x | 1.17x | 0.79x |
| 0.3794 | 29.941 | 9 | abs | 0.01x | 0.56x | 0.93x | 1.15x | 0.63x |
| 0.3905 | 31.440 | 9 | abs | 0.01x | 0.51x | 0.94x | 1.12x | 0.78x |
| 0.4019 | 32.450 | 9 | abs | 0.01x | 0.45x | 0.95x | 1.08x | 0.83x |
| 0.4257 | 37.730 | 10 | abs | 0.01x | 0.47x | 0.97x | 1.08x | 12.1x |
| 0.4509 | 42.766 | 14 | abs | 0.01x | 0.50x | 1.03x | 0.65x | **33.1x** |

## The unpinned band is one node wide

**Exactly one accepted node has no cell on a bound: C = 0.2845 S**, chi2 26.088.
At that node every unit sits within a factor of 2.2 of truth, from 0.45x on the
cover to 1.89x on the deep base. Nothing is pinned and nothing is absurd.

Every other accepted node carries 4 to 14 pinned cells. Above 0.3582 S the cover
is at **0.01x of truth** — the global floor — at every single node.

Ceiling at that unpinned node:

| | widest | central | tightest |
|---|---|---|---|
| C = 0.2845 S, prior theta >= 0.25 | **1585** | 769 | 450 |
| C = 0.2845 S, SNMR theta >= 0.265 | **1419** | — | — |

So requiring no cell at a clamp — with no geological limits at all — already
brings the widest-corner ceiling from **2560 to 1585 mg/L**, well clear of both
the 3,000 and the 1,500 mg/L lines.

**One caveat on that number.** "Largest C with no pinned cell" is limited by the
grid: the next node up is 0.3192 S, which has 4 pinned cells, and nothing was
solved between the two. So 0.2845 S is a **lower bound** on the largest unpinned
conductance, not the value itself. Resolving it would need a few solves in that
gap, and it would only move the ceiling upward from 1585.

The largest C accepted under the **relative** rule is 0.3634 S, with 9 pinned
cells — so the relative rule's edge is also in the pinned region.
