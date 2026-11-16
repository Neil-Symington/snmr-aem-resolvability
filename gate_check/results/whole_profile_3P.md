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
