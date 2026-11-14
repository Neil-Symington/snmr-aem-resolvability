# The consistency check fires on the first column it runs on

Release 2P, config `synthetics/config/confined_aquifer.toml`, truth
`data/truth_murray2P.jld2` hash `33bf02a46442711e`, canonical window
17.5–33.5 m, bounds sigma 1e-4 to 10 S/m, lambda2 = 1.

```
CONSISTENCY: chiref 28.4587  <=  min over nodes 27.1749  (+1e-6) : **FAIL**
```

The **constrained** minimum, 27.1749 at C = 2.8452 S, is **1.28 lower** than the
**unconstrained** minimum, 28.4587. That cannot be true of a real optimum:
adding a constraint cannot improve the fit. So the unconstrained solve is not
the minimum on this column.

**Multi-start does not cure it.** The three unconstrained starts returned
28.459, 28.459, 28.459 — they agree to three decimals. The unconstrained
Gauss-Newton is robustly stuck at a point a constrained solve beats. The
constraint evidently conditions the problem better than the free one.

This is the same pathology the field work recorded
(`ek_realdata/REPORT_v1_firstfive.md:559-563`, chi2 30.5 unconstrained against
19.9 constrained), now reproduced on a synthetic column, and it is exactly what
the consistency assertion was added to catch.

## What it means for the two rules

The relative rule is `chi2 <= chiref + 3.84`. With chiref 1.28 too high the
threshold is 32.30 where it should be 31.02, so the relative band is **too wide**
— a looser bound, which is the conservative direction, but wrong.

It also vindicates anchoring the band to N: that choice was made precisely
because chiref could not be trusted, and here it cannot.

**Proposed fix, for M1.1:** define the reference as
`chiref = min(unconstrained multi-start, min over constrained nodes)`, which is
what the field tables already did in practice ("the tables quote the profile
minimum rather than the unconstrained value"). The assertion then becomes a
report of how far the unconstrained solve missed, rather than a failure.

## Two other results from the same pass

**The backward sweep wins the majority.** Winning starts over 11 nodes:
`warm_bwd` 6, `warm_fwd` 3, and the two direction seeds 1 each. Sweeping only
forward, as the current gate does, leaves chi2 on the table at more than half
the nodes. That is the both-directions change earning its cost.

**Neither crossing is inside the search range.** chi2 runs 29.17 down to 27.17
and back to 32.04 at C = 5.06 S, all far below the absolute threshold of 45.84,
so both edges are open on this grid and the upper crossing lies beyond
C = 5.06 S. The fine nodes for 2P therefore go **above** −0.5, not around the
retired column's band. The profile is smooth throughout — no node flagged.

Truth at 0.6856 S is inside both bands, comfortably.

## Cost note

lambda2 = 1 ran at **57 s per solve** against **2.8 min** at lambda2 = 0.01, so
the weak-smoothing solves are about three times dearer as well as less stable.
