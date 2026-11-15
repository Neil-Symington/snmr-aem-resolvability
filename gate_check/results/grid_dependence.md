# The recovered minimum depends on grid density — and that moves the edge

16 Nov 2026. 3P, config `synthetics/config/paper_3P.toml`, truth
`d9de86c430ab0165`, window 7.5–23.5 m, bounds sigma 1e-4 to 10 S/m.

Three runs, the same forced conductances, min chi2 in each:

| C (S) | 12-node warm | cold refits | 5-node bracket (0.05 dex) | best |
|---|---|---|---|---|
| 0.2845 | 26.089 | 26.088 | 26.088 | 26.088 |
| 0.3192 | — | — | 26.239 | 26.239 |
| 0.3582 | 29.075 | 29.075 | **28.341** | 28.341 |
| 0.4019 | — | — | 32.450 | 32.450 |
| 0.4509 | **47.878** | **47.876** | **42.943** | **42.943** |

Threshold N + 3.84 = **45.84**.

## What this says

**The denser grid found a better fit at two shared nodes**, by 0.73 at 0.3582 S
and by **4.94** at 0.4509 S. At the second of those the improvement crosses the
threshold: 47.878 is rejected, 42.943 is accepted. **So 3P's absolute upper edge
moves from "between 0.3582 and 0.4509 S" to "beyond 0.4509 S".**

**The cold refits did not find it.** At 0.4509 S the cold starts gave 47.876
against the warm 47.878 — a difference of 0.002. Starting far away, from `mref`
or from a uniform halfspace, does not help. What helped was a warm start from a
*closer node*: in the bracket the forward chain reached 0.4509 S via 0.4019 S,
0.05 dex away, instead of via 0.3582 S, 0.1 dex away.

**So the mechanism is step size, not start diversity.** The warm chain is a
continuation method: each solve begins at the previous node's answer, and the
smaller the step the better that answer serves. Grid density and solve quality
are the same knob, which I had been treating as independent.

## Consequences, including for what I proposed

**The band edge is grid-dependent, and a coarse grid errs the unsafe way.** A
sparse grid over-estimates chi2 near the edge, which makes the band narrower,
which is the direction that can exclude the truth. On 3P it did not matter —
truth at 0.2063 S is well inside either band — but the sign is the wrong one.

**The audit passed and validated the wrong thing.** Improvements at the three
audit nodes were 0.001, 0.026 and 0.000, worst 0.026, far under the 1.0
escalation trigger. That is a true statement about cold starts and says nothing
about grid density, because audit and warm came from the same grid and were
stuck in the same place. **The audit cannot detect this failure**, so its passing
is weaker evidence than I presented it as.

**The original gate's instinct was right.** Its `bisect_edge` solves at a new C
each step, warm-started, so it is dense by construction exactly where the edge
is — cheap and well-conditioned at once. Its defect was that it sometimes did
not run, not how it worked.

## Recommendation for M1.1

Locate each edge by **bisection warm-started from the nearest already-solved
node**, not from a fixed grid, and iterate until the edge stops moving by more
than a stated tolerance. Keep a coarse grid for the profile's shape and the
smoothness tests. That is the original design with the bracketing bug fixed, and
it addresses the grid dependence that neither extra lambda2 values, nor both
sweep directions, nor cold starts could reach.

Cost: about 12 warm solves per edge at the winning lambda2, against 20 solves
for the 5-node bracket that still left the edge unresolved.
