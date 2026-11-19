# Band edges, 3S and 2S — 20 Nov 2026

Window 7.5–23.5 m (3S) and 17.5–33.5 m (2S), both 16.000 m thick. Stated limits
0.001–2 S/m, λ² set {0.01, 0.1}, both sweep directions. Absolute rule
χ² ≤ N + 3.84 = 45.840; relative rule χ²ref + 3.84.

Artefacts: `edges_3S_sed.jld2` (upper), `edges_3S_sedlow.jld2` (lower, after the
reflection fix `edge-lower-reflect-20261119b`), `edges_2S_sed.jld2`.

## 3S — test one, brackish

| rule | lower edge | upper edge | admissible range |
|---|---|---|---|
| absolute | **0.2950 S** converged | **1.7944–1.7964 S** not converged | 0.2950 – 1.7944 S |
| relative | **0.4817 S** converged | **1.6641 S** converged | 0.4817 – 1.6641 S |

Truth conductance over the window is **1.1246 S**, INSIDE under both rules.

The absolute upper edge is reported as the bracket from the solved nodes
(highest accepted −0.950195 dex, χ² 45.414; lowest rejected −0.949707 dex,
χ² 45.860), per Neil's ruling of 19 Nov. The bracket is 0.0005 dex wide, so the
edge is pinned to about 0.1%; NOT CONVERGED refers to the χ²-stability
criterion, not to bracket width.

## 2S — test two, brackish

All four edges OPEN. No rejected node exists on either side of the accepted
range within the stated limits, so neither edge is bisected and no ceiling can
be stated. Truth conductance 1.3747 S. The gate does not constrain this column.

## Salinity ceilings

Clay-free ceiling, all of the window's conductance given to pore water:
σ_w = (C/L)/φ^m, TDS = σ_w × 6400. Porosity 0.25–0.40, m 1.50–1.90.

| case | C | least cautious | central | most cautious |
|---|---|---|---|---|
| 3S absolute, upper | 1.7944 S | 2,837 | 4,850 | 9,998 mg/L |
| 3S relative, upper | 1.6641 S | 2,631 | 4,498 | 9,272 mg/L |
| 3S absolute, lower | 0.2950 S | 466 | 797 | 1,644 mg/L |
| 3S relative, lower | 0.4817 S | 762 | 1,302 | 2,684 mg/L |
| 2S, bound not edge | 32 S | 50,596 | 86,498 | 178,289 mg/L |

3S's true window pore water is 2,707 mg/L. Every ceiling listed for 3S is
clay-free, so it is an upper bound on salinity, not an estimate; the window's
mean clay fraction is 0.058.

## Does 3S's allowed range overlap test one's?

Test one (3P) has an OPEN lower edge, so its range runs from the floor to its
upper edge.

| rule | 3P range | 3S range | overlap |
|---|---|---|---|
| absolute | 0 – 0.4320 S | 0.2950 – 1.7944 S | **YES**, 0.2950 – 0.4320 S |
| relative | 0 – 0.3605 S | 0.4817 – 1.6641 S | **NO**, disjoint by a factor of 1.34 |

**The rule decides it.** Under the absolute rule the fresh and brackish columns
of test one are not separated by window conductance alone: a conductance between
0.2950 and 0.4320 S is admissible for both. Under the relative rule they are
disjoint, with a gap from 0.3605 to 0.4817 S.

## Caveats

1. The two rules' upper and lower edges were computed in different runs with
   slightly different χ²ref: 25.3500 for the upper (threshold 29.190) and
   25.2903 for the lower (29.130), because the lower run was seeded with the
   larger solved set and that set contains a lower minimum. The relative-rule
   range therefore mixes two thresholds 0.06 apart. The absolute rule is
   unaffected, its threshold being fixed at N + 3.84.
2. Three warnings in the lower-edge run where a warm re-solve RAISED a node's
   χ²: node −1.75 48.209 → 48.993 (absolute); nodes −1.625 37.786 → 38.559 and
   −1.5 27.922 → 28.440 (relative). Same family as the five logged for the
   upper edge. Stored minima are never overwritten by a worse solve, so the
   edges stand.
3. The refinement audit passed: 3 intervals, improvements 0.012, 0.028, −0.014,
   none above the 0.5 trigger.
