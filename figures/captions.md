# Figure captions — tier-P paper set

One line per panel, for writing final captions from. All figures are from the
tier-P runs on the **P lineage** (3P, 2P, 3S, 2S) unless stated.

**Conventions.** PNG at 300 dpi; 170 mm wide for two-column figures, 85 mm for
one-column. Panel labels (a), (b), … top-left. Fresh column **blue**, saline
control **orange**, truth **black**, manager lines **grey dashed**. *Exception:*
on the black-background density figures (3 and 4) the truth is drawn in a bright
**cyan-blue** and the posterior summary in white — a black line on a black field
is only its halo, and the brief's blue is too dark against inferno's low end. Fonts ≥ 8 pt
at print size. No suptitles; no truth values in any title — truths are drawn as
a line and named in the legend.

**These supersede the October set**, which used a different convention (test one
red, test two blue, truth cyan) and is archived unchanged in
`superseded_2026-10/`. One script per figure, `figNN_<slug>.jl`, sharing
`figstyle.jl`. Every value is in `numbers.csv`, which carries a `figure` column.

---

**Figure 1 — workflow diagram.** Not produced here; drawn by NS from Section 3.

**Figure 2 — `fig02_columns.png`.** The two synthetic columns, row per column
(top: test one, unconfined; bottom: test two, confined), all panels on one depth
axis 0–60 m, aquifers shaded.
(a, e) total porosity and SNMR water content. (b, f) clay fraction.
(c, g) pore-water salinity against the 1,000 mg/L line, fresh and saline control
overlaid. (d, h) bulk conductivity as a stack: the clay term as the shared grey
base, each column's water term on top of it to the total.
*The clay term is the larger part of σ_bulk on both columns, and it is identical
in a column and its control — only the water term moves.*

**Figure 3 — `fig03_aem_posteriors.png`.** The full Stage 1 AEM conductivity
posterior as a 2D density, 0–80 m, log σ in S/m, colour the number of draws per
cell on a log scale, **inferno on a black field**, so an empty cell and a
barely-sampled one read alike. **Depth bins are the model cells**, so the density
is the sampler's and not a smoothing artefact. p10, p50 and p90 over it as one
white family (p50 solid, p10/p90 dashed), truth in cyan-blue, aquifer tinted
white with dashed green boundaries.
*The percentile lines are a readout of the posterior, not models: pointwise
percentiles are depth-wise order statistics and are never propagated into any
downstream forward in this workflow.*
A grey wash is laid **over the whole panel — background, density and percentile
lines together** — in proportion to how little a cell at that depth could move
the data: transparent where Δlog σ for one noise σ is ≤ 0.1 dex, opaque at
≥ 1 dex, linear between. *Covering the background is deliberate: dimming only
the coloured density would make an unresolved depth look identical to a depth
where the posterior simply has no mass, which is the opposite of the intended
reading.* The weight comes from the operator's own Jacobian at the true model,
so nothing about the sampler enters it; the green sidebar plots it directly as
"resolved", 0 to 1. (a) test one, (b) test two.
*The density fades into grey exactly where the data stop constraining it. On
test two the wash is darkest over the aquifer itself — the visual counterpart of
the 2.66× conductance bias reported in Section 5.3.*

**Figure 4 — `fig04_snmr_posteriors.png`.** The full Stage 2 SNMR
water-content posterior as a 2D density, built to read exactly as Figure 3:
depth bins on the model cells, inferno on a black field, p10/p50/p90 as one
white family, truth in cyan-blue, aquifer tinted white with dashed green
boundaries. No detectability sidebar: the wash carries it, and the sidebar only
added clutter. Legend below the panels, clear of the data.
The wash is a real detectability curve, not a ramp below a line: per cell, the
smallest θ change the instrument could register, dmin = 1/√Σ_q(w·|K|/σ̂)²,
transparent at dmin ≤ 0.05 v/v and opaque at dmin ≥ 0.45 v/v (the θ prior's
ceiling). It covers the whole panel including the lines, for the same reason as
Figure 3. The depth where dmin crosses 0.05 — the **sensitivity DOI**, 59.2 m on test one
and 49.4 m on test two — is **not drawn as a line**: the wash already fades the
panel out there and a line on top of a fade states it twice. The values are in
`numbers.csv`.
*This is the sensitivity DOI (Christiansen & Auken: cumulative noise-normalised
kernel sensitivity, the depth where a Δθ = 0.05 v/v change stops moving the data
by 1σ), never the trans-D posterior-width DOI, which tracks nucleus count rather
than depth sensitivity and must not be quoted.* (a) test one, (b) test two.
***The depth limit is per column, not per tier:* 59.2 m on test one against
49.4 m on test two, both at tier "typical", because the SNMR kernel is
conductivity-sensitive and each column's own conductivity sets it. An earlier
draft of this caption gave 70.1 m as a tier-level number; that figure belongs to
the bare lineage on an 88.6 m loop and is not these columns.**

**Figure 5 — `fig05_interval.png`.** What Stage 3 chooses, row per column.
(a, d) histograms of the window top and base over all draws, **log counts**.
(b, e) histogram of W for each draw's best window, with the W ≥ 2 m line, fresh
and saline control. (c, f) the window-pass fraction as a bar, fresh beside its
control.
*Test two's window posterior is bimodal — an upper mode at 12–20 m on 95% of
draws and a lower one at 42–44 m on 5% — which a p50 with a p5–p95 hides. The
pass fraction is 54% on test one against 8% on test two; P for the best window
saturates at 1.0 on all four columns and discriminates nothing.*

**Figure 6 — `fig06_clay_salinity_tradeoff.png`.** The clay–salinity
trade-off. (a) test one, (b) test two: accepted draws as (clay term, TDS) pairs
with the deterministic trade-off curve through them, the 1,000 mg/L line
horizontal and the truth marked. The clay term is the clay conductivity times
its volume fraction, the quantity the Stage 0 prior is on.
*Every point on the curve fits the AEM data equally well; it is the clay prior,
not the data, that decides where on it the answer sits.*

**Figure 7 — `fig07_gate_sigmin1e-3.png`.** The gate, two panels on one
conductance axis. **(a)** χ² of the forced fit against window conductance for
the two fresh columns, interpolated through the grid nodes the gate solved at.
The solid line is N + 3.84 = 45.84 with the region above it excluded by the
data; the dotted line is N = 42. The open circle marks the true window
conductance on each curve. Saline controls are not drawn. **(b)** the range of
window conductance the data admit, one bar per column, controls faded. An
arrowhead marks an end set by the parameter bounds rather than by the data —
χ² never reaches N + 3.84 there — and a flat end an end the data close. The open
circle is the true window conductance. The three ticks are 1,000, 1,500 and
3,000 mg/L converted to window conductance at that column's median water
content and cementation exponent, with no clay term; they are labelled on the
two test columns and repeated unlabelled on their saline controls.

The file name carries σ_min = 1e-3 S/m because that is the run the edges come
from; the σ_min = 1e-4 re-run has not been made, and the targeted open-edge walk
found no χ² crossing above either upper edge, so a lower floor can only widen
these bars. The profile λ² is 1 on the test-one columns and 10 on the test-two
columns, where λ̂ = 100 and no λ² = 1 sweep exists; each curve is labelled with
its own λ². Per-node χ² values are in `numbers_chi2.csv`.
*Test one's χ² never rises to the exclusion line, so both its ends are limits of
the model rather than measurements; its bar nonetheless stops between the 1,000
and 1,500 mg/L ticks — fresh at 1,500 and 3,000, undetermined at 1,000. Test
two's χ² does reach the line, closing its upper end at 7.08 S, which is past
3,000 mg/L: all three thresholds undetermined.*

**Figure S8 — `figS8_saline_control_curves.png`.** L(C) at λ̂ for each column
and its saline control, with both truths and the 1,000 mg/L line in conductance.
(a) test one, (b) test two.
*Test one separates completely, with the 1,000 mg/L line falling between the two
curves. Test two does not separate at all: the two curves lie on top of each
other, both to the right of the line.*

**Figure 8 — `fig08_priors_narrowed.png`.** Narrowing the Stage 0 priors with
borehole logs, test two only. One panel, six rows on a log TDS axis. The top
four are the prior states at the bore — survey only, gamma, induction, both.
After a gap, the transfer test: the clay term fixed at its true value for the
window and carried to a survey-only neighbour, fresh and saline. Each bar is
the Stage 5 TDS p5–p95, the filled circle its p50, the open ring the truth
(561 mg/L, and 1,736 mg/L on the saline-neighbour row), and the dashed vertical
the 1,000 mg/L line. P(TDS < 1,000) is printed at the right of each bar. All
bars are grey because nothing on this column is certified; the orange PASS rule
is kept in the code for the field version.

The gamma log replaces the Stage 0 clay prior with a log-normal centred on the
true window clay term, 95% interval ×/÷ 2.5. The induction log is applied as a
likelihood factor on each Stage 1 draw's window conductance, centred on the
truth, ×/÷ 1.30 — a reweighting of the survey's own draws, not a replacement.
Stage 0 has no prior on clay fraction, only one on the product σ_clay·C, so that
is what the gamma log narrows. Every Stage 3 window is exactly 16.00 m thick,
so a likelihood on the window's mean conductivity and one on its conductance
are the same function. The survey-only row reproduces the Stage 5 result of
Section 5.2 exactly, being the same 40,000 draws.
*Neither log certifies this column and nor do both together: P(fresh) reaches
only 0.71 against the 0.95 the rule requires. The two transfer rows are the
sharper result — a clay term fixed at exactly the right value puts the fresh
neighbour at 4,883 mg/L against a truth of 561 and the saline one at 4,711
against 1,736, so the two are indistinguishable. What defeats the neighbour is
the survey's own window conductivity, not the clay.*

**Figure 9 — field pair.** `bon_mid_4` and `SMR_UD24`. **Not yet produced** —
adapts `pair_report/pair_figures2.jl` (`figsite`) with suptitles removed.

**Figure 10 — conductance in the field.** **Not yet produced** — the pair
report's Figure 3, `pair_report/pair_figures.jl` (`figP3`).

**Figure 11 — the clay prior in the field** (85 mm). **Not yet produced** —
`pair_report/pair_figures2.jl` (`figclay`), to carry the acceptance rate in the
legend and the log-uniform product prior drawn.

---

**S1–S3.** As they exist, in `superseded_2026-10/`.

**S4 — dropped.** The bound-decomposition waterfall is withdrawn: it walked
from the old gate to the new one and so mixed a method change with a column
change. In its place, three sensitivity factors for test one, computed on the
final run with everything else held fixed and one unknown swapped at a time
(`synthetics/sensitivity_factors.jl`). Against the data-only bound of
1,451 mg/L: clay term known **0.57×** (822 mg/L), porosity known **0.74×**
(1,079 mg/L), m known **0.83×** (1,200 mg/L). No figure.

**S7 — dropped.** The profile and Laplace curves are no longer needed now that
panel 7(a) carries χ²(C) in the main text.

**S5 — tier-F certified bound against truth per realisation.** Pending the
tier-F runs.

**S6 — test-two MCMC against the Laplace.** Pending the validation run.
