# Overnight status

Started 2026-11-13T22:01:18+11:00. This file is rewritten as each step finishes.

## EARLY RESULT — B1 test two (2P) diagnostics, run ahead of the queue

All four 2P ensembles finished clean (187, 188, 202, 202 min). I ran the
diagnostics immediately rather than leaving them to the queue, so that a failure
in 3P overnight could not leave you with nothing.

**This is the first between-ensemble R-hat Stage 1 has ever had.** Production
ran one ensemble with `nchainsatone = 1`, so only a segment R-hat was possible.

### The decision quantities

| | ens 1 | ens 2 | ens 3 | ens 4 | spread | R-hat | ESS bulk | gate |
|---|---|---|---|---|---|---|---|---|
| **Q2** window conductance, p50 (S) | 1.946 | 1.879 | 1.972 | 1.597 | **1.23x** | 1.127 | 41 | FAIL |
| Q2 p05 | 0.474 | 0.712 | 0.826 | 0.717 | **1.74x** | | | |
| **Q3** package conductance, p50 (S) | 6.420 | 6.765 | 6.696 | 6.457 | **1.05x** | 1.120 | 36 | FAIL |
| Q3 p05 | 5.430 | 6.010 | 5.811 | 5.812 | 1.11x | | | |

Four things come out of this.

**It fails, but mildly, and about as much as Stage 2 does.** R-hat 1.12–1.13
against a 1.05 threshold, ESS 36–41 against 100. Stage 2's Q1 was R-hat 1.081,
ESS 43. The two stages are in the same condition on the quantities that decide
anything — not one healthy and one broken.

**The sampler's mechanics are healthy, which rules out the obvious diagnosis.**
Birth 38–43% and death 39–43% accepted at T = 1; k changes on 90–93% of steps;
k spans 5–40, the entire prior; node-count R-hat is **1.041, which passes the
gate**. Nothing is stuck. The model parameters move freely and it is the derived
*integrals* that decorrelate slowly. More tempering rungs do not obviously
address that.

**The integral is four to five times more reproducible than the window.** Q3
medians spread 1.05x, Q2 medians 1.23x and Q2's lower band edge 1.74x. That is
the same effect already on the record for detection — integrate over a package
and the answer is stable, narrow the window and it is not.

**B1 reproduces the production run.** Production's single ensemble gave Q3
6.647 S; the four independent ensembles give 6.420–6.765. Its Q2 was 2.073 S
against 1.597–1.972 here, slightly above all four, which is about what an
unconverged quantity with a 1.23x spread should look like. The production number
was not a freak realisation.

### Against truth

Truth is 5.602 S for the package and 0.686 S for the window.

- Package: every ensemble sits **15–21% high**.
- Window: every ensemble sits **2.3–2.9x high**, confirming the smearing of the
  thin resistor into the conductors around it.

So the window conductance is wrong by a factor of three in a way no amount of
sampling fixes, and that error dwarfs the 1.23x reproducibility problem.

### One per-depth result that looks like good news and is not

| depth | R-hat | ESS bulk | gate |
|---|---|---|---|
| 5 m | 1.274 | 19 | FAIL |
| 14 m | 1.303 | 13 | FAIL |
| 29 m | 1.072 | 60 | FAIL |
| 59 m | 1.057 | 76 | FAIL |
| 102 m | 1.008 | 1445 | **PASS** |

The only depth that passes is the deepest. The ensembles agree at 102 m because
the AEM data do not constrain it and all four return the same prior-dominated
answer. That is agreement about nothing, and it should not be read as the deep
section being well resolved.

