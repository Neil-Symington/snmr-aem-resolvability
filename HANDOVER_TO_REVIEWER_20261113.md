# Handover to the reviewer — 13 Nov 2026

From the writer's machine. Branch **`occam-workflow-test`** carries everything
(it descends from `sampler-convergence-test`). `git fetch && git checkout
occam-workflow-test`.

Nothing here has touched `numbers.csv`, `data/aem2P/`, `data/cond_C2_typical*`
or `outputs/synth_*`.

---

## 1. Two things only your machine can give

**a. Stage 1 model binaries, about 7 MB.**

```bash
tar czf /tmp/aem2P_models.tar.gz data/aem2P/models_* data/aem2P/points_*
```

`assembleTat1` needs these to rebuild sigma(z); the first collector deliberately
excluded them, so the autocorrelation time of Q2 and Q3 cannot be computed here.
That is the one outstanding item of the five you sent.

**b. `outputs/synth_2P/` from the 13 Nov production run.**

There is a sampler-side verdict table on this machine from 1–2 Nov, but it is a
different vintage: its Stage 5 TDS p50 is 5519 mg/L against your 4971. The Occam
comparison has to be against the run the paper traces to, not that one.

---

## 2. Where our numbers differ, and who is right

**Your R-hat figures reproduce exactly.** All five gate depths (1.664, 1.574,
1.235, 1.093, 1.957), the median 1.523, the maximum 1.961, and Q1 at 1.081. No
dispute on any of it.

**Your ESS figures were right and mine were optimistic.** `prerun_lib.jl:65`
computed an autocorrelation time *within* each chain, averaged them, and divided
total draws by the result — so it never saw the disagreement *between* chains.
Four chains mixing briskly around four different values each get a short IAT and
the sample looks large. That is exactly the case the gate exists to catch.

Replaced with the multi-chain rank-normalised estimator, Geyer paired monotone
sequence. On your chains:

| quantity | old | corrected | yours |
|---|---|---|---|
| Q1 window water content | 65 | **43** | **43** |
| theta at 5 m | 42 | 8 | 14 |
| theta at 15 m | 44 | 8 | 15 |
| theta at 100 m | 38 | 6 | 13 |
| all 60 depths, median | 44 | 9 | 15 |

Q1 now matches you exactly. The per-depth figures now come out *below* yours by
about a factor of two, so one of us is still conservative — worth settling,
though it changes no conclusion, since both are one to two orders below 100.
Validation: iid normal 4x10000 gives 40000 exact; AR(1) 0.90 gives 2244 against
a theoretical 2105; 0.99 gives 190 against 201.

**One framing of yours I think is too strong.** "Each ensemble has settled at a
different dimension and stayed there" is not what the chains do. The medians
9/16/20/24 reproduce exactly, but they are medians of broad, overlapping
distributions, and every ensemble crosses between k = 9 and k = 24 repeatedly:

| ensemble | p05 | p50 | p95 | IAT(k) | ESS(k) | crossings 9↔24 |
|---|---|---|---|---|---|---|
| e1 | 4 | 9 | 20 | 182 | 55 | 100 (1 per 100 draws) |
| e2 | 12 | 16 | 23 | 59 | 170 | 27 (1 per 370) |
| e3 | 4 | 20 | 35 | 308 | 33 | **1450 (1 per 7)** |
| e4 | 13 | 24 | 37 | 590 | 17 | 413 (1 per 24) |

What has not converged is the **occupancy**, not the reachability: e1 spends
5097 of 10001 draws at k <= 9 and e4 spends 268. This matters for the fix.
Tempering helps when modes are separated and hard to reach; here they are
reached constantly, and the slow thing is the weighting. It argues for reducing
`nmax` over adding ladder rungs. Figure: `fig_stage2_k_trace.png`.

**The conductance reconciliation does not close, and the truth is the tell.**
Truth is deterministic, so it pins the depth range. Over 4–34 m:

| interval | truth, native 0.25 m grid | truth, AEM grid |
|---|---|---|
| **4–34 m** | **5.60 S** | 5.84 S |
| 4–38 m | 7.54 S | — |
| 4–40 m | 8.84 S | 8.29 S |

So the stated truth of 8.0 S is not the 4–34 m figure; it implies about 4–38.7 m.
Your run's 6.65 S *is* right for 4–34 m, and 5.60 S sits inside its p05–p95 of
5.44–8.32. The paper's 9.1 S is not reachable over 4–34 m by any statistic of
that run: median 6.49, mean 6.56, integral of sigma50 4.89, integral of sigma90
14.4. Sensitivity there is 0.35–0.65 S per metre because 34–46 m is the
conductive confining clay. **Please check the integration limits the paper
actually used.** My reading is a labelling error, not a sampling difference.

**The window number is worse than a discrepancy and it is in both.** Truth over
17.5–33.5 m is 0.686 S. The paper says 2.32 S and your run 2.07 S — both about
**3x truth**. That matches your own bug report item 4 (AEM median sigma 0.133
against a true 0.046) and is the thin resistor smeared into the conductors
around it.

**The encouraging one.** Q1 against a true 0.2445 v/v: ensemble medians 0.2451,
0.2440, 0.2427, 0.2393; pooled 0.2433; truth at the 58th percentile of the
pooled posterior. Medians agree between ensembles to 2.4%, the 5th percentile to
7.3%. The decision quantity is accurate and reproducible even though per-depth
theta is not. That supports reframing the claims rather than discarding them.

---

## 3. Changes I made to shipped code — please review

| file | change | risk |
|---|---|---|
| `stages/stage3_hypothesis.jl` | `mctol = nd == 1 ? 0.0 : 1/nd` in three places | **identical arithmetic for nd > 1** |
| `tools/aem_HYP.jl` | `SCT_NMAX`, `SCT_TMAX`, `SCT_NCHAINSATONE` env overrides | defaults are the production values |
| `tools/quicklook.jl` | also writes `data/smr_occam<TAG>.jld2` | additive; pass/fail contract unchanged |
| `convergence_test/prerun_lib.jl` | ESS estimator replaced | vendored copy only |

The Stage 3 one is the one to look at hardest. `1/nd` is a Monte Carlo
tolerance; at one draw it evaluates to 1.0, the whole probability range, so the
tie tolerance swallowed every window and the stage reported "P(W >= 2 m) = 1 for
100.0% of 5709 windows" on the same page as "the worst window gives W(p5) =
0.460 m", which is below the 2 m threshold. The contradiction is what exposed
it. It would bite any very small ensemble, not only a deterministic one.

---

## 4. One decision that is yours, not mine

`ek_prerun/prerun_lib.jl` is the **original** of the ESS function and carries the
same defect. It is the copy wired into the live gate in `ek_prerun`,
`ek_realdata`, `udf_realdata` and `ssc_prerun`, all of which test
`ess_bulk_min >= GATE_MIN_ESS` with `GATE_MIN_ESS = 100`.

Those gates have therefore been testing against an estimate that can be several
times too high, and will have passed runs that should have failed. I have **not**
touched it, because correcting it changes past pass/fail outcomes and that is a
decision with re-runs attached. `ssc_prerun` already carries a second threshold,
`GATE_RHAT_PROD = 1.10`, commented as "the one to quote", which is a precedent
for how to handle it.

---

## 5. State of play

Running here: B1, the control — production settings, four independent
ensembles, 2P then 3P. 2P's e1 and e2 took 188 minutes each; e3 and e4 are in
flight.

Done: Occam on all four columns (2P, 2S, 3P, 3S), 208–242 s each, of which about
two minutes is solving. Every column reaches SNMR chi2/N 1.000 exactly and AEM
chi2/N 0.66–0.74. Stages 3–6 are being driven from those models through a
one-draw ensemble adapter, so the existing rules run unmodified.

Two corrections of my own, for the record. Test one is **3P**, not the `1P` tag I
invented from `unconfined_aquifer.toml` — that config differs materially
(`p_clay` 3.0 against 1.5, different variability, the retired `clay_mineral`
dropdown), with sigma_bulk differing by up to 0.79 S/m. Caught before B1's second
leg started. And the saline controls had no truth builder anywhere; the recipe
was recovered by differencing the old 2P and 2S truths and is a water-only change
in the target aquifer to `sigma_w = 0.30 S/m`, about 1,920 mg/L.
