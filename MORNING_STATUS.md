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

---

## OCCAM vs SAMPLER — test two (2P), complete verdict

The Occam site run finished. Sampler column is the 1–2 Nov run (the 13 Nov
production `outputs/` is still on the reviewer's machine).

| row | sampler | Occam | factor |
|---|---|---|---|
| Interval: top | 17.00 m (p5–p95 15.75–42.00) **multimodal** | **17.25 m** (17.25–17.25) | — |
| Interval: base | 33.00 m (p5–p95 31.75–58.00) | **33.25 m** | — |
| Interval: length | 16.00 m, pinned at L_max | 16.00 m, pinned at L_max | 1.00x |
| Windows passing water criterion | 7.6% | **47.3%** | 6.2x |
| Admissible conductance band | 0.03192 – 7.079 S | 0.03049 – 8.82 S | lower 1.05x, upper **1.25x** |
| Band width as a factor | 222x | 289x | 1.30x |
| Data-only bound (clay-free) | 49,... see note | 49,256 mg/L | 0.82x |
| Stage 5 acceptance rate | 77.4% | 76.6% | 1.01x |
| Stage 5 estimate, TDS p50 | 5,519 mg/L (816–14,157) | **4,238 mg/L** (1,376–7,308) | 0.77x |
| Verdict: interval | PASS (multimodal) | **PASS** | same |
| Verdict: water | PASS (uninformative) | **PASS (uninformative)** | same |
| Verdict: salinity | UNDETERMINED | **UNDETERMINED** | same |
| **Verdict: OVERALL** | **UNDETERMINED** | **UNDETERMINED** | **same** |

(The sampler's data-only bound was 59,926 mg/L.)

### What this says

**Every verdict is identical**, including all three sub-verdicts, and the
quantitative rows differ by **1.2–1.3x**. On this column the trans-dimensional
sampler is not what produces the answer.

**The located window agrees to a quarter of a metre** — 17.25–33.25 m against
17.00–33.00 m, with truth at 18–34 m. Occam finds the aquifer.

**What the sampler actually buys is the bound, not the location.** It returns a
5th-to-95th range on the interval top of 15.75–42.00 m and flags it
**multimodal** (12.2–20.0 m at 95%, 42.0–44.0 m at 5%). Occam returns a point
and cannot say anything about either. Same for the water criterion: 47.3% of
windows pass against 7.6%, because P(W >= w) is driven by the 5th percentile of
W, and with no posterior spread that collapses onto the point value.

**One row to treat with care.** The Stage 4 floor is zero in both, but for
different stated reasons — the sampler says the band's lower edge is open, the
Occam run says there is no `stage4_floor` record for the column, because
`stage4_floor.jl` is not in the `run_site.jl` chain. That row is a non-result in
both and should not be read as agreement.

### Caveat on this table

The sampler column is a 1–2 Nov run, not the 13 Nov production one. Its Stage 5
TDS p50 is 5,519 mg/L against production's 4,971 mg/L, so the vintages differ by
about 1.11x on their own — comparable to several of the differences above.
**This comparison needs `outputs/synth_2P/` from the reviewer's machine before
any of these factors should be quoted.**

## B1 (the control) — FINISHED

```
2P ENSEMBLE_DONE e=1 rc=0 minutes=187
2P ENSEMBLE_DONE e=2 rc=0 minutes=188
2P ENSEMBLE_DONE e=3 rc=0 minutes=202
2P ENSEMBLE_DONE e=4 rc=0 minutes=202
3P ENSEMBLE_DONE e=1 rc=0 minutes=188
3P ENSEMBLE_DONE e=2 rc=0 minutes=188
3P ENSEMBLE_DONE e=3 rc=0 minutes=185
3P ENSEMBLE_DONE e=4 rc=0 minutes=186
```

## B1 diagnostics, 2P (exit 0)

```
1              0.4736     1.9461     3.0929
2              0.7121     1.8790     3.0465
3              0.8262     1.9718     3.4940
4              0.7168     1.5970     2.4536

-- Q3 package conductance (S) --
ensemble          p05        p50        p95
1              5.4299     6.4202     7.5398
2              6.0103     6.7651     7.6480
3              5.8106     6.6960     7.5470
4              5.8122     6.4571     7.2366

-- convergence across ensembles --
quantity                        R-hat   ESS bulk   ESS tail med spread     pass
Q2 window conductance           1.127         41         98       1.23x     FAIL
Q3 package conductance          1.120         36        154       1.05x     FAIL

-- per-depth log10 conductivity (reported separately; not decisive) --
depth                           R-hat   ESS bulk   ESS tail     pass
5 m                             1.274         19         15     FAIL
14 m                            1.303         13         28     FAIL
29 m                            1.072         60        230     FAIL
59 m                            1.057         76        588     FAIL
102 m                           1.008       1445      44141     PASS

-- dimension mixing at T=1 (the mechanism) --
ensemble       median k      k range k changes/step
1                  26.0         5-40        93.053%
2                  19.0         6-40        90.205%
3                  19.0         5-40        92.104%
4                  20.0         5-40        92.354%
node count across ensembles: R-hat 1.041, medians [26.0, 19.0, 19.0, 20.0], spread 1.37x

-- windowed acceptance rates, %% (mean over stat windows) --
ensemble       birth     death  position  property        dc   scope
1              41.90     42.13     49.32     69.31     49.10   all T
               42.85     43.17     48.87     68.25     53.33   T=1 windows
2              37.88     38.63     49.53     70.81     27.94   all T
               37.81     38.55     48.54     70.78     25.05   T=1 windows
3              41.14     41.35     49.87     70.61     42.54   all T
               42.26     42.44     48.71     69.28     46.94   T=1 windows
4              39.22     39.42     49.46     72.74     29.37   all T
               39.90     39.65     49.22     72.18     29.15   T=1 windows

wrote /home/symo/Documents/GA/AEM_SMR_Resolvability_Quickrun/release/convergence_test/results/stage1_2P_b1.csv
```

## B1 diagnostics, 3P (exit 0)

```
1              0.0964     0.1757     0.2466
2              0.0659     0.1374     0.2365
3              0.0616     0.1522     0.2472
4              0.0807     0.1642     0.2225

-- Q3 package conductance (S) --
ensemble          p05        p50        p95
1              0.3553     0.5041     0.6533
2              0.3389     0.4980     0.7562
3              0.3905     0.5613     0.7817
4              0.2701     0.4477     0.6321

-- convergence across ensembles --
quantity                        R-hat   ESS bulk   ESS tail med spread     pass
Q2 window conductance           1.038         94        183       1.28x     FAIL
Q3 package conductance          1.089         49        116       1.25x     FAIL

-- per-depth log10 conductivity (reported separately; not decisive) --
depth                           R-hat   ESS bulk   ESS tail     pass
5 m                             1.064        159         62     FAIL
14 m                            1.031        179        172     PASS
29 m                            1.019        194        296     PASS
59 m                            1.044         63        103     FAIL
102 m                           1.015        927        746     PASS

-- dimension mixing at T=1 (the mechanism) --
ensemble       median k      k range k changes/step
1                  18.0         4-40        92.454%
2                  22.0         6-40        92.754%
3                  17.0         4-40        92.004%
4                  18.0         4-40        92.004%
node count across ensembles: R-hat 1.041, medians [18.0, 22.0, 17.0, 18.0], spread 1.29x

-- windowed acceptance rates, %% (mean over stat windows) --
ensemble       birth     death  position  property        dc   scope
1              40.13     40.99     47.27     62.67     42.41   all T
               39.64     40.19     46.53     61.48     40.20   T=1 windows
2              41.32     41.66     47.10     63.33     46.13   all T
               40.99     41.53     46.75     63.08     44.72   T=1 windows
3              40.85     41.35     47.11     60.98     46.59   all T
               40.43     41.15     46.04     60.64     46.24   T=1 windows
4              41.25     41.64     47.55     61.77     43.22   all T
               40.55     40.63     46.33     61.44     42.35   T=1 windows

wrote /home/symo/Documents/GA/AEM_SMR_Resolvability_Quickrun/release/convergence_test/results/stage1_3P_b1.csv
```

## Occam site run, 2S (exit 1, 103 min)

```
"occam_2S","2S","2P","saline_control","Interval: top","17.25 m  (p5-p95 17.25 - 17.25)","unimodal"
"occam_2S","2S","2P","saline_control","Interval: base","33.25 m  (p5-p95 33.25 - 33.25)",""
"occam_2S","2S","2P","saline_control","Interval: length","16.00 m  (p5-p95 16.00 - 16.00)","PINNED at L_max in 100% of draws: spec-determined, not inferred"
"occam_2S","2S","2P","saline_control","Water criterion: windows passing","45.6% of 5709 allowed windows","the discriminating statistic; P for the best window sa
"occam_2S","2S","2P","saline_control","Water criterion","P(W >= 2.0 m) = 1.0000","45.6% of 5709 allowed windows satisfy it; SATURATED -- not a pass earned"
"occam_2S","2S","2P","saline_control","Admissible conductance band","0.03134 - 8.82 S, width 2.449 dex","lower edge open = false; upper edge open = false"
"occam_2S","2S","2P","saline_control","Verdict: interval","PASS",""
"occam_2S","2S","2P","saline_control","Verdict: water","PASS (uninformative)",""
"occam_2S","2S","2P","saline_control","Stage 4 floor under stated ceilings","TDS >= 0 mg/L","zero: no stage4_floor record for this column"
"occam_2S","2S","2P","saline_control","Verdict: salinity","UNDETERMINED","clay-free bound cannot clear the limit; clay and salt are not separable here, so this 
"occam_2S","2S","2P","saline_control","Verdict: OVERALL","UNDETERMINED",""
```

