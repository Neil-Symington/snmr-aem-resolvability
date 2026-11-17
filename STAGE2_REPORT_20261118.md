# Stage 2, columns 3S and 2S — 18 Nov 2026

Data from tag `snmr-data-truth-20261118` (generating kernel on the true
conductivity). Inversion kernel from the pooled median of each column's four
Stage 1 ensembles, 8008 post-burn-in T=1 draws each. 100,001 samples, 4 PT
chains, 4 ensembles, Tmax 2.5, burn-in 0.5. 24 moments, σ̂ median 4.63 nV.
Run in 11 min per column. Reproduce with `tools/stage2_report.jl --tag <TAG>`.

**Nothing has been tuned in response to any of this.**

## The four numbers

| | 3S (window 7.5–23.5 m) | 2S (window 17.5–33.5 m) |
|---|---|---|
| **(3) kernel error in data space** | **median 7.96σ, max 41.06σ, 79.2% of 24 above 3σ** | **median 3.23σ, max 20.27σ, 54.2% above 3σ** |
| … in amplitude | median 3.05%, max 13.01% | median 4.24%, max 36.93% |
| **(1) inferred noise factor s, p50 per ensemble** | **0.839, 1.775, 1.991, 1.816** | **1.152, 1.416, 1.265, 1.290** |
| … p05–p95 across ensembles | 0.51–4.45 | 0.65–2.34 |
| **(2) chi2 per datum vs the TRUE noise, p50** | **2.96, 1.54, 1.81, 1.95** | **2.12, 0.99, 0.97, 2.23** |
| **(4) window water content, pooled p50** | **0.2710 v/v** (truth 0.2704) → **1.002×** | **0.2474 v/v** (truth 0.2438) → **1.015×** |
| … truth inside pooled p05–p95 | **yes** (0.2631–0.2763) | **yes** (0.2358–0.2618) |
| … per-ensemble p50 / truth | 1.014, 1.007, 0.981, 0.997 | 1.000, 1.020, 1.057, 1.000 |
| **R-hat, window average** | **1.507** | **1.545** |
| **ESS bulk / tail, window average** | **9 / 19** | **8 / 9** |
| gate (R-hat ≤ 1.05, ESS ≥ 100) | **FAILS** | **FAILS** |
| ensemble-median spread | 1.033× | 1.057× |

## What these say

**The hierarchical noise is absorbing only a small part of the kernel error.**
`s` runs 0.84–1.99 on 3S and 1.15–1.42 on 2S — at most a factor of two on the
stated noise, and on 3S ensemble 1 it is *below* one. Absorbing a median 7.96σ
kernel error would take `s` near 8. So the mechanism is present but far from
sufficient.

**The error shows up instead as unexplained misfit.** chi2 per datum against the
true noise is **above 1 on seven of the eight ensembles**, 0.97 to 2.96. The fit
is worse than the true noise permits, which is what remains when `s` does not
stretch far enough.

**The kernel mismatch is comparable to the fault just removed.** 7.96σ median on
3S and 3.23σ on 2S, against the 24–71σ by which the corrected data differ from
the old inversion-generated data. The same *kind* of error, an order smaller but
the same order as the noise many times over, now sits between the data (true
conductivity) and the inversion (pooled median). This is the thing the brief
anticipated, and it is real.

**The window integral survives it.** Pooled window water content is within
0.2% (3S) and 1.5% (2S) of truth, with the truth inside p05–p95 on both columns,
despite 79% and 54% of data points being more than 3σ from what the inversion's
kernel predicts for the true water content. That is consistent with the standing
finding that integrals are far better determined than pointwise values.

**But the four ensembles are not sampling the same posterior, so do not read
0.2% as precision.** R-hat is 1.51 and 1.55 for the window average and ESS is
8–19 against a gate of ≤1.05 and ≥100. With ESS near 9 the Monte Carlo error on
a median is large. **The honest uncertainty on the window average is the spread
of the ensemble medians, 3.3% (3S) and 5.7% (2S)**, not the width of a pooled
posterior that has not converged. On 3S ensemble 1 the pattern is worse than
that: it has the highest misfit (chi2/N = 2.96) and the *lowest* noise factor
(s = 0.84), fitting worse while claiming less noise, which is the signature of a
chain that has not mixed.

## The control is what separates kernel error from the rest

Running now on both columns (`--truekernel`, inversion kernel on the true
conductivity, so (3) is identically zero by construction). Until it returns,
none of the above can be attributed between kernel error, noise and
non-uniqueness.

## A trap worth recording: the temperature swaps between chains

Getting `s` right took three attempts and the first two were wrong in opposite
directions.

1. The nuisance values file has columns **chain, sample, ln s**. Reading the
   first two the wrong way round kept the hot chains and dropped T = 1
   altogether, giving s ≈ 2–4.6.
2. "Chain 1 is the T = 1 chain" is **false**. Every one of the four chains
   visits every one of the four temperatures (1.0, 1.357, 1.842, 2.5); exactly a
   quarter of all rows are at T = 1. That is what parallel tempering does.
3. "The T = 1 chain is the tightest in ln s" is also false — an assertion to
   that effect fired on two of the eight ensembles.

The temperature is recorded per (chain, sample) in the nuisance **misfits** file,
whose columns are chain, sample, misfit, T, accept-rate. The values are now
joined to it on (chain, sample) and only T = 1 rows are kept. The hot chains
tolerate worse fits and inflate `s`: on 3S ensemble 1, T = 1 gives 0.84 where
chains at higher temperatures give 1.9–3.1. Any future reader of these nuisance
histories needs the join.
