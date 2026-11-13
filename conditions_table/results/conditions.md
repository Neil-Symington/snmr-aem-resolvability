# Conditions table and water thresholds

Relation, exactly as `stages/stage5_quality.jl:17-18` uses it:

```
sigma_bulk = sigma_w * theta_total^m + clay_term     clay_term = sigma_clay * C
sigma_w    = (sigma_bulk - clay_term) / theta_total^m
TDS        = sigma_w * 6400            [priors].mg_per_sm
```

Limits: theta_total ~ U(0.25, 0.40); m ~ U(1.50, 1.90); clay product prior 1e-04-1.0 S/m.
Clay ceiling = 2.0 x 0.5 x (1 - theta/0.45) S/m, so 0.444 S/m at theta=0.25 and 0.111 at theta=0.40.

Classes as the brief sets them: fresh < 1,500; brackish 1,500-3,000; saline > 3,000 mg/L.
**The code uses a single threshold of 1,000 mg/L** (`[hypothesis.salinity].value`), not these three.

## Column 2P

Window 17.5-33.5 m (thickness 16.0 m). Gate band 0.00199-0.44244 S/m (conductance 0.0319-7.0790 S); lower edge OPEN, upper edge closed. Source: production 1-2 Nov, verdict table.
SNMR water content across 4 ensembles: 0.220 to 0.264 v/v (lowest p05 to highest p95), from cond_C2_typical2PT.jld2.
Truth: sigma_bulk 0.0428 S/m, sigma_w 0.0848 S/m (TDS 543 mg/L), clay term 0.0064 S/m, theta_total 0.303, theta_SMR 0.244.

### Part 1 — what each class requires

| class | supported | clay term required (S/m) | at central theta, m | binds low | binds high |
|---|---|---|---|---|---|
| fresh | **supported** | 0.0000 - 0.4424 | 0.0000 - 0.2778 | clay term >= 0 | gate upper edge |
| brackish | **supported** | 0.0000 - 0.4256 | 0.0000 - 0.2778 | clay term >= 0 | gate upper edge |
| saline | **supported** | 0.0000 - 0.4088 | 0.0000 - 0.2778 | clay term >= 0 | gate upper edge |

Truth is **fresh** (543 mg/L) with clay term 0.0064 S/m; that class is supported and the true clay term is **inside** its required range.

## Column 3P

Window 7.5-23.5 m (thickness 16.0 m). Gate band 0.00185-0.02154 S/m (conductance 0.0296-0.3447 S); lower edge OPEN, upper edge OPEN. Source: production 2 Nov, stage4_synth_3P.jld2.
SNMR water content across 4 ensembles: 0.265 to 0.286 v/v (lowest p05 to highest p95), from cond_C2_typical3PT.jld2.
Truth: sigma_bulk 0.0129 S/m, sigma_w 0.0292 S/m (TDS 187 mg/L), clay term 0.0005 S/m, theta_total 0.320, theta_SMR 0.271.

### Part 1 — what each class requires

| class | supported | clay term required (S/m) | at central theta, m | binds low | binds high |
|---|---|---|---|---|---|
| fresh | **supported** | 0.0000 - 0.0215 | 0.0000 - 0.0215 | clay term >= 0 | gate upper edge |
| brackish | **supported** | 0.0000 - 0.0047 | - | clay term >= 0 | gate upper edge |
| saline | **RULED OUT** | - | - |  |  |

Truth is **fresh** (187 mg/L) with clay term 0.0005 S/m; that class is supported and the true clay term is **inside** its required range.

## Column 2S

Window 17.5-33.5 m (thickness 16.0 m). Gate band 0.00199-0.44906 S/m (conductance 0.0319-7.1850 S); lower edge OPEN, upper edge closed. Source: production 2 Nov, verdict table.
SNMR water content across 4 ensembles: 0.172 to 0.194 v/v (lowest p05 to highest p95), from cond_C2_typical2ST.jld2.
Truth: sigma_bulk 0.0692 S/m, sigma_w 0.2865 S/m (TDS 1834 mg/L), clay term 0.0064 S/m, theta_total 0.303, theta_SMR 0.244.

### Part 1 — what each class requires

| class | supported | clay term required (S/m) | at central theta, m | binds low | binds high |
|---|---|---|---|---|---|
| fresh | **supported** | 0.0000 - 0.4444 | 0.0000 - 0.2778 | clay term >= 0 | clay ceiling (theta) |
| brackish | **supported** | 0.0000 - 0.4322 | 0.0000 - 0.2778 | clay term >= 0 | gate upper edge |
| saline | **supported** | 0.0000 - 0.4154 | 0.0000 - 0.2778 | clay term >= 0 | gate upper edge |

Truth is **brackish** (1834 mg/L) with clay term 0.0064 S/m; that class is supported and the true clay term is **inside** its required range.

## Column 3S

Window 7.5-23.5 m (thickness 16.0 m). Gate band 0.00794-0.07944 S/m (conductance 0.1271-1.2710 S); lower edge OPEN, upper edge OPEN. Source: production 2 Nov, verdict table (SUSPECT: exactly a factor of 10 wide).
SNMR water content across 4 ensembles: 0.262 to 0.282 v/v (lowest p05 to highest p95), from cond_C2_typical3ST.jld2.
Truth: sigma_bulk 0.0510 S/m, sigma_w 0.2901 S/m (TDS 1856 mg/L), clay term 0.0005 S/m, theta_total 0.320, theta_SMR 0.271.

### Part 1 — what each class requires

| class | supported | clay term required (S/m) | at central theta, m | binds low | binds high |
|---|---|---|---|---|---|
| fresh | **supported** | 0.0000 - 0.0794 | 0.0000 - 0.0794 | clay term >= 0 | gate upper edge |
| brackish | **supported** | 0.0000 - 0.0626 | 0.0000 - 0.0448 | clay term >= 0 | gate upper edge |
| saline | **supported** | 0.0000 - 0.0458 | 0.0000 - 0.0101 | clay term >= 0 | gate upper edge |

Truth is **brackish** (1856 mg/L) with clay term 0.0005 S/m; that class is supported and the true clay term is **inside** its required range.

## Part 2 — what one extra measurement would rule out (2P and 3P)

An induction log pins the window conductivity to the truth +/- 10%. A gamma log
pins the clay term to within a factor of 2 of the truth. Neither is used as an
input to anything; they are hypothetical constraints applied to the same grid.

### 2P

| scenario | fresh | brackish | saline |
|---|---|---|---|
| survey alone | 0.0000-0.4424 | 0.0000-0.4256 | 0.0000-0.4088 |
| + induction log | 0.0000-0.0471 | 0.0000-0.0303 | 0.0000-0.0135 |
| + gamma log | 0.0032-0.0128 | 0.0032-0.0128 | 0.0032-0.0128 |
| + both | 0.0032-0.0128 | 0.0032-0.0128 | 0.0032-0.0128 |

### 3P

| scenario | fresh | brackish | saline |
|---|---|---|---|
| survey alone | 0.0000-0.0215 | 0.0000-0.0047 | **RULED OUT** |
| + induction log | 0.0000-0.0142 | **RULED OUT** | **RULED OUT** |
| + gamma log | 0.0003-0.0011 | 0.0003-0.0011 | **RULED OUT** |
| + both | 0.0003-0.0011 | **RULED OUT** | **RULED OUT** |

## Part 3 — water-content thresholds in the representative window

P(theta > x) computed separately in each ensemble, then summarised across them:
all > 0.95 -> "greater than 0.95"; all < 0.05 -> "less than 0.05"; otherwise the range.

| column | ensembles | P(theta > 0.20) | P(theta > 0.25) | P(theta > 0.30) | true theta_SMR | vintage |
|---|---|---|---|---|---|---|
| 2P | 4 | greater than 0.95 | 0.1 - 0.3 | less than 0.05 | 0.244 | 13 Nov (200k) |
| 3P | 4 | greater than 0.95 | greater than 0.95 | less than 0.05 | 0.271 | 22 Oct (100k) |
| 2S | 4 | less than 0.05 | less than 0.05 | less than 0.05 | 0.244 | 22 Oct (100k) |
| 3S | 4 | greater than 0.95 | greater than 0.95 | less than 0.05 | 0.271 | 22 Oct (100k) |


![conditions](fig_conditions.png)

Accepted salinity span, survey alone: 2P **50–38,414 mg/L**; 3P **50–1,852**;
2S 50–38,414; 3S 50–7,056.

## The five questions

**1. On 3P, is saline ruled out, and by which limit? Yes — by the gate's upper
edge.** The accepted conductivity tops out at 0.0215 S/m. Saline water means
sigma_w >= 0.469 S/m, which contributes at least 0.469 x 0.25^1.9 = 0.036 S/m
to the bulk even at the most favourable porosity and exponent. That already
exceeds the entire accepted band before any clay is added, and clay can only
ADD conductance, so nothing rescues it. The formation is simply too resistive
to hold salty water.

**2. On 2P, is any class ruled out from the survey alone? No.** All three are
supported, and their required clay ranges are nearly identical: 0–0.4424 S/m
for fresh, 0–0.4256 for brackish, 0–0.4088 for saline. Those differ by 8%
across the entire salinity axis, which is another way of saying the survey
alone separates nothing on this column.

**3. On 2P, which single measurement rules out the most? Neither rules out
anything.** The induction log narrows the required clay term about ninefold
(0.4424 to 0.0471 S/m) and the gamma log pins it to 0.0032–0.0128 S/m, but all
three classes survive both, and survive the two together. The reason is in the
two priors: porosity U(0.25, 0.40) and exponent U(1.5, 1.9) alone move
theta^m by a factor of **5.6**, which is wider than the whole fresh-to-saline
range. Until porosity and the exponent are pinned, no conductivity measurement
can decide salinity on this column.

**4. Is fresh ruled out anywhere? No** — fresh is supported on all four
columns, as it must be, since nothing stops the clay term going to zero and
the water being arbitrarily dilute within the band.

**5. Are the conditions tight enough to be useful? On 3P yes, on 2P no.** On 2P
every class allows a clay term anywhere from zero to about 0.41–0.44 S/m, which
is essentially the whole physically allowed range, so the table carries no
information. On 3P the same table gives 0–0.0215 S/m for fresh, 0–0.0047 for
brackish and rules saline out outright — a real constraint. The difference is
the gate band: 2P's spans a factor of 222 in conductivity, 3P's a factor of 12.

## Things found while doing this, which you should know

**The code does not use these class thresholds.** `[hypothesis.salinity].value`
is a single **1,000 mg/L**, and the verdict table tests only that. The three
classes here are the brief's, not the code's.

**The measured water content never enters the salinity calculation.**
`stage5_quality.jl:158` draws `theta_total` from `[priors].porosity_range`,
U(0.25, 0.40), declared in the config. The SNMR ensembles are read by Stage 3
for the interval and by Stage 4's floor, but the salinity posterior uses a
prior, not a measurement. So the SNMR half of the survey contributes nothing to
the salinity verdict as the pipeline currently stands.

That is not a simple oversight to fix. SNMR measures **free** water, theta_SMR,
which is a lower bound on theta_total: on 2P the truth is 0.244 against 0.303,
and substituting one for the other would bias sigma_w high by a factor of
(0.303/0.244)^1.7 = 1.44. Using the measurement properly needs the clay
correction, which is the quantity being solved for.

**The 22 Oct SNMR vintage is biased low and is all that exists for three of the
four columns.** On 2P, where both vintages exist, the 22 Oct run (100,001
samples) gives window water content medians 0.184/0.184/0.193/0.183 against a
truth of 0.2445 — **25% low on all four ensembles**, which is systematic. The
13 Nov production run (200,001) gives 0.245/0.244/0.243/0.239. Part 3 uses
13 Nov for 2P and 22 Oct for the rest, so **3P, 2S and 3S water thresholds
should be treated as provisional**. It shows: 2S reports P(theta > 0.20) below
0.05 while its true theta_SMR is 0.244.

**3S's gate band looks like an artefact.** It is 0.1271–1.271 S, exactly a
factor of ten, with both edges open. A band exactly one decade wide with both
edges open is the shape of a search bound, not a measurement. Its 3S numbers
above should not be quoted until that is checked.
