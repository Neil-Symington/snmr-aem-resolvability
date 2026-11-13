# M1.1 — is each gate edge a data limit?

chi2 against forced window conductance, continued past both reported edges.
Threshold is chiref + 3.84. Read only from stored Stage 4 output; no solves.

## 3P  (coarse grid)

chiref 25.374, threshold 29.214, window thickness 16.0 m, sigma search bounds [0.001, 10.0] S/m.
Reported band 0.0296 - 0.3447 S (edge_solves = 0).

| C (S) | sigma (S/m) | chi2 | vs threshold | nbound | solver good |
|---|---|---|---|---|---|
| 0.0160 | 0.00100 | 27.60 | **below** | 7 | **NO** |
| 0.0296 | 0.00185 | 27.46 | **below** | 0 | yes |
| 0.0546 | 0.00341 | 27.46 | **below** | 0 | yes |
| 0.1010 | 0.00631 | 27.77 | **below** | 0 | yes |
| 0.1865 | 0.01166 | 29.46 | above | 0 | yes |
| 0.3447 | 0.02154 | 42.02 | above | 0 | yes |
| 0.6370 | 0.03981 | 455.53 | above | 12 | **NO** |
| 1.1770 | 0.07356 | 3140.31 | above | 20 | **NO** |
| 2.1750 | 0.13594 | 7530.26 | above | 22 | **NO** |
| 4.0190 | 0.25119 | 14275.56 | above | 21 | **NO** |
| 7.4265 | 0.46416 | 20872.14 | above | 15 | **NO** |
| 13.7231 | 0.85770 | 28669.79 | above | 15 | **NO** |
| 25.3583 | 1.58489 | 38349.41 | above | 11 | **NO** |
| 46.8583 | 2.92864 | 83334.03 | above | 15 | **NO** |
| 86.5871 | 5.41170 | 156312.26 | above | 19 | **NO** |
| 160.0000 | 10.00000 | 295355.66 | above | 21 | **NO** |

Nodes below threshold AND solver-good: C = 0.0296 to 0.1010 S.
At the LOWEST search node (0.0160 S) chi2 = 27.60, which is BELOW the threshold -> low edge is **a search limit, not a data limit**.
At the HIGHEST search node (160.0000 S) chi2 = 295355.66, which is above the threshold -> high edge is **a data limit**.

## 3S  (sweep lambda2=10.0)

chiref 27.367, threshold 31.207, window thickness 16.0 m, sigma search bounds [0.001, 10.0] S/m.

| C (S) | sigma (S/m) | chi2 | vs threshold | nbound | solver good |
|---|---|---|---|---|---|
| 0.0160 | 0.00100 | 53.32 | above | 7 | **NO** |
| 0.0201 | 0.00126 | 55.02 | above | 7 | **NO** |
| 0.0254 | 0.00158 | 54.25 | above | 5 | **NO** |
| 0.0319 | 0.00200 | 53.44 | above | 4 | **NO** |
| 0.0402 | 0.00251 | 52.62 | above | 2 | **NO** |
| 0.0506 | 0.00316 | 51.80 | above | 0 | yes |
| 0.0637 | 0.00398 | 50.84 | above | 0 | yes |
| 0.0802 | 0.00501 | 49.80 | above | 0 | yes |
| 0.1010 | 0.00631 | 48.55 | above | 0 | yes |
| 0.1271 | 0.00794 | 47.09 | above | 0 | yes |
| 0.1600 | 0.01000 | 45.36 | above | 0 | yes |
| 0.2014 | 0.01259 | 43.26 | above | 0 | yes |
| 0.2536 | 0.01585 | 40.90 | above | 0 | yes |
| 0.3192 | 0.01995 | 38.25 | above | 0 | yes |
| 0.4019 | 0.02512 | 35.48 | above | 0 | yes |
| 0.5060 | 0.03162 | 32.75 | above | 0 | yes |
| 0.6370 | 0.03981 | 30.34 | **below** | 0 | yes |
| 0.8019 | 0.05012 | 28.85 | **below** | 0 | yes |
| 1.0095 | 0.06310 | 28.51 | **below** | 0 | yes |
| 1.2709 | 0.07943 | 50.52 | above | 0 | yes |
| 1.6000 | 0.10000 | 221.74 | above | 0 | yes |
| 2.0143 | 0.12589 | 346.70 | above | 7 | **NO** |

Nodes below threshold AND solver-good: C = 0.6370 to 1.0095 S.
At the LOWEST search node (0.0160 S) chi2 = 53.32, which is above the threshold -> low edge is **a data limit**.
At the HIGHEST search node (160.0000 S) chi2 = 211350.68, which is above the threshold -> high edge is **a data limit**.

