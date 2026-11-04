# Choosing the Stage 0 priors

Everything here is a *declared* choice. None of it is measured by the survey, and all
of it changes the answer, which is why it lives in `config.toml` where it can be seen
and argued with rather than buried in a script.

> **Citations are pending and are marked `[REF]`.** Every range below is the one the
> runs actually used, and the reasoning for each is given in full; what is missing is
> the supporting literature, which has deliberately not been invented. Treat a `[REF]`
> as "this choice is argued here but not yet cited", not as a gap in the reasoning.

## The rock-physics model

```
sigma_bulk = sigma_w * theta_total^m  +  sigma_clay * C
```

Bulk conductivity is pore water conducting through the connected pore space plus a
surface term on the clay. Both terms are S/m. `theta_total` is total water content,
not the SNMR water content — see below. TDS is converted at `[priors].mg_per_sm`,
6,400 mg/L per S/m, the usual mid-range figure for NaCl-dominated groundwater at
25 °C `[REF]`.

## Cementation exponent `m`

Keyed to lithology rather than given as a number, because a single value is a
measurement nobody made:

| keyword | range | for |
|---|---|---|
| `clean_sand` | U(1.3, 1.6) | well-sorted quartz sand, negligible clay |
| `clayey_sand` | U(1.5, 1.9) | sand with a sampled clay term — **the synthetic columns** |
| `consolidated_sandstone` | U(1.8, 2.1) | cemented, lower porosity |

`[REF]` for each band.

**Why the synthetics use `clayey_sand` and not `clean_sand`**, which an earlier draft
specified: these columns carry 6–12% clay with a sampled clay term, so the project's
own lithology rule puts them in the clayey band. The truth's m = 1.7 sits inside
U(1.5, 1.9) and *outside* U(1.3, 1.6) — clean sand would have excluded the true value
from the prior's support, which is a prior that cannot contain the answer.

Nothing in this project ever infers `m`. Every synthetic script written before
November 2026 used the truth's 1.7 directly, which is the truth firewall leaking in
through the rock physics; that is the bug this keyword exists to prevent.

## Total porosity

`porosity_range = [0.25, 0.40]`, uniform, for the synthetic sands. A **formation
range**, declared, not a measurement `[REF]`.

This is not the SNMR water content. At the instrument's 10.5 ms dead time the SNMR
sees free water plus roughly 0.63–0.72 of the capillary water — a fraction predicted
from T2, not fitted — so `theta_SMR` is neither porosity nor specific yield. The
porosity prior and the clay term are what reconstruct `theta_total` from it.

## The clay term

One log-uniform prior on the **product** `sigma_clay * C`, over 1e-4 to 1 S/m.

The gate cannot separate clay conductance from pore-water conductance, so the prior
goes on the product that actually enters `sigma_bulk` — the only combination the data
could ever constrain. Putting it anywhere else invents structure the measurement
cannot see.

**What this replaced, and why.** A logit-normal clay fraction (median 0.15, logit-SD
1.2, 6 m correlation) times a log-uniform `sigma_clay` over 0.10–2.0 S/m put **98.1%
of its mass on column 3P in physically impossible territory**: that column's entire
measured conductance is 0.0129 S/m, while a clay fraction of 0.15 times a median
`sigma_clay` of 0.45 S/m gives a clay term of 0.067 S/m — five times everything
observed. Stage 5's estimate then rested on the surviving 1.9% of draws. The product
prior covers the same physics without concentrating mass above what was measured.

**Acceptance rate is a result, not a diagnostic.** Draws whose clay term alone exceeds
the measured conductance give a negative `sigma_w`: the data are consistent with
arbitrarily fresh water. Those draws are counted, never clamped, and the acceptance
rate is a row in the verdict table. Below 10% the Stage 5 estimate is flagged FRAGILE,
because it is then a statement about the prior's tail rather than about the ground.

## Ceilings, for the salinity floor

`[ceilings]` is a different kind of object: not a prior but a set of stated maxima used
to compute a *floor* on salinity. Cap everything that could be explaining the measured
conductance instead of salt — clay at `sigma_clay_max` times `C_max(theta)`, porosity
at `porosity_max`, `m` at its prior's lower bound — and whatever conductance is left
must be water.

Three things to know before using it:

* **An open lower band edge gives a floor of zero by construction.** An open edge sits
  at the solver's `sigma_min`, not at the data, so there is no measured minimum
  conductance to subtract from. All four synthetic columns are open at the lower edge.
* **At `sigma_clay_max` = 2 S/m the ceiling is not a constraint on this ground.** It
  lands at 0.40–0.59 S/m against window-mean conductivities of 0.0095–0.155 S/m — four
  to forty times the whole measurement. A ceiling above everything observed constrains
  nothing, which is the same failure the old clay prior had.
* **Resistive ground cannot be floored at all.** At `porosity_max` = 0.45 and m = 1.50
  a 1,000 mg/L floor needs `sigma_bulk` ≥ 0.047 S/m, so ground more resistive than that
  will return zero whatever the clay ceiling is.

`porosity_max` = 0.45 is **not** the Stage 0 cap of 0.40. It is a looser physical
ceiling, and looser is the conservative direction for a floor: a larger porosity
divides the leftover conductance by more and so lowers `sigma_w`.

`m_for_floor = "lower"` because `sigma_w = X / porosity_max^m` with `porosity_max < 1`,
so `phi^m` shrinks as m grows and `sigma_w` *rises* with m. The lowest `sigma_w` — the
actual floor — is at the m prior's lower bound, not its upper one.
