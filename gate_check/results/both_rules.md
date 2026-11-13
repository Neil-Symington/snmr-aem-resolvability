# M1.1 item 5 — accepted band under both rules

`absolute` is the current rule and drives the verdict. `relative` is shown for
comparison only. Bands are the widest contiguous run of nodes that are both
solver-good and below the threshold; no new solves.

| column | truth hash | N | chiref | true window C (S) | absolute band (S) | truth | relative band (S) | truth |
|---|---|---|---|---|---|---|---|---|
| 3P | `1cb14adbbb77` | 42 | 25.37 | 0.2063 | 0.0296 - 0.3447 | **in** | 0.0296 - 0.1010 | **OUT** |
| 3S | `e27701d7ea4b` | 42 | 27.37 | 0.8153 | 0.1600 - 1.0095 | **in** | 0.6370 - 1.0095 | **in** |
| 2P | `6421d57fb629` | 42 | 29.47 | 0.8727 | 0.3192 - 2.0143 | **in** | NaN - NaN | none |
| 2S | `ada389f1b7dd` | 42 | 29.53 | 1.2716 | 0.5060 - 2.0143 | **in** | NaN - NaN | none |

Windows: test one 7.5-23.5 m, test two 17.5-33.5 m.

## What the two rules show

**Under the absolute rule, truth is inside on all four columns.** Under the
relative rule it falls **outside on 3P** — band 0.0296–0.1010 S against a truth
of 0.2063 S, a factor of two above the upper edge.

**On 2P and 2S the relative band is empty.** No node is both solver-good and
within 3.84 of chiref. That should not happen: as the forced conductance sweeps
past the unconstrained optimum's own window conductance, one constrained fit
ought to come arbitrarily close to chiref. An empty relative band means the
constrained refits never get within 3.84 of the unconstrained solve anywhere on
the grid.

That is direct evidence for item 3. It is the same pathology the field work
recorded, seen from the other side: there the *unconstrained* solve stalled
*above* the constrained sweep; here the *constrained* sweep stalls above the
unconstrained solve. Either way the two are not being minimised to the same
standard, and the gap is larger than the 3.84 the band is built on.

## Item 6 — where the +3.84 comes from

**The 3.84 itself is stated.** `synthetics/PLAN_CHAINED_HYPOTHESES.md:237` lists
the reading as

> **admissible range** — where Δχ² ≤ 1.00 / 3.84 / 6.63 | the data's own bound
> on this unit's conductivity, one parameter, prior-free

which are the 68%, 95% and 99% points of chi-squared on one degree of freedom.
So 3.84 is the one-parameter 95% increment, and the plan specifies it as
**Δχ², relative to the minimum**.

**Anchoring it to N instead is also stated, and the reason is given — in the
field report, not in the synthetic code.**
`ek_realdata/REPORT_v1_firstfive.md:559-563`:

> At some sites the unconstrained Gauss–Newton stalled above the constrained
> sweep (KRP_t5_s2: χ² 30.5 unconstrained against 19.9 constrained), which
> cannot be true of a real optimum and means the single-start unconstrained
> solve is not the minimum. The band is referenced to `ndata + 3.84`, not to
> that minimum, so it is unaffected; the tables quote the profile minimum
> rather than the unconstrained value.

and, on the direction of the error:

> an under-converged profile overstates Δχ², narrows the band and makes
> certification too *easy*.

So the rationale is a **deliberate robustness guard against an unreliable
chiref**, not a statistical argument for an absolute criterion. The 3.84 keeps
its profile-likelihood meaning while losing the reference point that gives it
that meaning.

`config.toml:385` records the choice ("referenced to the TARGET misfit, not to
the best fit found") but gives no reason; the reason is only in the field
report.

**The consequence for this work.** The guard was chosen because the refits could
not be trusted. Item 3 is an attempt to make them trustworthy. If it succeeds,
the reason for anchoring to N goes away and the relative rule becomes available;
if it fails, the absolute rule is still the right call. So items 3 and 5 should
be read together, and the rule should not be changed before item 3 reports.
