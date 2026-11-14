# M1.1 item 3 — why the refits stall. Diagnosis from the code and stored output.

## The mechanism

`stage4_gate_v2.jl:179-195`. Each node's constrained solve is **warm-started
from the previous node's model**, sweeping low conductance to high:

```julia
warm = copy(mref)                      # only the first node starts cold
for (i, C) in enumerate(Cgv)
    m, chi, H, nf = occam(aem, warm, m0, R, ndata; fixval = C, mode = :map, ...)
    if nf.nbound > 0 || nf.conres >= 1e-8 || !nf.settled
        ... retry from mref, keep whichever has lower nf.phi ...
    end
    warm = copy(m)                     # the next node inherits this model
end
```

Two consequences, and both are visible in the stored output.

**The profile is path-dependent.** One sweep direction, each node inheriting the
last. Once a node lands in a poor basin every node after it starts there.

**The retry fires on solver flags, never on a bad chi2.** A solve that converges
cleanly to a far worse optimum passes `nbound == 0 && conres < 1e-8 && settled`,
is marked `good`, and is never retried. Nothing in the gate can detect it.

## The evidence: 3S, lambda2 = 10

| node | C (S) | chi2 | phi | conres | nbound | good |
|---|---|---|---|---|---|---|
| 16 | 0.5060 | 32.751 | 89.98 | 0 | 0 | true |
| 17 | 0.6370 | 30.340 | 82.69 | 3e-16 | 0 | true |
| 18 | 0.8019 | 28.845 | 78.91 | 3e-16 | 0 | true |
| 19 | 1.0095 | **28.505** | 86.97 | 0 | 0 | true |
| 20 | 1.2709 | **50.521** | 142.14 | 2e-16 | 0 | **true** |
| 21 | 1.6000 | 221.737 | 368.65 | 0 | 0 | true |

Approaching the minimum from below, chi2 falls by 2.41, 1.50, 0.34 per node. One
node past it, it rises by **22.0**, then by 171. An asymmetry of about seventy
between adjacent 0.1-dex nodes is not a profile likelihood; near its minimum
chi2(C) should be roughly quadratic in log C, and the left-hand curvature
predicts about 28.8 at node 20, not 50.5.

`phi` jumps with it, 86.97 to 142.14, so the solver has converged to a genuinely
worse optimum rather than mis-reporting chi2. And `good = true` throughout, with
`nbound = 0` and `conres` at 2e-16 — every quality flag the gate has says this
node is fine.

## Why the relative band is empty on 2P and 2S

Same cause. If the warm-start chain follows a poor path across the whole sweep,
every node's chi2 sits above `chiref + 3.84` even though better constrained fits
exist. That is what an empty relative band means, and it is why anchoring the
band to N hid the problem instead of fixing it.

## The fix, for item 3

1. **Sweep both directions** and keep the pointwise lower result. This breaks the
   path dependence for one extra sweep per lambda2, the cheapest useful change.
2. **Cold start from `mref` at every node**, in addition to the warm starts, not
   only when a flag trips.
3. **Retry on a chi2 criterion as well as the flags** — a node whose chi2 exceeds
   a local quadratic fit through its neighbours by more than a stated tolerance
   is refit, whatever its flags say.
4. **Flag non-smooth profiles** on that same test and carry the flag into the
   output, so a stalled node can never pass silently again.
5. **Keep the lambda2 sweep** as the evidence: a genuine minimisation should give
   a profile whose shape varies smoothly with lambda2.

## One definitional question to settle first

The retry keeps the lower **`phi`**, the penalised objective. `phi` is what
`mode = :map` minimises, so for the MAP solution at a given lambda2 that is
correct. But the band is defined on **chi2**. If two starts give similar phi and
different chi2, the band wants the lower chi2 and the MAP wants the lower phi,
and they are not the same point. Whether the band should be
"min chi2 subject to the constraint" or "chi2 at the MAP, at fixed lambda2"
needs deciding before the multi-start selection rule is written, because the two
choices select different solves.
