# Review brief — snmr-aem-resolvability

For a reviewer on another machine. Written 24 Nov 2026.

## What this is

A joint AEM + surface-NMR resolvability study: can airborne EM and SNMR
together say whether a groundwater interval is fresh, when neither can alone?
Synthetic columns with known truth, inverted through the same workflow a field
survey would use, then a **gate** that asks which window conductances are
admissible against the AEM data.

## Get it

    git clone https://github.com/Neil-Symington/snmr-aem-resolvability.git
    cd snmr-aem-resolvability
    git checkout m1-production-rerun          # NOT main; see "branches"
    julia --project=. -e 'using Pkg; Pkg.instantiate()'

Julia **1.12.6**. Everything below was measured on a 6-core box, 31 GB.

### Branches

- `m1-production-rerun` — the study, the paper's numbers, the new gate. Review this.
- `main` — has field-survey work (`fieldrun/`) that `m1-production-rerun` does
  NOT, merged there on 24 Nov. The two have diverged; merging them is in progress.

## Run the tests first

    julia --project=. test/runtests.jl        # ~2-3 min, 150+ assertions

They should all pass. Several encode bugs that cost real time and must not
return — read their headers; they are the most informative files here.

## The thing to review hardest: the gate

`gate_check/uniform_gate.jl`, `uniform_bisect.jl`, `uniform_check.jl`.

**The formulation.** For a trial conductance C the whole window is held at
log10(C/L) and REMOVED from the unknowns; only cells outside it are inverted.
No constraint row, no projection, so the achieved conductance is exactly C.
This REPLACED a constrained formulation that could not hold its constraint at
high conductance and reported "open" edges that were really solver failures.

**What to check:**

1. Is the uniform window the right question to ask of the data? It removes
   shape freedom inside the window, so its chi2 is an UPPER bound and its edges
   are conservative. Is that the right conservatism for the claim being made?
2. The verdict rules: a solve stopped by the stall rule counts only if its chi2
   is more than 2.0 from the threshold; inside that band three starts must
   agree to 0.5; a solve that stopped on the forward-call budget is an upper
   bound that can prove admissible but never rejected. Are those sound?
3. Bracket validity: an edge is only meaningful if the admissible end is below
   its threshold and the rejected end above it. `uniform_check.jl` asserts this
   because seeding one bracket for two different thresholds produced a wrong
   relative edge that the first check passed as complete.

## Reproduce a number

The cheapest meaningful check is one gate value on 3P. Its admissible range is
0.016-0.36 S (absolute rule, threshold chi2 <= 45.84):

    julia --project=. gate_check/uniform_gate.jl --col 3P \
      --cs 0.33,0.40 --out /tmp/review --shard r

Expect **0.33 S admissible, 0.40 S rejected**. About 10-20 minutes.

## Known limits — please check these rather than rediscover them

- **Stage 1 does not meet its convergence gate** on any column: R-hat 1.07-1.29
  against a 1.05 target, ESS 13-72 against 100. The gate does not read Stage 1
  draws, so edges are unaffected, but the posterior itself is under-sampled.
- **Three conventions for the SNMR depth grid** existed in the repository. One
  is now shared (`src/window.jl`); `stages/stage3_hypothesis.jl` still uses its
  own, which differs by -0.47% to +1.22% on the window mean.
- **3S's gate models** were lost when a process was stopped; they are being
  regenerated. Its edges are sound.
- **`PARKED.md`** lists seven known fixes not yet made, with evidence.
- **`RUN_REPORT.md`** lists every result lost or superseded, and why.

## What I am least sure of

1. Whether the uniform window is too conservative to be the paper's gate.
2. The solver: 88.5% of each solve is Jacobian evaluation, solves run 90-180
   steps, and the step ignores the parameter bounds so `clampm!` can undo the
   constraint. An active-set treatment is parked, not done.
3. Stage 1's under-sampling, above.

## Report back

Anything that is wrong, anything the code does that the documentation does not
say, and whether the gate's formulation and verdict rules are defensible.
