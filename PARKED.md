# Parked

Out of scope for the production re-run. One line each; no work until unparked.

- **Occam comparison beyond what exists.** 2P and 2S are done; the verdict matched the sampler and the located window agreed to a quarter of a metre, while the interval did not.
- **Noise mis-specification test.** Whether chiref near 0.60N on test one reflects overstated noise or effective degrees of freedom; see the M1.1 report.
- **Relative gate threshold.** Using a fraction of chiref instead of a fixed chi2 + 3.84.
- **Other datasets.** Nothing outside the four synthetic columns.
- **SNMR decay times.** T2* and the dead-time correction.
- **Clay percentages from the clay term.** The clay term stays as a conductivity; converting it to a clay fraction is a later, illustrative step.
- **A more conductive unconfined column.** Test one is resistive; a conductive unconfined case would separate survey geometry from ground.
- **Field sites.** East Kimberley, Upper Darling and the retired field workflow.
- **Convergence checks in the other pipelines.** `ek_prerun`, `ek_realdata`, `udf_realdata` and `ssc_prerun` share the defective `ess_basic` and gate it at ESS >= 100; correcting it changes past pass/fail outcomes.
- **Clay limits from Upper Darling borehole logs.** Would replace the stated clay ceiling with measured bounds.
- **Per-draw SNMR kernel (condition C4).** The kernel is conductivity-sensitive; on earlier work the kernel choice moved P(2 aquifers) from 0.09 to 0.45, so this is a real effect, not a refinement.
- **Selectable window rules.** `[hypothesis.interval].select` with `robust_water` (the current rule, default), `mean_theta`, `shallowest` and `deepest`, applied across the whole sweep, plus a test for each and the rule recorded in the output. Specified in full at DECISIONS.md commit `9c2a5a2`.
- **Interval labels and the second-interval layout.** `[hypothesis.interval].label` defaulting to `"main"`, Stage 3-6 outputs under `[output].dir/<name>/<label>/` with Stage 0-2 shared, a loud `--force`-gated refusal to overwrite Stage 3-6 outputs in `run_site.jl`, the label and depth range and rule recorded in every output and the manifest, and the site report listing every label found. Needed before any second interval is run, because every Stage 3-6 filename is keyed on the site name alone and the overwrite is silent.
- **`report_level` and `secondary_level` from config.** Literals today at `stage3_hypothesis.jl:323` and `:356`.
- **Depth-dependent stated limits.** The per-cell machinery exists (`setcellbounds!`, `boundkind`) but no stored run uses it; every release store carries one global `mbounds` pair.
- **The JSON data file for the interactive page.** Schema supplied 17 Nov, recorded at DECISIONS.md commit `9c2a5a2`. Needs the smoothest-passing model per node (not stored; about 1.5 h of new solves), a stated mapping from the stores' `located`/`open` edge states to `data|limit|search`, and the classes work.
- **Ten-draw coverage test.** About 2.1 h a draw: 21 h serial, roughly 7 h three-up, under STATED.
- **Noise-consistency follow-up on the SNMR misfit.** The check itself was already run before this was parked and is recorded in DECISIONS.md under 18 Nov: the true water content forwarded through the true kernel gives **chi2/N = 0.806** against the noisy data on all four columns, consistent with 1 at -0.67 sigma, so sigma-hat does not understate the added noise and the inversions' chi2/N of 1.0-3.0 is the sampler failing to reach an achievable fit. What is parked is the follow-up Neil asked for conditionally on a value near 2 -- documenting how sigma-hat is derived and how the per-moment noise shape and the typical/poor tiers are drawn. Not needed: the value is near 1.

## Release tidy: the role label "saline_control" on a brackish column

Parked 21 Nov 2026, to rename after the paper.

`src/config.jl:210` accepts `[data].role` of `primary|saline_control` only, so
the brackish columns -- 3S, 2S and now 2Sv2 -- all carry `saline_control`.
That label is wrong twice over. The partner columns are BRACKISH, not saline;
and they are not controls of the fresh columns, they are the other arm of a
paired comparison. `role = "variant"` was tried on 2Sv2 and rejected by the
validator.

The rename touches the validator's allowed set, the five site configs that use
it, and anything that reads `role` downstream (Stage 0 seals it, so every
stage0 artefact carries the string). It is cosmetic and it is not worth doing
while runs are in flight.

## Salinity class limits live in the figure code, not the config

Parked 21 Nov 2026, to do after the rebuild and the gate runs land.

`tools/make_figures.py:63` holds `CLASS_LINES = [(1500.0, "1,500 mg/L"),
(3000.0, "3,000 mg/L")]`. Those boundaries appear on the misfit and ceiling
figures and on the report pages as though they were study parameters, but they
are in the plotting code and in no config. `config.toml:322-324`
`[hypothesis.salinity]` carries a single `value = 1000.0` instead.

To do: replace that single value with a list of class boundaries under
`[hypothesis.salinity]`, 1500 and 3000 mg/L for now, and have the figures and
the report pages read them from the config. No result changes -- nothing in
Stages 0-6 reads CLASS_LINES, and the gate's own threshold is unrelated.

## Apply the SNMR water content as a lower limit on porosity in Stage 5

Parked 22 Nov 2026, to do after the rebuild and the gate runs land.

`stages/stage5_quality.jl:98` reads `plo, phi_ = Float64.(pr["porosity_range"])`
and draws theta_total ~ U(0.25, 0.40) from `config.toml:213`. Nothing in that
file refers to theta_SMR, so the salinity ceilings are computed over porosities
the SNMR data already exclude: the pooled SNMR p05 over the window is 0.2726 on
3P, 0.2634 on 2P and 0.2602 on 3S, all above the 0.25 floor.

To do: set the lower limit to `max(config value, pooled SNMR p05 over the
window)`. The most cautious ceiling falls 5-13%, and the direction is right --
if the rock holds more water than 0.25, less of the window's conductance need
be attributed to salt. The least cautious ceiling is unchanged, being computed
at phi = 0.40.

## State the temperature convention beside mg_per_sm

Same date, same batch. `config.toml:190` gives `mg_per_sm = 6400.0` with no
reference temperature anywhere in the repository, and the conversion is a bare
multiplication (`stage5_quality.jl:95`, `make_figures.py:324`). 6,400 mg/L per
S/m is the usual 25 C figure, but the code does not say so and a reviewer asked.
Add the convention as a comment beside the key.

## Stage 3 still reads zsmr as layer centres

Parked 22 Nov 2026. Do after the chain lands, in one go.

`stages/stage3_hypothesis.jl:91-96` builds cells from `zsmr` as layer CENTRES,
edges at the midpoints, first top 0, last bottom extrapolated by half the final
spacing. The forward model treats `zsmr` as layer BOTTOMS with the first layer
from 0 to 1.0 m, which is the convention now in `src/window.jl` and used by
`tools/stage2_report.jl`.

The difference on the truth window mean is -0.47% on test one and +1.22% on
test two. That is enough to move p_pass and the tie-break window, so it is a
RESULT CHANGE, not a tidy-up, and it must not happen while columns are being
run through the stage as it stands.

To do, in one pass over every column: switch Stage 3 to `Window.window_mean`,
re-run Stage 3 on 3P, 2P, 3S, 2S, 2Pv2 and 2Sv2, and report the passing
intervals and the tie-break window BEFORE and AFTER for each. Stages 4 and 5
have already run on the current tie-break windows, so any column whose window
moves needs its downstream outputs re-examined before the change is adopted.

## Active-set treatment of the parameter bounds in occam_core.jl

Parked 23 Nov 2026, BEHIND the sweep-splitting validation. It changes solves,
so it needs the full check against `gate-of-record-20261119`: every node's
min chi2 to 1e-6 and every edge unchanged.

THE PROBLEM. The Gauss-Newton step is computed with no knowledge of the
parameter bounds; `clampm!` then pushes the iterate back inside them, which
breaks the conductance constraint the step had just been projected onto.
`projcond!` (`src/occam_core.jl:187-197`) re-projects exactly -- scaling every
window cell by a common factor is closed-form -- but it cannot fix a cell
pinned at the 2 S/m ceiling, and it gives up after 12 attempts.

TWO SYMPTOMS, BOTH MEASURED.

  * Solves run ~90-150 accepted steps. The traced solve at C = 0.8997 S,
    lambda2 0.1, took 118 steps and 968 forward calls, 88.5% of it Jacobian
    evaluation, and chi2 stopped moving by more than 0.01 at step 88.
  * At C = 8.9975 S on 2Pv2 the constraint is LOST: conres 4.2e-3 at
    lambda2 0.01 and 1.5e-4 at 0.1, even from a seed rescaled to meet the
    constraint exactly and even at a step cap of 400. A model meeting that
    conductance within the bounds demonstrably exists -- the seed is one -- so
    this is the optimiser losing the constraint, not infeasibility. It is why
    2Pv2's upper edge reads OPEN.

THE CANDIDATE. Hold the clamped cells fixed, and solve and project over the
free cells only, which is what the seed-builder in the node-12 check already
does. That should both keep the constraint and cut the step count.

NOT a speed-up in the sense the sweep-splitting is: splitting changes no
number, this changes every solve. Order: validate splitting first, then this,
then the stopping rule and Jacobian reuse.

## The edge finder reports OPEN when it should report UNRESOLVED

Parked 23 Nov 2026. Fix WITH the active-set work, which is what makes the
distinction decidable.

`EdgeFinder.bracket` looks for the lowest REJECTED node above the highest
admissible one. A node that could not be solved at all -- the constrained
optimiser never reached the constraint, so chi2 is Inf -- is not rejected, so
it is not a bracket, and the edge is reported OPEN. "Open" means the data
admit arbitrarily large conductance. It is the wrong word for "the solver
stopped working here", and the two are reported identically.

MEASURED, on 2Pv2. The profile's grid runs 5.0596 S (admissible, chi2 27.79),
then 8.9975 S (Inf), 16.0 S (Inf), 28.45 S (Inf). Both edges came back OPEN.
A hand probe at 6.50 S, lambda2 0.01, solved cleanly -- conres 1.4e-16, no cell
on the 2 S/m ceiling -- at chi2 49.43, which is REJECTED against 45.84. So a
real upper edge exists between 5.0596 and 6.50 S and the column was reported
unconstrained when it is not. 2Sv2 has the same shape (five unsolved nodes) and
is under the same suspicion.

WHAT IT SHOULD DO. When the node beyond the last admissible one is unsolved,
probe INWARD from it -- halving toward the admissible node -- until it either
brackets or runs out of room, and report `unresolved` if it still cannot. Never
`open` on the strength of an unsolved node.

Note the lower end is a different case and is correct as it stands: the lowest
grid node is 0.0160 S on both columns, which is exactly the stated floor
0.001 S/m times the 16 m window, so nothing below it exists to solve. That end
is set by the stated limit, and saying so is accurate.

