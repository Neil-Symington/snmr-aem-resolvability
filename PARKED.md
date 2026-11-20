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

