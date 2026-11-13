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
