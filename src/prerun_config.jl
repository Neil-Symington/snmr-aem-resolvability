# prerun_config.jl — every setting of the EFTF survey pre-run, in one place.
#
# Source: "Survey pre-run — parameter specification (26 Sep 2026)", henceforth
# THE SPEC, section numbers below refer to it. Constraints marked RUNBOOK come
# from section4/RUNBOOK.md and override a naive reading of the spec where the
# environment cannot deliver what the spec assumes.
#
# Nothing here is tuned to a site (spec §0). Departures are listed at the bottom
# and are reproduced verbatim in RESULTS_PRERUN.md.

module PrerunConfig

using HiQGA.transD_GP

# ── §0 Scope ────────────────────────────────────────────────────────────────
# The spec says "all 280 compatible FID soundings". The EK archive holds 333
# compatible FID soundings (260 sites, 432 files incl. EMI pairs) and no filter
# on pulse moments, misfit, noise, depth or pulse length reproduces 280.
# Operator decision (26 Sep 2026): run all 333; the analysis table carries every
# covariate so a 280-subset can be taken afterwards without re-running.
const RUN_ALL           = true

# ── §1 Model and priors ─────────────────────────────────────────────────────
const NZ                = 60           # cells
const ZSTART            = 1.0          # m, top cell
const EXTENDFRAC        = 1.028        # → zbnd[end] = 147.4 m ≈ "bottom ≈150 m"
const DZ                = 1.0
const FBOUNDS           = [0.001 0.45] # water-content bounds, same for Occam
const NMIN, NMAX        = 2, 40        # k ~ U{2,40}
const LAMBDA_GP         = [3.0]        # GP correlation length, in depth cells
const DELTA             = 1e-4         # GP additive variance  (see DEPARTURES)
const GP_KERNEL         = transD_GP.GP.OrstUhn()

# ── §2 Noise model ──────────────────────────────────────────────────────────
# σᵢ² = s·σ̂ᵢ², σ̂ⱼ the per-pulse-moment shape, s the sampled level.
#
# The spec's σ̂ⱼ hierarchy is  stack scatter → GMR-reported → flat+FLAGGED.
# The EK .mat files carry neither raw stacks nor GMR per-q errors, so a literal
# reading lands on "flat" for all 333. RUNBOOK §4 (open item 2) records that the
# real per-moment σ̂ *is* recoverable from the FID residual, which the usual
# workflow throws away when it collapses each record to a single V₀ — so that is
# what is used here, as an explicitly-named fourth branch, never silently.
# Measured: the per-q shape has median CV 0.61 across 12 soundings, i.e. real
# structure, not a flat shape dressed up.
const SIGMA_HAT_BRANCH  = :fid_residual_per_q
const S_PRIOR_CENTER    = 1.0          # s₀ — no WP1 tail σ, so "trust σ̂ⱼ level"
const S_PRIOR_TAU       = 0.5          # log-normal width
const S_SDEV            = [0.3]        # Gaussian RW in ln s
const S_BOUNDS          = [log(1e-2) log(1e2)]   # wide; the prior does the shaping
const S_INIT            = 1.0

# ── §3 Sampler ──────────────────────────────────────────────────────────────
const NCHAINS           = 4            # PT chains per ensemble
const TMAX              = 2.5          # log-spaced ladder (HiQGA builds it so)
const NSAMPLES          = 150_001      # spec: 150 000 (100 000 if compute-limited)
const BURNINFRAC        = 0.5
const SAVE_FREQ         = 10           # "keep every 10th post-burn-in sample"
const NENSEMBLES        = 4            # RUNBOOK §5: m ≥ 4 independent ensembles
const SEED_BASE         = 20260926     # per-sounding stream = SEED_BASE + hash(site)

# ── §4 Convergence gates ────────────────────────────────────────────────────
const GATE_RHAT         = 1.05
const GATE_DEPTHS       = [5.0, 15.0, 30.0, 60.0, 100.0]   # m
const GATE_SPLITHALF_DW = 0.03         # v/v agreement between split halves
const GATE_MIN_ESS      = 100          # RUNBOOK §5: R̂ alone must not certify
const TREND_SIGMA       = 2.0          # misfit-trace stationarity tolerance

# ── §5 Occam companion ──────────────────────────────────────────────────────
const OCCAM_REGTYPE     = :R2
const OCCAM_W0          = 0.02         # uniform start / reference
const OCCAM_LAM2MIN     = -1.0         # log10 μ² ∈ [−1, 6] — already the defaults
const OCCAM_LAM2MAX     = 6.0
const OCCAM_UNCERTAINTY = :laplace     # linearised Cₘ; bootstrap skipped (spec §5)
const DOI_REFS          = (0.02, 0.10) # two-reference-model DOI, Christiansen & Auken 2012

# Reference-anchoring weight. NOT in the spec, and without it the spec's §5 is
# self-contradictory: SMRPInversion applies `wref` ONLY when β² > 0 (its own
# docstring says so), and more fundamentally a difference operator annihilates a
# constant, so R·(m − m_ref) = R·m for any uniform reference. With β² = 0 the
# two DOI reference inversions are mathematically identical — measured max
# separation 0.0000 for 330 of 333 soundings, giving a meaningless "resolved to
# the base" DOI everywhere.
#
# Chosen by sweep on 104201_trans_st10 (β² ∈ {0, 0.01, 0.1, 1, 10}): 0.01 is the
# smallest value that makes the references bite — separation 0.972, DOI 111.3 m
# instead of 147.4 — at a fit cost of 2–5% in χ²/N (47.13 → 47.75/49.55). Larger
# values buy no further separation and degrade the fit steeply (β²=1 → χ²/N 122).
const OCCAM_BETA2       = 0.01

# ── §6 Derived-output settings ──────────────────────────────────────────────
const DOI_WIDTHFRAC     = 0.6          # primary; 0.5 and 0.7 for sensitivity
const DOI_WIDTHFRAC_ALT = (0.5, 0.7)
const DETECT_THETA      = 0.10         # v/v
const DETECT_MINTHICK   = 5.0          # m
const DETECT_OCCAM_CHI2 = 1.5          # Occam detection also requires χ²/N ≤ this
const KL_NBINS          = 40           # histogram bins for ∫KL dz

"Depth grid: identical for every sounding so posteriors share a common grid (§6)."
function grid()
    transD_GP.setupz(ZSTART, EXTENDFRAC; dz = DZ, n = NZ, showplot = false)
end

# ── Departures from the spec, recorded so the production run can differ ──────
const DEPARTURES = [
 ("§0 scope",
  "Spec says 280 compatible FID soundings; the archive holds 333 and no filter " *
  "reproduces 280. Running all 333."),
 ("§1 δ",
  "Spec gives GP additive variance δ = 1e-4 as 'the SMRPI value', but " *
  "synthetics/common.jl — the repo's SMRPI convention file — uses 0.01. The spec " *
  "value is used; the discrepancy is unresolved and matters (two orders of " *
  "magnitude in the nugget changes GP smoothness)."),
 ("§2 σ̂ⱼ branch",
  "No raw stacks and no GMR-reported per-q errors exist in the EK .mat files, so " *
  "the spec's hierarchy falls through to 'flat and FLAGGED'. Instead the per-q σ̂ " *
  "is recovered from the FID residual (RUNBOOK §4 open item 2), recorded as " *
  "branch :fid_residual_per_q. Measured per-q shape CV = 0.61, so the shape is real."),
 ("§2 s proposal fraction",
  "Spec asks for nuisance moves on ~15% of iterations. HiQGA's stat+nuisance " *
  "driver (MCMC_Driver.jl:625) performs exactly one nuisance move and one model " *
  "move every iteration — 50%, with no knob. Not adjustable; recorded."),
 ("§3/§4 R̂ across chains",
  "Spec gates R̂ < 1.05 across the 4 PT chains. Under PT the chains target " *
  "p^(1/T), and exactly one chain holds T = 1 at a time (RUNBOOK §1.5, §5), so a " *
  "within-ensemble R̂ measures ladder spacing, not convergence. Instead 4 " *
  "independent ensembles are run per sounding with distinct seeds and prefixes, " *
  "and rank-normalised split-R̂ (Vehtari et al. 2021) is computed across their " *
  "T=1 streams. Cost: 4× the nominal sampler budget."),
 ("§4 ESS gate added",
  "RUNBOOK §5 records that ESS binds long before R̂ (IAT of log s: median 8, max " *
  "51 in stored samples). An ESS ≥ $(GATE_MIN_ESS) condition is gated alongside R̂; " *
  "bulk- and tail-ESS are reported. Not in the spec, but the spec's gate is not " *
  "meaningful without it."),
 ("§5 two-reference DOI needs β² > 0",
  "Spec §5 pairs an R2 regulariser with the Christiansen & Auken two-reference " *
  "DOI. These are incompatible as written: a difference operator annihilates a " *
  "constant reference, and SMRPInversion applies wref only when β² > 0 (default " *
  "0). Measured: max separation 0.0000 for 330/333 soundings, DOI = grid bottom " *
  "for 99%. β² = 0.01 is used, chosen by sweep as the smallest value that makes " *
  "the references influential at negligible fit cost."),
 ("§5 χ² target unreachable under σ̂",
  "Occam targets χ² = N under σ̂ⱼ² at s₀ = 1, but σ̂ understates the noise ~4.5×, " *
  "so the target cannot be met: median χ²/N = 20.28 and 306/333 soundings end " *
  "with log₁₀λ² pinned at the 6.0 ceiling. Every model in the spec-literal Occam " *
  "run is therefore the smoothest the bounds allow, not the smoothest that fits. " *
  "Independent confirmation of the noise level: Occam √20.28 = 4.50 vs trans-D " *
  "hierarchical s (p50 ≈ 20) = 4.47, agreeing to better than 1%."),
 ("§6 missing covariates",
  "TMI gradient, distance to infrastructure, month, stack count and " *
  "coherent-residue ratio do not exist for EK — the archive has no site " *
  "coordinates, no dates and no raw stacks. Ref_loop IS available " *
  "(EK_SNMR_data_spreadsheet.csv). §7's 's/σ̂ tail listed with TMI gradient and " *
  "AEM distance alongside' cannot be produced as written."),
]

end # module
