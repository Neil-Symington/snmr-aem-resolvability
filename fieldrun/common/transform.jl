# transform.jl — the petrophysical transform, sampled.
#
# SMR gives θ_SMR. AEM gives σ_bulk. Neither gives σ_w, the pore-water
# conductivity the hydrogeological hypotheses are actually about. Getting there
# needs the forward transform inverted:
#
#   θ_SMR  = θ_total·(1 − clay)^p        p = 1.5
#   σ_bulk = σ_w·θ_total^m + σ_clay·clay  σ_clay = 0.5 S/m
#
#   =>  θ_total = θ_SMR / (1 − clay)^p
#       σ_w     = (σ_bulk − σ_clay·clay) / θ_total^m
#
# clay and m are unknown, so they are DRAWN from priors per sample and the
# spread they induce is part of the answer. That is the whole point: a
# conductive layer is conductive because of salt or because of clay, and
# σ_bulk alone cannot say which. Panel f of Figure 3 is this ambiguity.
#
# ADMISSIBILITY IS ENFORCED OVER A WINDOW, NOT THE WHOLE PROFILE. Requiring
# every cell to be simultaneously admissible gives a zero acceptance rate: the
# resistive basement (σ_bulk ≈ 0.017 S/m) forces clay ≤ 0.034 there, and a
# correlated prior centred on 0.15 essentially never delivers that at every
# deep cell at once. It is also the wrong question — a hypothesis about 0–60 m
# should not be rejected because of a cell at 150 m. Each predicate therefore
# conditions on its OWN depth window, and the acceptance rate for that window is
# reported, because a low one means the prior and σ_bulk disagree there.
#
# Not every (clay, m) is possible given the data. Two physical
# constraints bite:
#   * θ_total = θ_SMR/(1−clay)^p must not exceed plausible porosity (0.5)
#   * σ_w must be ≥ 0: if σ_clay·clay already exceeds σ_bulk, the clay alone
#     over-explains the conductivity and that clay value is ruled out
# Draws violating either are rejected, which makes this rejection sampling from
# the prior restricted to the admissible set — i.e. the posterior of (clay, m)
# under a uniform likelihood on what is physically reachable. Rejection rates
# are returned, because a high one means the prior and the data disagree.
#
# PRIOR SUBSTITUTION, stated rather than hidden: §5 asks for "clay ~ Beta
# centred 0.15, broad". Neither Distributions nor SpecialFunctions is in this
# environment, and a Beta quantile is needed for the correlated (copula) draw.
# A logit-normal centred on 0.15 is used instead: same support, same centre,
# broad in the same sense, and it needs only `randn`. Its P(clay > 0.4) is
# ~0.07, so S-F still has to be carried by the data.

module Transform

using Random, Statistics, LinearAlgebra

const P_CLAY  = 1.5
const σ_CLAY  = 0.5      # S/m
const φ_MAX   = 0.5      # porosity ceiling for admissibility
const σW_FLOOR= 1e-4     # S/m
# σ_bulk is an AEM INVERSION RESULT, not a measurement: its p10–p90 spans ~3x
# the median on these runs. Treating it as exact in the admissibility test was
# the bug. When σ_clay·clay exceeds σ_bulk the clay alone over-explains the
# conductivity — which is the FRESH-WATER end member, not an impossible model.
# Rejecting it threw away precisely the interpretations that make S-C/S-E/S-F
# correct, and biased every transform hypothesis toward "saline".
#
# Now: within SIGMA_TOL of the AEM value, treat the excess as fresh water and
# floor σ_w. Beyond it, reject. SIGMA_TOL = 2.0 is set by the AEM's own spread,
# not chosen to make the answer come out.
const SIGMA_TOL = 2.0

# clay prior: logit-normal, median 0.15, broad (5–95 % ≈ 0.039–0.437)
const CLAY_MU = log(0.15 / 0.85)
# Widened from 0.9 to 1.2 after testing: at 0.9 the prior p95 was 0.36, which
# excludes the Shepparton truth of 0.46 altogether. A prior that cannot reach
# the truth makes the broad-vs-informed comparison in §5 trivially rigged. At
# 1.2 the 5–95 % range is 0.023–0.559 and P(clay > 0.4) = 0.134, so S-F still
# has to be carried by the data.
const CLAY_SD = 1.2
# m prior: one value per sample — a rock-physics exponent, not a depth profile
const M_LO, M_HI = 1.5, 2.2
# vertical correlation of the clay draw, matching the truth generator
const CLAY_CORR = 6.0

logistic(x) = 1 / (1 + exp(-x))

"""
Cholesky factor of a squared-exponential correlation matrix on `z`. Cached by
the caller — this is the expensive part and the grid never changes within a run.
"""
function corr_chol(z; corr = CLAY_CORR)
    n = length(z)
    K = [exp(-0.5 * ((z[i] - z[j]) / corr)^2) for i in 1:n, j in 1:n] + 1e-8I
    cholesky(Symmetric(K)).L
end

"""
One correlated clay profile from the prior. Correlated rather than
independent-per-cell on purpose: independent draws would average out over any
depth window, driving the window-mean to the prior mean and making an interval
hypothesis like S-F answerable by arithmetic rather than by data.
"""
draw_clay(L, rng; mu = CLAY_MU, sd = CLAY_SD) =
    logistic.(mu .+ sd .* (L * randn(rng, size(L, 1))))

"""
Logit-normal location vector for the §5 "lithology-informed" prior: centred on
a known clay profile instead of a single uninformative value. `sd` is tightened
by the caller. This is the comparison that shows how much the prior carries.
"""
informed_mu(clay_profile; floor = 0.01, ceil = 0.95) =
    [log(c / (1 - c)) for c in clamp.(clay_profile, floor, ceil)]

draw_m(rng) = M_LO + (M_HI - M_LO) * rand(rng)

"""
    invert(θ_SMR, σ_bulk, clay, m) -> (σ_w, θ_total, ok)

The transform, per cell. `ok` is false where the draw is inadmissible.
"""
function invert(θ_SMR, σ_bulk, clay, m; idx = eachindex(θ_SMR))
    n = length(θ_SMR)
    σw = fill(NaN, n); θt = fill(NaN, n); ok = true
    @inbounds for i in idx
        f = (1 - clay[i])^P_CLAY
        t = f <= 1e-6 ? Inf : θ_SMR[i] / f
        θt[i] = t
        if !(t > 0) || t > φ_MAX
            ok = false; σw[i] = NaN; continue
        end
        s = (σ_bulk[i] - σ_CLAY * clay[i]) / t^m
        if s < 0
            # clay over-explains: admissible as "fresh" while the shortfall is
            # inside the AEM's own uncertainty, rejected beyond it
            if σ_CLAY * clay[i] <= SIGMA_TOL * σ_bulk[i]
                σw[i] = σW_FLOOR; continue
            end
            ok = false; σw[i] = NaN; continue
        end
        σw[i] = max(s, σW_FLOOR)
    end
    σw, θt, ok
end

"""
    sample(θ_SMR, σ_bulk, L; rng, maxtries) -> (σ_w, clay, m, tries) or nothing

One admissible (σ_w, clay, m) for one posterior θ sample and one σ_bulk, by
rejection. Returns `nothing` if no admissible draw is found, which the caller
must count rather than silently drop.
"""
function sample(θ_SMR, σ_bulk, L; rng = Random.GLOBAL_RNG, maxtries = 200,
                idx = eachindex(θ_SMR), mu = CLAY_MU, sd = CLAY_SD)
    for t in 1:maxtries
        clay = draw_clay(L, rng; mu = mu, sd = sd); m = draw_m(rng)
        σw, θt, ok = invert(θ_SMR, σ_bulk, clay, m; idx = idx)
        ok && return (σw = σw, clay = clay, m = m, tries = t)
    end
    nothing
end

"""
    ensemble(Θ, σ_bulk, z; nper, seed)

Transform a whole posterior: for each column of `Θ` (θ_SMR draws) draw `nper`
admissible (clay, m) pairs against `σ_bulk`. Returns matrices of σ_w and clay
with one column per (θ draw × transform draw), plus the rejection rate.

`σ_bulk` may be a vector (fixed, conditions C3/C3b) or a matrix with one column
per θ draw (nested, condition C4).
"""
function ensemble(Θ, σ_bulk, z; nper = 4, seed = 20261012,
                  window = (0.0, Inf), L = nothing,
                  mu = CLAY_MU, sd = CLAY_SD)
    L === nothing && (L = corr_chol(z))
    idx = findall(zi -> window[1] <= zi <= window[2], z)
    rng = MersenneTwister(seed)
    nz, nd = size(Θ)
    σW = Matrix{Float32}(undef, nz, nd * nper)
    CL = Matrix{Float32}(undef, nz, nd * nper)
    MM = Vector{Float32}(undef, nd * nper)
    k = 0; failed = 0; tries = 0
    for j in 1:nd
        sb = σ_bulk isa AbstractMatrix ? view(σ_bulk, :, j) : σ_bulk
        for _ in 1:nper
            s = sample(view(Θ, :, j), sb, L; rng = rng, idx = idx, mu = mu, sd = sd)
            if s === nothing; failed += 1; continue; end
            k += 1; tries += s.tries
            @views σW[:, k] .= Float32.(s.σw)
            @views CL[:, k] .= Float32.(s.clay)
            MM[k] = Float32(s.m)
        end
    end
    (σw = σW[:, 1:k], clay = CL[:, 1:k], m = MM[1:k], idx = idx,
     n = k, failed = failed, mean_tries = k == 0 ? NaN : tries / k,
     accept = k == 0 ? 0.0 : k / (k + failed))
end

end # module
