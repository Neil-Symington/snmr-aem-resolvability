# rockphys.jl — the petrophysical transform of hypotheses/transform.jl, with the
# clay surface conductivity lifted from a hardcoded constant to a prior.
#
# The forward model is unchanged:
#
#     theta_SMR = theta_total * (1 - clay)^p          p = P_CLAY
#     sigma_bulk = sigma_w * theta_total^m + sigma_clay * clay
#
# inverted for the pore-water conductivity:
#
#     theta_total = theta_SMR / (1 - clay)^p
#     sigma_w     = (sigma_bulk - sigma_clay * clay) / theta_total^m
#
# WHY sigma_clay IS A PRIOR HERE. transform.jl fixes sigma_clay = 0.5 S/m, an
# illite-ish mid value. stage2_salinity.jl made it settable and made the PRIOR
# the default, for the reason given in its own header: cation exchange capacity
# is essentially never measured, natural gamma responds to illite far more than
# to kaolinite or smectite, and a fixed value is "an assumption dressed as a
# measurement". On the synthetics it turned a wrong answer into a confident one.
# There is no mineralogy at all for these East Kimberley sites, so the prior is
# the only defensible choice.
#
# Everything else -- the correlated logit-normal clay prior (median 0.15,
# sd 1.2, 6 m correlation length), m ~ U(1.5, 2.2), the porosity ceiling of 0.5,
# the SIGMA_TOL = 2.0 fresh-end handling -- is transform.jl's, unchanged, and is
# imported from it rather than copied so the two cannot drift apart.

module RockPhys

using Random, Statistics, LinearAlgebra

# Vendored beside this file. The original include pointed at
# ../hypotheses/transform.jl, which is OUTSIDE the repository and was never
# shipped -- the same defect that made Stage 2 unrunnable from a clone for two
# rounds, found here before anyone tried to run it. transform.jl is a
# self-contained module with no includes and no absolute paths, so it moves as
# one file.
include(joinpath(@__DIR__, "transform.jl"))
using .Transform: corr_chol, draw_clay, draw_m, logistic,
                  P_CLAY, φ_MAX, σW_FLOOR, SIGMA_TOL, CLAY_MU, CLAY_SD, CLAY_CORR

export ClayPrior, invert_cell, sample_one, ensemble, MGL

# mg/L of total dissolved solids per S/m of pore-water conductivity -- the
# convention used throughout this project (gate_threshold.jl).
const MGL = 6400.0

"""
How the surface-conduction term is drawn.

`mode = :sigclay`  log-uniform sigma_clay over [lo, hi] S/m, multiplied by the
                   drawn clay fraction. Default 0.10-2.0 spans kaolinite to
                   smectite, i.e. no knowledge of clay type.
`mode = :product`  log-uniform over the PRODUCT sigma_clay * clay directly, on
                   [lo, hi] S/m. This is what the report's verdict tables quote
                   (the broad product prior, 0.002-0.8): the product is what the
                   conductivity equation actually contains, and it is what a
                   single calibration bore measures as a residual.
`mode = :none`     surface conduction set to zero, but the clay fraction is
                   still drawn, so it still inflates theta_total through the
                   bound-water term. This is "no clay CREDIT for conductivity",
                   not the upper limit on sigma_w.
`mode = :zero`     clay identically zero: theta_total = theta_SMR AND surf = 0.
                   This is the conservative, prior-free limit and the one the
                   gate uses. Both substitutions push sigma_w UP -- a smaller
                   theta_total means a larger sigma_w = sigma_bulk/theta^m, and
                   no surface term means all of sigma_bulk is pore water -- so a
                   site that looks fresh here is fresh under any clay model.
"""
struct ClayPrior
    mode::Symbol
    lo::Float64
    hi::Float64
end
ClayPrior() = ClayPrior(:sigclay, 0.10, 2.0)

loguniform(rng, lo, hi) = exp(log(lo) + rand(rng) * (log(hi) - log(lo)))

"Surface-conduction term per cell for one draw, given a clay profile."
function surface_term(cp::ClayPrior, clay, rng)
    if cp.mode === :none || cp.mode === :zero
        return zeros(length(clay)), 0.0
    elseif cp.mode === :product
        s = loguniform(rng, cp.lo, cp.hi)
        return fill(s, length(clay)), s
    else
        sc = loguniform(rng, cp.lo, cp.hi)
        return sc .* clay, sc
    end
end

"""
    invert_cell(θ_SMR, σ_bulk, clay, surf, m; idx, p)

transform.jl's `invert`, with the surface-conduction term passed in per cell
instead of computed from a fixed sigma_clay. Returns (sigma_w, theta_total, ok).
"""
function invert_cell(θ_SMR, σ_bulk, clay, surf, m; idx = eachindex(θ_SMR),
                     p = P_CLAY)
    n = length(θ_SMR)
    σw = fill(NaN, n); θt = fill(NaN, n); ok = true
    @inbounds for i in idx
        f = (1 - clay[i])^p
        t = f <= 1e-6 ? Inf : θ_SMR[i] / f
        θt[i] = t
        if !(t > 0) || t > φ_MAX
            ok = false; σw[i] = NaN; continue
        end
        s = (σ_bulk[i] - surf[i]) / t^m
        if s < 0
            # The clay term over-explains the conductivity. That is the
            # FRESH-WATER end member, not an impossible model, so while the
            # shortfall sits inside the AEM's own spread it is floored rather
            # than rejected (transform.jl header: rejecting it biased every
            # hypothesis toward saline).
            if surf[i] <= SIGMA_TOL * σ_bulk[i]
                σw[i] = σW_FLOOR; continue
            end
            ok = false; σw[i] = NaN; continue
        end
        σw[i] = max(s, σW_FLOOR)
    end
    σw, θt, ok
end

"One admissible draw of (sigma_w, clay, surf, m) by rejection."
function sample_one(θ_SMR, σ_bulk, L, cp; rng = Random.GLOBAL_RNG, maxtries = 200,
                    idx = eachindex(θ_SMR), p = P_CLAY)
    for t in 1:maxtries
        clay = cp.mode === :zero ? zeros(length(θ_SMR)) : draw_clay(L, rng)
        m = draw_m(rng)
        surf, sc = surface_term(cp, clay, rng)
        σw, θt, ok = invert_cell(θ_SMR, σ_bulk, clay, surf, m; idx = idx, p = p)
        ok && return (σw = σw, θt = θt, clay = clay, surf = surf, sc = sc,
                      m = m, tries = t)
    end
    nothing
end

"""
    ensemble(Θ, σ_bulk, z, cp; nper, seed, window, p)

Transform a whole SNMR posterior. `σ_bulk` is either a vector (one conductivity
profile, e.g. the Occam model) or a matrix with one column per theta draw, which
is the correct pairing when AEM posterior draws are available: an AEM percentile
profile is not an AEM draw, and propagating one in place of the other was a
measured defect in this project.

Returns per-draw matrices and the acceptance rate. Rejections are counted, not
silently dropped: a low acceptance rate means the prior and sigma_bulk disagree
over that window, which is itself a result.
"""
function ensemble(Θ, σ_bulk, z, cp::ClayPrior; nper = 4, seed = 20261026,
                  window = (0.0, Inf), L = nothing, p = P_CLAY)
    L === nothing && (L = corr_chol(z))
    idx = findall(zi -> window[1] <= zi <= window[2], z)
    rng = MersenneTwister(seed)
    nz, nd = size(Θ)
    σW = Matrix{Float32}(undef, nz, nd * nper)
    CL = Matrix{Float32}(undef, nz, nd * nper)
    TT = Matrix{Float32}(undef, nz, nd * nper)
    SF = Matrix{Float32}(undef, nz, nd * nper)
    MM = Vector{Float32}(undef, nd * nper)
    SC = Vector{Float32}(undef, nd * nper)
    k = 0; failed = 0; tries = 0
    for j in 1:nd
        sb = σ_bulk isa AbstractMatrix ? view(σ_bulk, :, j) : σ_bulk
        for _ in 1:nper
            s = sample_one(view(Θ, :, j), sb, L, cp; rng = rng, idx = idx, p = p)
            if s === nothing; failed += 1; continue; end
            k += 1; tries += s.tries
            @views σW[:, k] .= Float32.(s.σw)
            @views CL[:, k] .= Float32.(s.clay)
            @views TT[:, k] .= Float32.(s.θt)
            @views SF[:, k] .= Float32.(s.surf)
            MM[k] = Float32(s.m); SC[k] = Float32(s.sc)
        end
    end
    (σw = σW[:, 1:k], clay = CL[:, 1:k], θt = TT[:, 1:k], surf = SF[:, 1:k],
     m = MM[1:k], sigclay = SC[1:k], idx = idx, n = k, failed = failed,
     mean_tries = k == 0 ? NaN : tries / k,
     accept = k == 0 ? 0.0 : k / (k + failed))
end

"""
    unit_sigma_w(E, idx)

Geometric mean pore-water conductivity over a window, per draw. Geometric to
match the functional the AEM profile constrains (the forced constraint in
occam_coverage.jl fixes the window's MEAN LOG10 conductivity), so the two sides
of the gate are the same kind of average.
"""
unit_sigma_w(E, idx) = vec(10 .^ mean(log10.(max.(E.σw[idx, :], σW_FLOOR)), dims = 1))

"Geometric-mean total water content over a window, per draw."
unit_theta_total(E, idx) = vec(10 .^ mean(log10.(max.(E.θt[idx, :], 1e-6)), dims = 1))

end # module
