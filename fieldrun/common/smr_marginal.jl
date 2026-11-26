# smr_marginal.jl -- the SNMR noise scale s integrated out of the likelihood,
# instead of sampled. Selected by [noise] s_treatment = "marginal".
#
# THE MODEL IS UNCHANGED. SMRPInversion's SMRSoundingHierarchical (amplitude
# only) has, per data residual r_i and noise shape sigma-hat_i,
#
#     U(w, u) = 1/2 [ R(w) e^-u + N u + (u - mu)^2 / tau^2 ],   u = log s,
#     R(w)    = sum_i r_i^2 / sigma-hat_i^2,
#
# and HiQGA samples u as a nuisance, uniform within [u_lo, u_hi], alongside w.
# Here u is integrated out instead, so the sampler sees only w:
#
#     U_marg(w) = -log  integral_{u_lo}^{u_hi} exp(-U(w, u)) du .
#
# The posterior on w is the same; what changes is the sampler. Sampling u adds
# a coordinate tied to the misfit and the layer count (more layers fit better,
# so s falls), which a random walk in u crawls along, and an extra round of
# worker calls every iteration. Both are what the replicate check pointed at:
# s, misfit and k were the slowest scalars, and s drifted in 26 of 32 ladders.
#
# The integral is a trapezoid sum on a fixed grid in u (spacing 0.01, against a
# conditional width of about sqrt(2/N) ~ 0.3), evaluated with logsumexp. With a
# flat prior on u and wide bounds it is analytic -- Gamma(N/2) (2/R)^(N/2), the
# Jeffreys-prior marginal Ray and Taylor's SMRSoundingUnknown uses -- and the
# tests check the quadrature against that.
#
# s AFTER THE RUN. Given w, p(u | R) is one-dimensional, so every stored model
# gets a posterior draw of u from the same grid (draw_log_s). Detection and the
# report read those draws exactly as they read the sampled ones. The chains
# store U_marg, not R; U_marg is strictly increasing in R, so R is recovered
# from a lookup table (rss_from_U).
#
# HiQGA tempers the whole misfit, so the heated chains here target
# exp(-U_marg / T), not the sampled version's exp(-U / T). Only T = 1 is
# inference, and that is the same distribution in both.

module SMRMarginal

using SMRPInversion, SNMRForward, Random
import SMRPInversion: get_misfit

export SMRSoundingMarginal, newSMRSoundingMarginal, MarginalQuad,
       marginal_U, rss_from_U, draw_log_s, mean_log_s

"The u = log s quadrature: grid, log weights (trapezoid + log prior), N."
struct MarginalQuad
    u    :: Vector{Float64}
    logw :: Vector{Float64}
    h    :: Float64
    N    :: Int
    # lookup table for inverting U_marg(R): log R grid and U on it
    logR :: Vector{Float64}
    Utab :: Vector{Float64}
end

logsumexp(x) = (m = maximum(x); isfinite(m) ? m + log(sum(exp.(x .- m))) : m)

_lterms(Q, R) = Q.logw .- 0.5 .* (R .* exp.(-Q.u) .+ Q.N .* Q.u)

function _U(u, logw, N, R)
    -logsumexp(logw .- 0.5 .* (R .* exp.(-u) .+ N .* u))
end

"""
    MarginalQuad(N, mu, tau, ubounds; h = 0.01)

`mu`, `tau`: the log-normal prior on s (mu = log of the centre). `tau = Inf`
gives a flat prior on u = log s. `ubounds` = (u_lo, u_hi), natural log.
"""
function MarginalQuad(N::Integer, mu::Real, tau::Real, ubounds; h::Real = 0.01)
    lo, hi = Float64.(ubounds)
    lo < hi || throw(ArgumentError("u bounds must be lo < hi, got $ubounds"))
    n = max(2, ceil(Int, (hi - lo) / h)) + 1
    u = collect(range(lo, hi, length = n)); hh = step(range(lo, hi, length = n))
    wt = fill(hh, n); wt[1] = wt[end] = hh / 2
    prior = isinf(tau) ? zeros(n) : -0.5 .* (u .- mu) .^ 2 ./ tau^2
    logw = log.(wt) .+ prior
    # R = N s spans s over the bounds; pad by a wide margin either side
    logR = collect(range(log(N) + lo - 15, log(N) + hi + 15, length = 6000))
    Utab = [_U(u, logw, N, exp(l)) for l in logR]
    all(diff(Utab) .> 0) || error("U_marg is not increasing on the lookup grid")
    MarginalQuad(u, logw, hh, N, logR, Utab)
end

"U_marg for a residual sum of squares R (in sigma-hat units)."
marginal_U(Q::MarginalQuad, R::Real) = _U(Q.u, Q.logw, Q.N, R)

"Invert U_marg: the R that gave a stored misfit."
function rss_from_U(Q::MarginalQuad, U::Real)
    T = Q.Utab
    (U <= T[1] || U >= T[end]) && error("misfit $U outside the lookup table")
    i = searchsortedlast(T, U)
    f = (U - T[i]) / (T[i+1] - T[i])
    l = Q.logR[i] + f * (Q.logR[i+1] - Q.logR[i])
    # one Newton step in log R tidies the linear interpolation
    R = exp(l)
    p = exp.(_lterms(Q, R) .- maximum(_lterms(Q, R))); p ./= sum(p)
    dUdlogR = 0.5 * R * sum(p .* exp.(-Q.u))
    exp(l - (marginal_U(Q, R) - U) / dUdlogR)
end
rss_from_U(Q::MarginalQuad, U::AbstractVector) = [rss_from_U(Q, x) for x in U]

"Posterior weights of the u grid given R."
function _post(Q, R)
    lt = _lterms(Q, R); p = exp.(lt .- maximum(lt)); p ./ sum(p)
end

"One posterior draw of u = log s given R (grid cell, then uniform within it)."
function draw_log_s(Q::MarginalQuad, R::Real, rng::AbstractRNG)
    p = _post(Q, R)
    c = cumsum(p); i = min(searchsortedfirst(c, rand(rng)), length(c))
    clamp(Q.u[i] + (rand(rng) - 0.5) * Q.h, Q.u[1], Q.u[end])
end
draw_log_s(Q::MarginalQuad, R::AbstractVector, rng::AbstractRNG) =
    [draw_log_s(Q, x, rng) for x in R]

"E[u | R], the Rao-Blackwellised log s: a function of w alone."
mean_log_s(Q::MarginalQuad, R::Real) = sum(_post(Q, R) .* Q.u)
mean_log_s(Q::MarginalQuad, R::AbstractVector) = [mean_log_s(Q, x) for x in R]

"""
Amplitude-only SMR sounding with s integrated out. The fields HiQGA and
SMRPInversion's generic code use (V0, ϕ, Fm, linearsat, amponly, offset_ϕ,
stretch_ϕ) are as in the other sounding types.
"""
mutable struct SMRSoundingMarginal <: SMRPInversion.SMRSounding
    V0        :: Vector{Float64}
    ϕ         :: Vector{Float64}
    Fm        :: SNMRForward.MRSForward
    σ̂²_V0     :: Vector{Float64}
    Q         :: MarginalQuad
    linearsat :: Bool
    amponly   :: Bool
    offset_ϕ  :: Float64
    stretch_ϕ :: Float64
end

"""
    newSMRSoundingMarginal(V0, ϕ, Fm; σ̂_V0, s_prior_center, s_prior_logwidth,
                           s_bounds_ln, linearsat = true)

The same arguments as `newSMRSoundingHierarchical`, plus the natural-log
bounds HiQGA would have sampled u within.
"""
function newSMRSoundingMarginal(V0, ϕ, Fm; σ̂_V0::AbstractVector,
                                s_prior_center::Real, s_prior_logwidth::Real,
                                s_bounds_ln, linearsat::Bool = true)
    length(σ̂_V0) == length(V0) || throw(ArgumentError("σ̂_V0 must match V0"))
    Q = MarginalQuad(length(V0), log(s_prior_center), s_prior_logwidth, s_bounds_ln)
    SMRSoundingMarginal(Float64.(V0), Float64.(ϕ), Fm, Float64.(σ̂_V0) .^ 2, Q,
                        linearsat, true, 0.0, 1.0)
end

function get_misfit(w::Vector{<:Real}, S::SMRSoundingMarginal)
    Vres, _ = SMRPInversion.applynu(S, w)
    R = sum((S.V0 .- Vres) .^ 2 ./ S.σ̂²_V0)
    isfinite(R) || return Inf
    marginal_U(S.Q, R)
end

end # module
