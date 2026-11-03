# aem_sensitivity.jl — per-cell AEM sensitivity, from the operator's own
# Jacobian. Nothing about the sampler enters.
#
# aem.J is d(data)/d(log10 sigma) and aem.W is 1/sigma_noise, so
#     s_j = || (W J)[:, j] ||
# is the number of noise standard deviations the data move for a ONE-DECADE
# change in cell j alone, and
#     dlog_j = 1 / s_j
# is the change in log10 conductivity that cell would need, by itself, to move
# the data by one standard deviation. It is the AEM twin of the water content
# `constraint_depth` asks for in the SNMR kernel.
#
# Evaluated at the TRUTH model: the linearisation point matters (the problem is
# not linear) and the truth is the one choice that does not borrow from any
# inversion.
#
#   julia --project=. synthetics/aem_sensitivity.jl <2M|3M>
using JLD2, Statistics, Printf, LinearAlgebra

const HERE = normpath(@__DIR__)
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(HERE, "..", "data")
const TAG  = length(ARGS) >= 1 ? ARGS[1] : "2M"
const SEED = 123
const DLOG_REF = 0.3        # decades of log10 conductivity = a factor of 2
include(joinpath(ROOT, "src", "aem_setup.jl"))
using HiQGA.transD_GP
const SK = HiQGA.transD_GP.SkyTEM1DInversion
include(joinpath(ROOT, "src", "aem_jac.jl"))

function main()
    T = load(joinpath(OUT, "truth_murray$(TAG).jld2"))
    B = build_aem_jac(T, SEED); aem, grid = B.aem, B.grid
    zall = grid.zall
    mtrue = log10.(B.σt[2:end])
    SK.getresidual(aem, mtrue; computeJ = true)
    WJ = aem.W * aem.J
    s = [norm(view(WJ, :, j)) for j in axes(WJ, 2)]
    dlog = 1 ./ s
    # CUMULATIVE sensitivity, the depth-of-investigation form (Christiansen &
    # Auken 2012): perturb EVERYTHING below cell j by DLOG_REF decades at once
    # and ask how many noise standard deviations of data change that produces.
    # Summed coherently, because the responses add in the data, not in quadrature.
    # Unlike the per-cell number this falls monotonically with depth, so it
    # measures how far the system sees — not whether one particular cell happens
    # to sit in a notch of the kernel.
    nc = size(WJ, 2)
    cum = [norm(sum(view(WJ, :, j:nc), dims = 2)) * DLOG_REF for j in 1:nc]
    for j in 2:nc; cum[j] = min(cum[j], cum[j-1]); end   # once lost, stays lost
    @printf("AEM sensitivity, column %s, %d data, %d cells\n", TAG, length(aem.res), length(s))
    @printf("%8s %12s %12s %14s\n", "z (m)", "sigma/decade", "decades for 1s", "cumulative (s)")
    for j in eachindex(s)
        zall[j] > 130 && break
        zall[j] > 0 && @printf("%8.1f %12.2f %12.3f %14.2f\n", zall[j], s[j], dlog[j], cum[j])
    end
    jldsave(joinpath(OUT, "aem_sens_$(TAG).jld2"); zall = zall, sigma_per_decade = s,
            dlog_for_1sigma = dlog, cum_sigma = cum, dlog_ref = DLOG_REF,
            tag = TAG, seed = SEED, ndata = length(aem.res))
    println("wrote aem_sens_$(TAG).jld2")
end
main()
