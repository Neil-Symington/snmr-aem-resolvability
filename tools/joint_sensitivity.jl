# joint_sensitivity.jl -- per-cell resolution for the JOINT data set.
#
#   julia --project=. joint_sensitivity.jl 2P
#
# The twin of aem_sensitivity.jl, with the sounding's rows stacked under the
# AEM's. Same definition throughout:
#     s_j    = || (W J)[:, j] ||        noise sigmas per decade of log10 sigma
#     dlog_j = 1 / s_j                  decades needed to move the data by one
# so panel (c) of Figure 3 is on the same scale as (a) and (b) and the strips
# are comparable.
#
# SIGNS MATTER in the cumulative form, which sums the Jacobian columns
# COHERENTLY. The AEM's J is d(data)/d(log10 sigma); the DC forward is written
# in log10 rho = -log10 sigma, so its numerical Jacobian is NEGATED before
# stacking. Get that wrong and the two data sets appear to cancel.
#
# Evaluated at the TRUTH model, like the AEM-only version: the linearisation
# point matters and the truth is the one choice that borrows from no inversion.
using JLD2, Statistics, Printf, LinearAlgebra
const HERE = normpath(@__DIR__); const ROOT = normpath(joinpath(HERE,".."))
const OUT  = joinpath(HERE, "..", "data")
const TAG  = length(ARGS) >= 1 ? ARGS[1] : "2P"
const SEED = 123
const DLOG_REF = 0.3
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "dc_operator.jl"))
using HiQGA.transD_GP
const SK = HiQGA.transD_GP.SkyTEM1DInversion
include(joinpath(ROOT, "src", "aem_jac.jl"))

"numerical d(rho_a)/d(log10 sigma_j), central differences"
function dc_jacobian(F::DCSchlum, m; h = 1e-4)
    n = length(m); J = zeros(F.ndata, n)
    for j in 1:n
        mp = copy(m); mp[j] += h; a = copy(forward!(F, mp))
        mm = copy(m); mm[j] -= h; b = copy(forward!(F, mm))
        J[:, j] = -(a .- b) ./ (2h)      # negated: m is log10 rho
    end
    J
end

function main()
    aemtag = TAG in ("2P","3P") ? TAG[1]*"M" : TAG
    T = load(joinpath(OUT, "truth_murray$(TAG).jld2"))
    B = build_aem_jac(T, SEED); aem, grid = B.aem, B.grid
    mtrue = log10.(B.σt[2:end])
    SK.getresidual(aem, mtrue; computeJ = true)
    WJa = aem.W * aem.J

    D = load(joinpath(OUT, "dc_sounding_$(TAG).jld2"))
    F = DCSchlum(D["rho_a"], D["sigma"], D["AB_2"], D["MN_2"], grid.z[2:end])
    Jd = dc_jacobian(F, -mtrue)                 # DC takes log10 rho
    WJd = Diagonal(1 ./ D["sigma"]) * Jd

    for (name, WJ) in (("AEM only", WJa), ("AEM + sounding", vcat(WJa, WJd)))
        s = [norm(view(WJ,:,j)) for j in axes(WJ,2)]
        dlog = 1 ./ s
        nc = size(WJ,2)
        cum = [norm(sum(view(WJ,:,j:nc), dims=2))*DLOG_REF for j in 1:nc]
        for j in 2:nc; cum[j] = min(cum[j], cum[j-1]); end
        k = findall(z -> 10 <= z <= 60, grid.zall)
        @printf("%-16s ndata %3d | dlog over 10-60 m: %.3f - %.3f (median %.3f)\n",
                name, size(WJ,1), minimum(dlog[k]), maximum(dlog[k]), median(dlog[k]))
        if name != "AEM only"
            jldsave(joinpath(OUT, "aem_sens_$(TAG)J.jld2"); zall = grid.zall,
                sigma_per_decade = s, dlog_for_1sigma = dlog, cum_sigma = cum,
                dlog_ref = DLOG_REF, tag = TAG*"J", seed = SEED,
                ndata = size(WJ,1), ndata_aem = size(WJa,1), ndata_dc = size(WJd,1),
                joint = true)
            println("wrote aem_sens_$(TAG)J.jld2")
        end
    end
end
main()
