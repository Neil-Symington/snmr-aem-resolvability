# joint_post.jl -- turn a finished joint chain into the files Stages 2-5 expect.
#
#   julia --project=. joint_post.jl 2P
#
# Stage 2 (run_conditions.jl) builds the SNMR kernel from
# aem<TAG>/aem_posterior.jld2, key "sigma50". The joint chain lives in
# aemdc<TAG>/, so this writes aem<TAG>J/aem_posterior.jld2 with the joint
# percentiles -- "Stage 2 on the joint median kernel" -- and leaves the
# AEM-only aem<TAG>/ untouched so every number in 5.2-5.5 still reproduces.
using Distributed, JLD2, Statistics, Random, Printf
const HERE = normpath(@__DIR__); const ROOT = normpath(joinpath(HERE,".."))
const OUT  = joinpath(HERE, "..", "data")
const TAG  = length(ARGS) >= 1 ? ARGS[1] : "2P"
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "dc_operator.jl"))
using HiQGA.transD_GP

function main()
    rundir = joinpath(OUT, "aemdc$(TAG)")
    S = load(joinpath(rundir, "aemdc_setup.jld2"))
    grid = build_aem_grid()
    opt = make_aem_opt(grid.znall, "hypsynth_aemdc_")
    old = pwd(); cd(rundir); opt.xall[:] .= grid.znall
    F = transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
    cd(old)
    log10ρ = hcat([f[:] for f in F]...); σ = 10 .^ (-log10ρ); n = size(σ,2)
    pct(p) = [quantile(view(σ,i,:), p) for i in axes(σ,1)]
    σ10, σ50, σ90 = pct(0.10), pct(0.50), pct(0.90)
    dst = joinpath(OUT, "aem$(TAG)J"); mkpath(dst)
    jldsave(joinpath(dst, "aem_posterior.jld2");
        σ10 = σ10, σ50 = σ50, σ90 = σ90, zall_aem = grid.zall,
        znall_aem = grid.znall, z_aem = grid.z, ndraws = n,
        σ_aem_truth = S["σ_aem_truth"], joint = true, parent_run = "aemdc$(TAG)")
    @printf("wrote aem%sJ/aem_posterior.jld2 from %d joint draws\n", TAG, n)
    # Figure 3 draws the FULL post-burn-in set, not the 600-draw subsample, so
    # the joint column needs the same file extract_full_draws.jl writes for the
    # AEM-only ones.
    jldsave(joinpath(OUT, "aem_draws_full_$(TAG)J.jld2");
        sigma = σ, zall_aem = grid.zall, column = TAG, joint = true,
        ndraws = n, parent_run = "aemdc$(TAG)")
    @printf("wrote aem_draws_full_%sJ.jld2 with %d draws\n", TAG, n)
    @printf("  sigma50 over 17.5-33.5 m: %.5f S/m\n",
            let k = findall(z -> 17.5 <= z <= 33.5, grid.zall); mean(σ50[k]); end)
end
main()
