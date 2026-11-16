# test_stage2_pool.jl -- pooling ONE Stage 1 folder must reproduce that folder's
# own sigma50 exactly.
#
#   julia --project=. test/test_stage2_pool.jl [folder]
#
# This is the check that the pooling path reads the same draws, with the same
# burn-in, and takes the same quantile as Stage 1 does. If it passes, pooling
# several folders is the same operation over a larger set.
using JLD2, Statistics, Printf
const HERE=@__DIR__; const ROOT=normpath(joinpath(HERE,".."))
const OUT = joinpath(ROOT,"data")
const GADI_SRC = joinpath(ROOT,"src")
include(joinpath(ROOT,"src","aem_setup.jl"))
using HiQGA.transD_GP

function pooled_sigma(folders)      # identical to stages/stage2_snmr.jl
    mats = Matrix{Float64}[]; zall = nothing; perfolder = Int[]
    for d in folders
        S = load(joinpath(d,"aem_setup.jld2")); znall = S["znall_aem"]
        zall === nothing && (zall = collect(Float64, S["zall_aem"]))
        opt = make_aem_opt(znall, "hypsynth_aem_"; nmax = get(S,"nmax_k",40))
        opt.xall[:] .= znall
        old = pwd(); cd(d)
        F = try transD_GP.assembleTat1(opt, :fstar; temperaturenum=1, burninfrac=0.5)
            finally cd(old) end
        m = 10 .^ (-hcat([f[:] for f in F]...))
        push!(mats, m); push!(perfolder, size(m,2))
    end
    σ = hcat(mats...)
    ([quantile(view(σ,i,:),0.50) for i in axes(σ,1)], zall, size(σ,2), perfolder)
end

folder = length(ARGS) >= 1 ? ARGS[1] :
         first(filter(d -> isdir(joinpath(OUT,d)) &&
                           isfile(joinpath(OUT,d,"aem_posterior.jld2")) &&
                           startswith(d,"aem"), readdir(OUT)))
d = isabspath(folder) ? folder : joinpath(OUT, folder)
@printf("pooling ONE folder: %s\n", basename(d))
σp, zp, ntot, perf = pooled_sigma([d])
A = load(joinpath(d,"aem_posterior.jld2")); σ50 = A["σ50"]
@printf("  pooled draws %d (Stage 1 recorded %d)\n", ntot, A["ndraws_total"])
@printf("  lengths: pooled %d, stored sigma50 %d\n", length(σp), length(σ50))
exact = length(σp)==length(σ50) && all(σp .== σ50)
@printf("  EXACTLY equal: %s\n", exact)
if !exact && length(σp)==length(σ50)
    @printf("  max abs diff %.3e   max rel diff %.3e\n",
            maximum(abs.(σp .- σ50)), maximum(abs.((σp .- σ50) ./ σ50)))
end
if exact
    println("\nPASS: pooling one folder reproduces the single-folder kernel exactly.")
else
    println("\nFAIL"); exit(1)
end
