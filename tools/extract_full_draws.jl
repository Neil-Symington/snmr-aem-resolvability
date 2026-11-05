# extract_full_draws.jl -- re-assemble ALL post-burn-in T=1 draws from the AEM
# chains, not the 600-draw subsample aem_draws_*.jld2 holds.
#
# The density figures were speckly because 600 draws over ~70 conductivity bins
# leaves single-figure counts per bin. The chains carry 2002 draws each and the
# files are still on disk, so the fix is more samples rather than smoothing the
# picture -- a smoothed histogram would be a drawing of the smoother, not of
# the sampler.
using HiQGA, JLD2, Printf, Statistics
using HiQGA.transD_GP
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "aem_setup.jl"))
const OUT = joinpath(HERE, "..", "data")

# An optional run label, matching aem_HYP.jl's 4th argument, which writes
# data/aem<TAG><LABEL>/. Without it a labelled Stage 1 run is invisible here.
const LABEL = length(ARGS) >= 1 ? ARGS[1] : ""

for tag in ("3M", "2M", "3S", "2S")
    rundir = joinpath(OUT, "aem$tag$LABEL")
    isdir(rundir) || (@printf("%-3s no chain directory\n", tag); continue)
    S = load(joinpath(rundir, "aem_setup.jld2"))
    znall = S["znall_aem"]; zall = S["zall_aem"]
    opt = make_aem_opt(znall, "hypsynth_aem_")
    old = pwd(); cd(rundir)
    F = try
        opt.xall[:] .= znall
        transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
    finally
        cd(old)
    end
    # fstar is log10(rho); the pipeline is conductivity-first
    σ = hcat([vec(10 .^ (-f)) for f in F]...)
    @printf("%-3s assembled %d draws over %d cells; sigma %.4g to %.4g S/m\n",
            tag, size(σ,2), size(σ,1), minimum(σ), maximum(σ))
    jldsave(joinpath(OUT, "aem_draws_full_$tag.jld2");
            sigma = σ, zall_aem = zall, ndraws = size(σ,2), source = "assembleTat1 T=1, burninfrac 0.5")
end
