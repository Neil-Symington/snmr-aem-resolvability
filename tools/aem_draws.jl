# aem_draws.jl — pull many intact conductivity profiles out of a finished AEM chain.
#
# Stage 2 has been reading the AEM posterior's p50 profile, which is a point
# estimate: it carries no uncertainty in sigma_bulk at all. That is harmless
# while nothing measures sigma_bulk, and wrong as soon as something does, because
# an induction log then has no posterior to constrain and can only SUBSTITUTE
# for the survey's answer. A substitution can make a posterior worse; a
# constraint cannot. This writes the draws so the log can be applied as a
# constraint instead.
#
# The production run stores only 30 draws, which is too coarse to survive the
# rejection a log imposes, so they are re-assembled here.
#
#   julia --project=. synthetics/aem_draws.jl [tag] [ndraws]

using Distributed, JLD2, Statistics, Random, Printf

const HERE = normpath(@__DIR__)
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(HERE, "..", "data")
const TAG    = length(ARGS) >= 1 ? ARGS[1] : "2M"
const NDRAW  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 600
# aem_HYP.jl's 4th argument writes to data/aem<TAG><LABEL>/ so a side run
# cannot overwrite the production posterior. This must take the same label or
# it looks in the wrong directory: a review hit exactly that, with aem_HYP.jl
# writing data/aemREV_smoke/ and this failing ENOENT on data/aemREV/.
const LABEL  = length(ARGS) >= 3 ? ARGS[3] : ""
const SEED   = 20261012
# The AEM noise realisation is pinned to HiQGA's makenoisydata! default, 123,
# NOT to SEED. Until build_skytem_operator passed rseed through, the caller's
# seed was silently ignored and every run drew the 123 realisation -- which is
# therefore the realisation behind every published number. Pinning it here
# keeps those reproducible now that the keyword is honoured; vary it only to
# test sensitivity to the noise draw, never to reproduce the paper.
const NOISE_RSEED = 123
include(joinpath(ROOT, "src", "aem_setup.jl"))
using HiQGA.transD_GP

function main()
    T = load(joinpath(OUT, "truth_murray$(TAG).jld2"))
    built = build_skytem_operator(T["σ_bulk"], T["z"]; noise_scale = 1.0, rseed_noise = NOISE_RSEED)
    grid = built.grid
    rundir = joinpath(OUT, "aem$(TAG)$(LABEL)")
    isdir(rundir) || error("no run directory $rundir -- if aem_HYP.jl was given " *
                           "a label, pass the same one here as the 3rd argument")
    opt = make_aem_opt(grid.znall, "hypsynth_aem_")
    opt.xall[:] .= grid.znall
    old = pwd(); cd(rundir)
    F = transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
    cd(old)
    n = length(F)
    Random.seed!(SEED)
    ix = n >= NDRAW ? sort(randperm(n)[1:NDRAW]) : collect(1:n)
    σ = hcat([10 .^ (-F[i][:]) for i in ix]...)        # S/m, per AEM cell
    @printf("%s: %d post-burn-in models, keeping %d  (%d cells)\n", TAG, n, length(ix), size(σ,1))
    jldsave(joinpath(OUT, "aem_draws_$(TAG).jld2");
            zall_aem = grid.zall, sigma = σ, draw_ix = ix, ndraws_total = n, seed = SEED)
    println("wrote aem_draws_$(TAG).jld2")
end

main()
