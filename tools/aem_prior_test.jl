# aem_prior_test.jl -- is the AEM's trans-D posterior prior-dominated too?
#
#   julia --project=. aem_prior_test.jl <nsamples> <tag> <nchains> <nmax>
#
# A copy of aem_HYP.jl with ONE knob exposed: nmax. Everything else -- grid,
# bounds, GP, proposals, temperature ladder, noise seed -- is byte-identical to
# the production run, so a difference between two of these is the nuclei count
# and nothing else.
#
# THE CONTROL MATTERS AS MUCH AS THE TEST. The production chain is 200,001
# samples; this is 50,001. Comparing a short nmax-10 chain against a long
# nmax-40 one would confound the prior with chain length, so BOTH are run here
# at 50,001 and only the pair is compared. The production number is quoted
# separately and never differenced against these.
#
# Bounds are left at the AEM's 0.32-1778 Ohm.m deliberately. The DC re-prior
# changed nuclei AND bounds together; here only nuclei move, so if the answer
# shifts we know which knob did it.
using Distributed, JLD2, Statistics, Random, Printf
const HERE = normpath(@__DIR__); const ROOT = normpath(joinpath(HERE,".."))
const OUT  = joinpath(HERE, "..", "data")
const NSAMPLES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 50_001
const TAG      = length(ARGS) >= 2 ? ARGS[2] : "2M"
const NCHAINS  = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 4
const NMAX     = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 40
const TMAX     = 2.5
const SEED     = 20261012
# The AEM noise realisation is pinned to HiQGA's makenoisydata! default, 123,
# NOT to SEED. Until build_skytem_operator passed rseed through, the caller's
# seed was silently ignored and every run drew the 123 realisation -- which is
# therefore the realisation behind every published number. Pinning it here
# keeps those reproducible now that the keyword is honoured; vary it only to
# test sensitivity to the noise draw, never to reproduce the paper.
const NOISE_RSEED = 123
include(joinpath(ROOT, "src", "aem_setup.jl"))
using HiQGA.transD_GP
addprocs(NCHAINS - (nworkers() == 1 ? 0 : nworkers()))
@assert nworkers() == NCHAINS
@everywhere using Distributed, HiQGA.transD_GP

function main()
    T = load(joinpath(OUT, "truth_murray$(TAG).jld2"))
    built = build_skytem_operator(T["σ_bulk"], T["z"]; noise_scale = 1.0, rseed_noise = NOISE_RSEED)
    aem, grid = built.aem, built.grid
    rundir = joinpath(OUT, "aem$(TAG)_n$(NMAX)"); mkpath(rundir)
    for f in readdir(rundir); occursin("hypsynth_aem_", f) && rm(joinpath(rundir,f), force=true); end
    @printf("AEM-only prior test: column %s, nmax %d, %d samples, %d chains\n",
            TAG, NMAX, NSAMPLES, NCHAINS)
    opt = make_aem_opt(grid.znall, "hypsynth_aem_"; nmax = NMAX)
    old = pwd(); cd(rundir); t0 = time()
    try
        transD_GP.main(opt, aem; Tmax = TMAX, nsamples = NSAMPLES,
                       nchains = NCHAINS, nchainsatone = 1)
    finally; cd(old); end
    @printf("MCMC done in %.1f min\n", (time()-t0)/60)
    cd(rundir); opt.xall[:] .= grid.znall
    F = transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
    cd(old)
    n = length(F); log10ρ = hcat([f[:] for f in F]...); σ = 10 .^ (-log10ρ)
    jldsave(joinpath(OUT, "aem_draws_$(TAG)_n$(NMAX).jld2");
        sigma = σ, zall_aem = grid.zall, column = TAG, nmax = NMAX,
        nsamples = NSAMPLES, ndraws = n)
    @printf("wrote aem_draws_%s_n%d.jld2 with %d draws\n", TAG, NMAX, n)
end
main()
