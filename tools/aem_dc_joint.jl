# aem_dc_joint.jl -- Stage 1, joint. SkyTEM AEM + Schlumberger sounding in one
# trans-dimensional chain.
#
#   julia --project=. aem_dc_joint.jl [nsamples] [tag]
#
# Identical to aem_HYP.jl in every setting -- 200,001 samples, 4 chains,
# Tmax 2.5, nchainsatone 1, same grid, same prior, same noise seed -- except
# that the operator handed to transD_GP.main is a JointOperator wrapping the
# SkyTEM dBzdt and the DCSchlum sounding. GenericJointInversion.get_misfit
# sums the two log-likelihoods; neither operator knows about the other.
#
# The AEM-only chains are NOT touched: this writes to aemdc<TAG>/ and the
# draws to aem_draws_<TAG>J.jld2, so every AEM-only number in 5.2-5.5 still
# comes from the files that produced it.
using Distributed, JLD2, Statistics, Random, Printf

const HERE = normpath(@__DIR__)
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(HERE, "..", "data")
const NSAMPLES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 200_001
const TAG      = length(ARGS) >= 2 ? ARGS[2] : "2P"
const NCHAINS  = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 4
const TMAX     = 2.5
const SEED     = 20261012        # the AEM noise seed, unchanged
const NDRAW    = 600

const GADI_SRC = joinpath(ROOT, "src")
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "dc_operator.jl"))
using HiQGA.transD_GP

addprocs(NCHAINS - (nworkers() == 1 ? 0 : nworkers()))
@assert nworkers() == NCHAINS "need $NCHAINS workers, got $(nworkers())"
@everywhere using Distributed, HiQGA, HiQGA.transD_GP
@everywhere include(joinpath(@__DIR__, "dc_operator.jl"))

function main()
    # the AEM column the P tier uses: 3P/2P share the M lineage's AEM data
    aemtag = TAG in ("2P","3P") ? TAG[1]*"M" : TAG
    T = load(joinpath(OUT, "truth_murray$(TAG).jld2"))
    built = build_skytem_operator(T["σ_bulk"], T["z"]; noise_scale = 1.0,
                                  rseed_noise = SEED)
    aem, grid = built.aem, built.grid
    D = load(joinpath(OUT, "dc_sounding_$(TAG).jld2"))
    dc = DCSchlum(D["rho_a"], D["sigma"], D["AB_2"], D["MN_2"], grid.z[2:end])
    mtrue = log10.(built.ρ_aem_truth[2:end])

    @printf("JOINT Stage 1, column %s\n", TAG)
    @printf("  AEM : %d LM + %d HM gates\n", length(aem.dlow), length(aem.dhigh))
    @printf("  DC  : %d Schlumberger points, AB/2 %.1f-%.1f m, MN/2 = AB/10\n",
            dc.ndata, minimum(D["AB_2"]), maximum(D["AB_2"]))
    @printf("  DC chi2/N at the truth on this grid: %.4f  (noise realisation alone %.4f)\n",
            chi2byN(dc, mtrue), mean(((D["rho_a"] .- D["rho_a_clean"])./D["sigma"]).^2))

    F = transD_GP.GenericJointInversion.JointOperator([aem, dc])

    rundir = joinpath(OUT, "aemdc$(TAG)"); mkpath(rundir)
    for f in readdir(rundir)
        occursin("hypsynth_aemdc_", f) && rm(joinpath(rundir, f), force = true)
    end
    jldsave(joinpath(rundir, "aemdc_setup.jld2");
        σ_aem_truth = built.σ_aem_truth, ρ_aem_truth = built.ρ_aem_truth,
        z_aem = grid.z, zboundaries_aem = grid.zbnd,
        znall_aem = grid.znall, zall_aem = grid.zall,
        dlow = aem.dlow, dhigh = aem.dhigh, σlow = aem.σlow, σhigh = aem.σhigh,
        dc_AB_2 = D["AB_2"], dc_MN_2 = D["MN_2"], dc_rho_a = D["rho_a"],
        dc_sigma = D["sigma"], dc_seed = D["seed"],
        nsamples = NSAMPLES, nchains = NCHAINS, Tmax = TMAX,
        noise_scale = 1.0, rseed_noise = SEED, column = TAG, joint = true)

    opt = make_aem_opt(grid.znall, "hypsynth_aemdc_")
    old = pwd(); cd(rundir)
    t0 = time()
    try
        transD_GP.main(opt, F; Tmax = TMAX, nsamples = NSAMPLES,
                       nchains = NCHAINS, nchainsatone = 1)
    finally
        cd(old)
    end
    @printf("MCMC done in %.1f min\n", (time() - t0) / 60)

    cd(rundir); opt.xall[:] .= grid.znall
    Fs = transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
    cd(old)
    n = length(Fs)
    @printf("assembled %d post-burn-in draws at T=1\n", n)
    log10ρ = hcat([f[:] for f in Fs]...)
    σall = 10 .^ (-log10ρ)

    # posterior-predictive chi2/N, each data set separately
    nsub = min(n, 400); ix = round.(Int, range(1, n, length = nsub))
    c2dc = Float64[]; c2aem = Float64[]
    for i in ix
        m = log10ρ[:, i]
        push!(c2dc, chi2byN(dc, m))
        transD_GP.SkyTEM1DInversion.getfield!(m, aem)
        rl = (aem.Flow.dBzdt .- aem.dlow) ./ aem.σlow
        rh = (aem.Fhigh.dBzdt .- aem.dhigh) ./ aem.σhigh
        nd = count(aem.selectlow) + count(aem.selecthigh)
        push!(c2aem, (sum(rl[aem.selectlow].^2) + sum(rh[aem.selecthigh].^2)) / nd)
    end
    @printf("POSTERIOR-PREDICTIVE chi2/N over %d draws\n", nsub)
    @printf("  AEM  p50 %.3f   p5-p95 %.3f-%.3f\n", median(c2aem), quantile(c2aem,.05), quantile(c2aem,.95))
    @printf("  DC   p50 %.3f   p5-p95 %.3f-%.3f\n", median(c2dc), quantile(c2dc,.05), quantile(c2dc,.95))

    Random.seed!(SEED)
    sel = n >= NDRAW ? sort(randperm(n)[1:NDRAW]) : collect(1:n)
    jldsave(joinpath(OUT, "aem_draws_$(TAG)J.jld2");
        sigma = σall[:, sel], zall_aem = grid.zall, column = TAG, joint = true,
        chi2byN_aem = c2aem, chi2byN_dc = c2dc, ndraws = length(sel),
        parent_run = "aemdc$(TAG)")
    @printf("wrote aem_draws_%sJ.jld2 with %d draws\n", TAG, length(sel))
end
main()
