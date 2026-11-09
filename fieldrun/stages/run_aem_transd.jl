# Ported from udf_realdata/ into the repository. The include paths and the
# project are the only things changed from the working version: COMMON is
# fieldrun/common/ and the project is the repository root, not ../ek_prerun.
# The control file first: support modules read it at load time.
include(joinpath(@__DIR__, "..", "common", "stage_prelude.jl"))
# run_aem_transd.jl — Stage 1: trans-dimensional AEM inversion, UDF
# returning an ENSEMBLE of conductivity profiles.
#
#   julia --project=. run_aem_transd.jl --site NAME [--nsamples N] [--nchains N]
#
# This replaces the smooth Occam model as the Stage-1 product. The paper is
# explicit that Stage 1 is trans-D and that its output is an ensemble, which
# Stage 2 then draws from per SNMR sample so kernel uncertainty is propagated.
# The Occam run is still needed first, for two things only: the calibrated noise
# scale, and a sanity reference for the posterior.
#
# Modelled on synthetics/aem_HYP.jl. Two deliberate differences:
#
#  1. `make_aem_opt` is reproduced here rather than imported from
#     gadi/src/aem_setup.jl, because that file includes HiQGA's
#     electronics_halt at top level, which would define LM_times, HM_times,
#     rTx and friends in Main and collide with the UDF system in
#     udf_system.jl. Only the OptionsStat construction is wanted.
#
#  2. The resistivity bounds are set for THIS ground, not Ord-Keep's. The UDF
#     floodplain is conductive: the smooth models span 0.4 to 104 ohm.m, and
#     SMR_UD25 reaches within 0.1 dex of the Occam conductive bound. So the
#     prior runs log10 rho in [-1.0, 3.0] (0.1 to 1,000 ohm.m), which leaves
#     headroom at the conductive end where this ground actually lives and does
#     not spend prior on the 10,000 ohm.m that Ord-Keep needed and this survey
#     never sees.

using Distributed, Printf

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        if startswith(ARGS[i], "--") && i < length(ARGS) && !startswith(ARGS[i+1], "--")
            ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2
        elseif startswith(ARGS[i], "--")
            ARG[ARGS[i][3:end]] = "true"; i += 1
        else; i += 1; end
    end
end
const SITE     = get(ARG, "site", "")
const NSAMPLES = NSAMPLES_AEM   # [mcmc] nsamples_aem, or --nsamples
const NCHAINS  = NCHAINS_CFG    # from [mcmc], or --nchains
const TMAX     = 2.5
const NDRAWS   = 200           # intact profiles retained for Stage 2 and Stage 5
const SEED     = 20261028
const FBOUNDS  = [-1.0 3.0]    # log10 resistivity, set from the smooth models

isempty(SITE) && error("--site is required")

addprocs(NCHAINS - (nworkers() == 1 ? 0 : nworkers()))
@assert nworkers() == NCHAINS "need $NCHAINS workers, got $(nworkers())"
@everywhere using Distributed, HiQGA.transD_GP

include(joinpath(COMMON, "aem_lib.jl"))
using JLD2, Statistics, Random, LinearAlgebra

"""
OptionsStat for the AEM model space, in log10 resistivity.

Reproduced from gadi/src/aem_setup.jl `make_aem_opt` with only `fbounds`
changed; see the header for why it is copied rather than imported.
"""
function make_aem_opt(znall, prefix; nmin = 2, nmax = 40, fbounds = FBOUNDS)
    K = transD_GP.GP.OrstUhn()
    sdev_pos  = [0.05 * abs(diff([extrema(znall)...])[1])]
    sdev_prop = 0.07 * diff(fbounds, dims = 2)[:]
    sdev_dc   = 0.008 * diff(fbounds, dims = 2)[:]
    transD_GP.OptionsStat(
        nmin = nmin, nmax = nmax,
        xbounds = permutedims([extrema(znall)...]), fbounds = fbounds,
        fdataname = prefix, xall = permutedims(collect(znall)),
        λ = [2.0], δ = 0.1,
        sdev_prop = sdev_prop, sdev_pos = sdev_pos, sdev_dc = sdev_dc,
        save_freq = 50, demean = false, sampledc = true,
        quasimultid = false, K = K)
end

function main()
    mkpath(OUT)
    occ = joinpath(OUT, "aem_occam_$(SITE).jld2")
    isfile(occ) || error("run run_aem_occam.jl --site $SITE first (needed for the noise scale)")
    A = load(occ)
    g = load_gates(SITE)
    grid = build_grid()
    # same operator the Occam and the gate use, at the same calibrated noise
    aem = build_operator(g.d_lm, g.σ_lm, g.d_hm, g.σ_hm, A["height_m"];
                         grid = grid, calcjacobian = false,
                         σ_scale = A["noise_scale"])
    ndata = length(aem.res)
    @printf("%s: %d LM + %d HM gates, ndata %d, noise x %.3f\n",
            SITE, length(g.d_lm), length(g.d_hm), ndata, A["noise_scale"])
    @printf("  %d samples x %d chains, Tmax %.1f, log10rho bounds [%.1f, %.1f]\n",
            NSAMPLES, NCHAINS, TMAX, FBOUNDS[1], FBOUNDS[2]); flush(stdout)

    rundir = joinpath(OUT, "aem_td_$(SITE)"); mkpath(rundir)
    prefix = "udftd_$(SITE)_"
    for f in readdir(rundir)
        occursin(prefix, f) && rm(joinpath(rundir, f), force = true)
    end

    opt = make_aem_opt(grid.znall, prefix)
    old = pwd(); cd(rundir)
    t0 = time()
    try
        transD_GP.main(opt, aem; Tmax = TMAX, nsamples = NSAMPLES,
                       nchains = NCHAINS, nchainsatone = 1)
    finally
        cd(old)
    end
    @printf("  MCMC done in %.1f min\n", (time() - t0) / 60); flush(stdout)

    cd(rundir)
    opt.xall[:] .= grid.znall
    F = transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
    cd(old)
    n = length(F)
    n == 0 && error("no post-burn-in draws assembled")
    log10ρ = hcat([f[:] for f in F]...)
    σ = 10 .^ (-log10ρ)                     # S/m per AEM cell, one column per draw

    pct(p) = [quantile(view(σ, i, :), p) for i in axes(σ, 1)]
    σ10, σ50, σ90 = pct(0.10), pct(0.50), pct(0.90)
    Random.seed!(SEED)
    ix = n >= NDRAWS ? randperm(n)[1:NDRAWS] : collect(1:n)

    # split-half drift on the unit-mean conductivity: the project's standing
    # convergence check, reported rather than gated
    h = n ÷ 2
    gm(S, k) = 10 .^ mean(log10.(S[k, :]), dims = 1) |> vec
    k60 = findall(z -> 20 <= z <= 80, grid.zall)
    d1, d2 = median(gm(σ[:, 1:h], k60)), median(gm(σ[:, h+1:2h], k60))
    @printf("  %d draws; 20-80 m geometric sigma: first half %.5f, second %.5f (drift %.2fx)\n",
            n, d1, d2, max(d1, d2) / min(d1, d2))
    @printf("  p10/p50/p90 over 20-80 m: %.5f / %.5f / %.5f S/m\n",
            10^mean(log10.(σ10[k60])), 10^mean(log10.(σ50[k60])), 10^mean(log10.(σ90[k60])))
    @printf("  Occam smooth model there: %.5f S/m\n",
            10^mean(log10.(A["sigma"][k60])))

    jldsave(joinpath(OUT, "aem_td_$(SITE).jld2");
        field_id = SITE, zall = grid.zall, zbnd = grid.zbnd,
        sigma10 = σ10, sigma50 = σ50, sigma90 = σ90,
        sigma_draws = σ[:, ix], draw_ix = ix, ndraws_total = n,
        splithalf_first = d1, splithalf_second = d2,
        nsamples = NSAMPLES, nchains = NCHAINS, Tmax = TMAX,
        noise_scale = A["noise_scale"], ndata = ndata,
        fbounds = FBOUNDS, seed = SEED, secs = time() - t0)
    println("  wrote aem_td_$(SITE).jld2 ($(length(ix)) intact profiles retained)")
end

main()
