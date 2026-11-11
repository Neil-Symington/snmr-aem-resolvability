# aem_HYP.jl — §2 of the hypothesis-paper synthetic: AEM forward, noise, invert.
#
# SkyTEM-class system and AusAEM-standard electronics_halt noise, both from
# gadi/src/aem_setup.jl (noise_scale = 1.0 is the AusAEM standard; not tuned).
# Trans-D inversion at 200k × 4 chains per §2.
#
# Products, both logged:
#   * p10/p50/p90 σ(z) — what AusAEM would ship, and what the SMR kernel uses
#     in conditions C1–C3b (§3 takes the kernel from the p50, as the EK pre-run did)
#   * 30 intact posterior draws — the nested tier-C4 input
#
# TAG selects the truth column and the output directory, so the two-aquifer
# trans-D run cannot overwrite the earlier one or the Occam result:
#   TAG=2 -> truth_murray2.jld2, data/aem2/
#
#   julia --project=. synthetics/aem_HYP.jl [nsamples] [tag]

using Distributed, JLD2, Statistics, Random, Printf

const HERE = normpath(joinpath(@__DIR__))
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(HERE, "..", "data")
const NSAMPLES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 200_001
const TAG      = length(ARGS) >= 2 ? ARGS[2] : ""
# Optional 3rd/4th args: chain count and a run label. Defaults reproduce the
# production run exactly. A label writes to aem<TAG><LABEL>/ so a longer or
# wider-tempered chain cannot overwrite the production posterior.
const NCHAINS  = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 4
const LABEL    = length(ARGS) >= 4 ? ARGS[4] : ""
const TMAX     = 2.5
const NDRAWS   = 30
const SEED     = 20261012
# The AEM noise realisation is pinned to HiQGA's makenoisydata! default, 123,
# NOT to SEED. Until build_skytem_operator passed rseed through, the caller's
# seed was silently ignored and every run drew the 123 realisation -- which is
# therefore the realisation behind every published number. Pinning it here
# keeps those reproducible now that the keyword is honoured; vary it only to
# test sensitivity to the noise draw, never to reproduce the paper.
const NOISE_RSEED = 123

# aem_setup.jl includes HiQGA's electronics_halt by absolute path, so this works
# from anywhere — but the include must stay at top level, because that nested
# include is itself a top-level expression.
const GADI_SRC = joinpath(ROOT, "src")
include(joinpath(ROOT, "src", "aem_setup.jl"))
using HiQGA.transD_GP

# Workers must be started and loaded at TOP LEVEL: `@everywhere using ...` is a
# toplevel-only expression and fails inside a function.
addprocs(NCHAINS - (nworkers() == 1 ? 0 : nworkers()))
@assert nworkers() == NCHAINS "need $NCHAINS workers, got $(nworkers())"
@everywhere using Distributed, HiQGA.transD_GP

function main()
    mkpath(OUT)
    T = load(joinpath(OUT, "truth_murray$(TAG).jld2"))
    # The truth runs to 200 m on a 0.25 m grid; pass it straight in so the AEM
    # sees the full column rather than the SMR grid's 147 m truncation.
    built = build_skytem_operator(T["σ_bulk"], T["z"]; noise_scale = 1.0,
                                  rseed_noise = NOISE_RSEED)
    aem, grid = built.aem, built.grid
    @printf("AEM operator: %d LM + %d HM gates; truth σ interpolated to %d AEM cells\n",
            length(aem.dlow), length(aem.dhigh), length(grid.z))

    rundir = joinpath(OUT, "aem$(TAG)$(LABEL)"); mkpath(rundir)
    # transD_GP.main asserts its history files do not exist, so a re-run must
    # clear the previous chains first (same pattern as prerun_lib.clear_tag).
    for f in readdir(rundir)
        occursin("hypsynth_aem_", f) && rm(joinpath(rundir, f), force = true)
    end
    jldsave(joinpath(rundir, "aem_setup.jld2");
        σ_aem_truth = built.σ_aem_truth, ρ_aem_truth = built.ρ_aem_truth,
        z_aem = grid.z, zboundaries_aem = grid.zbnd,
        znall_aem = grid.znall,    # index space — not a depth axis
        zall_aem  = grid.zall,     # physical depth (m)
        dlow = aem.dlow, dhigh = aem.dhigh,
        σlow = aem.σlow, σhigh = aem.σhigh,
        nsamples = NSAMPLES, nchains = NCHAINS, Tmax = TMAX,
        noise_scale = 1.0, rseed_noise = NOISE_RSEED)

    opt = make_aem_opt(grid.znall, "hypsynth_aem_")
    old = pwd(); cd(rundir)
    t0 = time()
    # THE CHAIN IS NOT SEEDED, AND CANNOT BE FROM HERE. `Random.seed!(SEED)`
    # used to live BELOW this call, next to the subsampling, where it seeds only
    # the randperm that picks which NDRAWS draws to keep -- so `seed = SEED`
    # went into every aem_posterior.jld2 describing the subsampling while the
    # MCMC ran from whatever RNG state its worker task held. Measured on the
    # field branch: two runs, one config, one seed, gave sigma_draws differing
    # by 8.05 S/m with draw_ix IDENTICAL, which is the tell that the seeded part
    # reproduced and the chain did not.
    #
    # Moving the seed here does NOT fix it, and these lines are kept only
    # because they are correct for anything the master draws. Julia's default
    # RNG is TASK-LOCAL and Distributed runs every remotecall in a fresh task,
    # so seeding a worker in one call does not reach the sampling done in the
    # next. Measured directly, on Julia 1.12.6:
    #
    #   seed-then-draw in separate remote calls reproduces: false
    #   seed-and-draw in the SAME remote call reproduces:   true
    #
    # Worse, HiQGA's MCMC_Driver.jl issues a fresh remotecall_fetch(do_mcmc_step)
    # for EVERY step, so each step is a new task and a once-per-chain seed would
    # be lost after the first one. It does no seeding of its own and main() takes
    # no seed argument, so the per-chain RNG is not reachable from a caller.
    # Fixing it needs HiQGA to carry an explicit RNG in the chain state, or to
    # reseed each step from (chain_seed, step). Until then this posterior is one
    # draw from a distribution of runs, not a reproducible result, and the
    # `seed` field below does not describe it. This mirrors the comment in
    # fieldrun/stages/run_aem_transd.jl: the defect is identical on both
    # branches and was fixed there first.
    let sd = SEED   # a local: the closure must carry the VALUE to the worker,
                    # not a reference to a Main binding the worker lacks
        @sync for p in workers()
            @async remotecall_fetch(() -> (Random.seed!(sd + myid()); nothing), p)
        end
    end
    Random.seed!(SEED)
    try
        transD_GP.main(opt, aem; Tmax = TMAX, nsamples = NSAMPLES,
                       nchains = NCHAINS, nchainsatone = 1)
    finally
        cd(old)
    end
    @printf("MCMC done in %.1f min\n", (time() - t0) / 60)

    cd(rundir)
    opt.xall[:] .= grid.znall
    F = transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
    cd(old)
    n = length(F)
    @printf("assembled %d post-burn-in draws at T=1\n", n)

    log10ρ = hcat([f[:] for f in F]...)
    σ = 10 .^ (-log10ρ)                       # S/m, per AEM cell
    pct(p) = [quantile(view(σ, i, :), p) for i in axes(σ, 1)]
    σ10, σ50, σ90 = pct(0.10), pct(0.50), pct(0.90)

    # A separate, deliberate stream for the subsampling, so which draws are
    # retained does not depend on how many proposals the chain happened to make.
    # This line used to be the ONLY seeding in the tool, which is the defect
    # above. Same offset as fieldrun/stages/run_aem_transd.jl so the two
    # branches agree.
    Random.seed!(SEED + 1_000_000)
    ix = n >= NDRAWS ? randperm(n)[1:NDRAWS] : collect(1:n)
    σ_draws = σ[:, ix]

    jldsave(joinpath(rundir, "aem_posterior.jld2");
        zall_aem = grid.zall, σ10 = σ10, σ50 = σ50, σ90 = σ90,
        σ_draws = σ_draws, draw_ix = ix, ndraws_total = n,
        # `seed` is kept under its old name so existing readers still load, but
        # it governs the SUBSAMPLING ONLY -- see the comment above the chain.
        # seed_subsample is the stream that actually produced draw_ix. Renaming
        # `seed` across every artefact is a release decision, not a fix.
        σ_aem_truth = built.σ_aem_truth, seed = SEED,
        seed_subsample = SEED + 1_000_000, chain_seeded = false)

    # Sanity: does the AEM posterior separate the fresh upper from the brackish
    # lower aquifer? That contrast is the design requirement in §1.
    zf = T["z"]; sbT = T["σ_bulk"]
    on(a,b) = mean(sbT[findall(zz -> a <= zz <= b, zf)])
    at(v,a,b) = mean(v[findall(zz -> a <= zz <= b, grid.zall)])
    println("\nrecovery by unit (p10 / p50 / p90 against truth):")
    @printf("%-32s %9s %9s %9s %9s %9s\n", "unit", "truth", "p10", "p50", "p90", "p50/truth")
    # Unit windows come from the truth file, trimmed 3 m at each end to measure
    # clear of the gradational transitions. Previously these were hardcoded to the
    # two-aquifer column, so any other TAG printed a meaningless summary.
    widths = Tuple{String,Float64}[]
    for u in T["units"]
        a, b = u[1] + 3.0, min(u[2], 115.0) - 3.0
        b - a < 2.0 && continue
        t = on(a, b)
        p10, p50, p90 = at(σ10, a, b), at(σ50, a, b), at(σ90, a, b)
        @printf("%-32s %9.4f %9.4f %9.4f %9.4f %9.2f\n", u[3], t, p10, p50, p90, p50 / t)
        push!(widths, (u[3], (p90 - p10) / p50))
    end
    println("\nfractional p10-p90 width (= (p90-p10)/p50), by unit:")
    for (nm, w) in widths
        @printf("  %-32s %6.2f\n", nm, w)
    end

    println("\nwrote ", rundir)
end

main()
