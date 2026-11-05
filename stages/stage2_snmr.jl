# run_conditions.jl — §4: the four conditions, same data.
#
#   C1   uniform nucleus values (the EK pre-run prior)   kernel: AEM p50
#   C2   corrected low-favouring prior (item C: log)     kernel: AEM p50
#   C3   = C2's posterior; σ enters only the PREDICATE   (no new inversion)
#   C3b  = C2's posterior; σ bracketed p10/p90           (no new inversion)
#   C4   corrected prior, kernel rebuilt per AEM draw    30 nested inversions
#
# C3 and C3b deliberately run no sampler: they differ from C2 only in what the
# predicate is told about σ, so re-inverting would be the same chain with a
# different seed and would confuse a predicate effect with a sampling one.
#
# Production settings throughout, matching the EK pre-run: λ_GP 3 cells,
# δ 1e-4, bounds [0.001, 0.45] (log10 of those under the corrected prior),
# hierarchical s centred on the per-q σ̂ level with τ 0.5, 4 independent
# ensembles, pooled. nsamples is 100k — the pre-run's actual figure, not the
# spec's 150k, which it was compute-limited out of.
#
# SHARDING. Each ensemble uses 4 workers (one per PT chain), so a single
# process leaves most of a 12-core box idle. `--shard i/N` splits the work
# across N concurrent processes, the same way the EK pre-run sharded soundings.
# For C1/C2 the split is over conditions; for C4 it is over the 30 AEM draws.
#
#   julia --project=. synthetics/run_conditions.jl --tier typical --conds C1,C2
#   julia --project=. synthetics/run_conditions.jl --conds C4 --shard 1/3

using Distributed, Printf

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        if startswith(ARGS[i], "--") && i < length(ARGS)
            ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2
        else; i += 1; end
    end
end
const TIER   = get(ARG, "tier", "typical")
# --tag 2 switches to the two-aquifer column and its Occam-derived SMR data,
# and writes cond_<C>_<tier>2.jld2 so the two experiments never collide.
const TAG    = get(ARG, "tag", "")
# The "T" suffix selects the trans-D SMR data and output, but the TRUTH column
# is the same two-aquifer one either way — truth_murray2.jld2, not 2T.
const TTAG   = replace(TAG, "T" => "")
const CONDS  = split(get(ARG, "conds", "C1,C2"), ",")
const NSAMP  = parse(Int, get(ARG, "nsamples", "100001"))
const NENS   = 4
const NCHAIN = 4          # PT chains per ensemble
# §2 saves 30 AEM draws; C4 uses the first NDRAWS of them. Cut to 10 on
# 2026-10-12 at the user's direction to halve the wall clock. This is the outer
# Monte Carlo sample size for ∫p(θ|d,σ)p(σ|d_AEM)dσ, so a smaller N makes the
# C4 width noisier, not biased — eval_hypotheses bootstraps over the draws so
# prediction 4 can still be scored with an interval rather than a point.
const NDRAWS = parse(Int, get(ARG, "ndraws", "30"))
const SHARD, NSHARD = let t = split(get(ARG, "shard", "1/1"), "/")
    parse(Int, t[1]), parse(Int, t[2])
end

addprocs(NCHAIN - (nworkers() == 1 ? 0 : nworkers()))
@everywhere using Statistics, Random, Printf
@everywhere using HiQGA.transD_GP, SMRPInversion, SNMRForward
using Statistics, Random, JLD2, LinearAlgebra

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(HERE, "..", "data")
const RUN  = joinpath(OUT, "smr_runs"); mkpath(RUN)
# Vendored into src/. This include used to point at ../ek_prerun/, which is
# outside the repository and was never shipped, so Stage 2 -- and therefore
# every stage after it -- could not start from a clone at all.
include(joinpath(ROOT, "src", "prerun_config.jl")); const CFG = PrerunConfig

const FB_LINEAR = [0.001 0.45]
const FB_LOG    = [log10(0.001) log10(0.45)]
const SBOUNDS   = [log(1e-2) log(1e4)]

"Production OptionsStat. `fb` is the ACTIVE bounds, so the step sizes scale with
 the parameterisation rather than being copied from the linear one."
function make_opt(prefix, fb, nz, zbnd)
    opt = transD_GP.OptionsStat(
        nmin = CFG.NMIN, nmax = CFG.NMAX,
        xbounds = [1.0 Float64(nz)], fbounds = fb, fdataname = prefix,
        xall = permutedims(collect(1:1.0:nz)), λ = CFG.LAMBDA_GP, δ = CFG.DELTA,
        sdev_prop = 0.07 * diff(fb, dims = 2)[:], sdev_pos = [0.05 * nz],
        sdev_dc = 0.008 * diff(fb, dims = 2)[:], save_freq = CFG.SAVE_FREQ,
        quasimultid = false, K = CFG.GP_KERNEL)
    opt.xall[:] .= zbnd
    opt
end

clear_tag(dir, p) = for f in readdir(dir); occursin(p, f) && rm(joinpath(dir, f), force = true); end

"""
One ensemble. Returns its post-burn-in T=1 draws in LINEAR water content —
the log parameterisation is undone here, before any statistic is taken, so no
downstream code has to know which prior produced the samples.
"""
function run_ensemble(tag, e, V, σ̂, F, islog, zbnd, pids)
    prefix = "$(tag)_e$(e)_"
    nz = length(zbnd); fb = islog ? FB_LOG : FB_LINEAR
    old = pwd(); cd(RUN)
    local samples, secs
    try
        clear_tag(pwd(), prefix)
        Random.seed!(hash((tag, e)) % 10^8)
        S = SMRPInversion.newSMRSoundingHierarchical(V, zeros(length(V)), F;
                σ̂_V0 = σ̂, s_prior_center = CFG.S_PRIOR_CENTER,
                s_prior_logwidth = CFG.S_PRIOR_TAU, s_init = CFG.S_INIT,
                linearsat = !islog, amponly = true)
        opt  = make_opt(prefix, fb, nz, zbnd)
        optn = transD_GP.OptionsNuisance(opt; sdev = copy(CFG.S_SDEV),
                bounds = copy(SBOUNDS), nsamples = NSAMP, updatenuisances = true)
        optn.fdataname = prefix
        secs = @elapsed transD_GP.main(opt, optn, S, pids;
                    nsamples = NSAMP, nchainsatone = 1, Tmax = CFG.TMAX)
        opt2 = make_opt(prefix, fb, nz, zbnd)
        opt2.fstar_filename    = "models_"  * prefix * "s.bin"
        opt2.costs_filename    = "misfits_" * prefix * "s.bin"
        opt2.x_ftrain_filename = "points_"  * prefix * "s.bin"
        Fs = transD_GP.assembleTat1(opt2, :fstar; temperaturenum = 1,
                                    burninfrac = CFG.BURNINFRAC)
        samples = isempty(Fs) ? zeros(nz, 0) : hcat([f[:] for f in Fs]...)
        islog && (samples = 10.0 .^ samples)     # back to θ BEFORE any statistic
        for f in readdir(pwd()); startswith(f, "points_") && rm(f, force = true); end
    finally
        cd(old)
    end
    Float32.(samples), secs
end

"All four ensembles for one condition, kept separate so the spread is reportable."
function run_condition(tag, V, σ̂, F, islog, zbnd, pids)
    ens = Vector{Matrix{Float32}}(undef, NENS); tot = 0.0
    for e in 1:NENS
        s, secs = run_ensemble(tag, e, V, σ̂, F, islog, zbnd, pids)
        ens[e] = s; tot += secs
        @printf("    %s e%d: %d draws, %.1f min\n", tag, e, size(s, 2), secs/60); flush(stdout)
    end
    ens, tot
end

"Kernel for a given σ(z) on the AEM layering."
function kernel_for(σ_aem, z_aem, zsmr, R, ϕ, Be)
    ct = SNMRForward.ConductivityModel(σ_aem, diff(z_aem))
    SNMRForward.MRSForward_square(R, zsmr, QGRID, ϕ, 0.0, Be, ct)
end

function main()
    T = load(joinpath(OUT, "truth_murray$(TTAG).jld2"))
    dpath = joinpath(OUT, "smr_data_murray$(TAG).jld2")
    isfile(dpath) || error("""
    no SNMR data at $dpath

      Forward-model it first. For a trans-D Stage 1 under tag <T> and label <L>:
        julia --project=. tools/aem_HYP.jl 200001 <T> 4 <L>
        julia --project=. tools/smr_data.jl transd --tag <T> --label <L>
      That writes smr_data_murray<T>T.jld2, so run this stage with --tag <T>T.

      This file used to have no producer in the release at all, which stopped an
      independent review dead after Stage 1.
    """)
    D = load(dpath)
    # Kernel source, by tag. This was a hardcoded three-way branch in which every
    # tag other than "" and "2T" fell through to aem_occam.jld2 — the TWO-AQUIFER
    # Occam profile. Running tag 3T therefore forward-modelled the SMR data with
    # the resistive column's conductivity and inverted it with the Murray
    # column's, which is far more conductive, so the inversion made up the missing
    # signal with water and railed theta at its 0.45 prior ceiling. Resolved
    # generically: a "T" suffix means the trans-D p50 of THAT column, no suffix
    # means that column's Occam profile.
    kernelfile = TAG == "" ? joinpath(OUT, "aem", "aem_posterior.jld2") :
                 endswith(TAG, "T") ? joinpath(OUT, "aem$(TTAG)", "aem_posterior.jld2") :
                 TTAG == "2" ? joinpath(OUT, "aem_occam.jld2") :
                               joinpath(OUT, "aem_occam$(TTAG).jld2")
    isfile(kernelfile) || error("kernel source not found: $kernelfile")
    A = load(kernelfile)
    global QGRID = D["qgrid"]

    zsmr = D["zsmr"]; R = D["R"]; Be = D["Be"]; ϕ = deg2rad(D["inc_deg"])
    V  = TIER == "poor" ? D["poor_V"] : D["typical_V"]
    σ̂  = TIER == "poor" ? D["poor_σ̂"] : D["typical_σ̂"]
    zA = A["zall_aem"]; keep = findall(z -> z > 0, zA); zA = zA[keep]
    σref = haskey(A, "σ50") ? A["σ50"][keep] : A["σ_occam"][keep]
    @printf("SMR kernel built from %s  (σ %.4f–%.4f S/m)\n",
            basename(dirname(kernelfile)) == "outputs" || isempty(basename(dirname(kernelfile))) ?
            basename(kernelfile) : joinpath(basename(dirname(kernelfile)), basename(kernelfile)),
            extrema(σref)...)

    pids = workers()
    @printf("run_conditions — tier %s, conds %s, %d samples × %d chains × %d ensembles\n",
            TIER, join(CONDS, ","), NSAMP, NCHAIN, NENS)
    @printf("workers %s; σ̂ median %.2f nV\n\n", string(pids), 1e9*median(σ̂))

    F50 = kernel_for(σref, zA, zsmr, R, ϕ, Be)

    todo = [c for (i, c) in enumerate(CONDS) if c == "C4" || mod(i - 1, NSHARD) == SHARD - 1]
    for cond in todo
        f = joinpath(OUT, "cond_$(cond)_$(TIER)$(TAG).jld2")
        cond == "C4" && (f = joinpath(OUT, "cond_C4_$(TIER)_shard$(SHARD).jld2"))
        if isfile(f); @printf("%s already done, skipping\n", cond); continue; end
        t0 = time()
        if cond in ("C1", "C2")
            islog = cond == "C2"
            # TAG belongs in the chain-file prefix: without it every run wrote
            # points_/models_/misfits_C2_typical_e1_* into the shared smr_runs
            # directory, so two tags running at once overwrote each other's
            # chains. Harmless sequentially, fatal for a batch.
            @printf("── %s  prior=%s  kernel=AEM p50\n", cond, islog ? "log (corrected)" : "uniform")
            ens, _ = run_condition("$(cond)_$(TIER)$(TAG)", V, σ̂, F50, islog, zsmr, pids)
            jldsave(f; ens = ens, cond = cond, tier = TIER, islog = islog,
                    kernel = basename(kernelfile),
                    nsamples = NSAMP, zsmr = zsmr)
        elseif cond == "C4" && TAG == "2"
            @warn "C4 needs AEM posterior draws; the Occam run has none. Skipping."
            continue
        elseif cond == "C4"
            @printf("── C4  prior=log (corrected)  kernel=per AEM draw (nested)\n")
            draws = A["σ_draws"][keep, :]
            nd = min(NDRAWS, size(draws, 2))
            draws = draws[:, 1:nd]
            mine = [k for k in 1:nd if mod(k - 1, NSHARD) == SHARD - 1]
            @printf("   shard %d/%d takes draws %s\n", SHARD, NSHARD, string(mine))
            allens = Dict{Int,Any}()
            for k in mine
                Fk = kernel_for(draws[:, k], zA, zsmr, R, ϕ, Be)
                e, _ = run_condition("C4_$(TIER)$(TAG)_d$(k)", V, σ̂, Fk, true, zsmr, pids)
                allens[k] = e
                @printf("  C4 draw %d (shard %d/%d) done, %.1f min elapsed\n",
                        k, SHARD, NSHARD, (time()-t0)/60)
                flush(stdout)
                jldsave(joinpath(OUT, "cond_C4_$(TIER)_shard$(SHARD).jld2");
                        ens_by_draw = allens, cond = "C4", tier = TIER,
                        kernel = "AEM_draws", nsamples = NSAMP, zsmr = zsmr,
                        σ_draws = draws, z_aem = zA, shard = SHARD, nshard = NSHARD)
            end
        else
            @warn "C3/C3b run no sampler — they reuse C2. Nothing to do." cond
            continue
        end
        @printf("   %s complete in %.1f min -> %s\n\n", cond, (time()-t0)/60, basename(f))
    end
    println("done")
end

main()
