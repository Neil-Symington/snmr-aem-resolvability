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

using Distributed, Printf, Dates

const ARG = Dict{String,String}()
# Valueless flags, listed so the pair-wise loop below does not eat the token
# AFTER them. `--force --tier typical` would otherwise store force="--tier" and
# then skip "typical", silently running the default tier -- the same shape of
# defect as everything else in this file's history.
const FLAGS = ("--force",)
let i = 1
    while i <= length(ARGS)
        if ARGS[i] in FLAGS
            ARG[ARGS[i][3:end]] = "true"; i += 1
        elseif startswith(ARGS[i], "--") && i < length(ARGS)
            ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2
        else; i += 1; end
    end
end
# --force recomputes a condition whose output file already exists. Without it
# this stage SKIPPED any existing output and said so in one line, with no check
# that the existing answer came from the inputs now on disk -- so re-running
# after rebuilding Stage 1, or after changing the tier or the sample count,
# quietly kept the old posterior and every downstream number derived from it.
const FORCE = haskey(ARG, "force")
# --truekernel inverts the SAME data with the inversion kernel built on the TRUE
# conductivity, which tools/smr_data.jl stores in the data file as sigma_gen.
# It is the CONTROL: the data are generated on the truth, so with --truekernel
# the generating and inversion kernels coincide and whatever the inversion then
# loses is noise and non-uniqueness, not kernel error. Comparing it against the
# --pool run separates the two. Required from 18 Nov; nothing downstream reads
# the control, and it writes its own output name so it cannot be mistaken for one.
const TRUEKERNEL = haskey(ARG, "truekernel")
const TIER   = get(ARG, "tier", "typical")
# --tag 2 switches to the two-aquifer column and its Occam-derived SMR data,
# and writes cond_<C>_<tier>2.jld2 so the two experiments never collide.
const TAG    = get(ARG, "tag", "")
# The "T" suffix selects the trans-D SMR data and output, but the TRUTH column
# is the same two-aquifer one either way — truth_murray2.jld2, not 2T.
const TTAG   = replace(TAG, "T" => "")
# aem_HYP.jl's 4th argument writes data/aem<TAG><LABEL>/, so a labelled Stage 1
# run must be named here too. Without this the kernel lookup below misses and
# the stage dies on "kernel source not found" -- the same defect a review found
# in aem_draws.jl, which was fixed there and nowhere else. Every path that
# builds a run directory from a tag now takes the label.
const LABEL  = get(ARG, "label", "")
const CONDS  = split(get(ARG, "conds", "C1,C2"), ",")
const NSAMP  = parse(Int, get(ARG, "nsamples", "100001"))
# EIGHT ensembles from 18 Nov, confirmed by Neil. The chains are not seeded in
# practice -- Random.seed! inside remotecall_fetch is inert, HiQGA issuing a
# fresh remote call per MCMC step and Julia's RNG being task-local -- so a
# re-run is four MORE independent ensembles rather than a reproduction of the
# first four. Eight in one run is the same information with one provenance
# record instead of two. Stage 1 stays at four.
const NENS   = parse(Int, get(ARG, "nens", "8"))
# Per-ensemble seed base, from SMR_SEED or the historical default. Ensemble e
# uses SEED_BASE + e, and each of its workers SEED_BASE + e + worker id.
const SEED_BASE = haskey(ENV, "SMR_SEED") ? parse(Int, ENV["SMR_SEED"]) : 20261101
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
# EVERY RUN WRITES TO ITS OWN FOLDER: column, kernel type, timestamp. Two
# full-length results were lost to name collisions on 18 Nov -- the control's
# chain files over the pooled run's, then the pooled re-run's posterior over the
# original's -- because the output name was a function of the column alone and
# --force was passed to get past the provenance guard. A fresh folder per run
# removes the collision class, so --force is no longer needed and must not be
# passed over a full-length result.
const STAMP  = Dates.format(Dates.now(), "yyyymmdd-HHMMSS")
const KIND   = haskey(ARG, "truekernel") ? "truekernel" : "pooled"
const RUNDIR = joinpath(OUT, "stage2", "$(get(ARG, "tag", "x"))_$(KIND)_$(STAMP)")
const RUN    = joinpath(RUNDIR, "chains"); mkpath(RUN)
# Vendored into src/. This include used to point at ../ek_prerun/, which is
# outside the repository and was never shipped, so Stage 2 -- and therefore
# every stage after it -- could not start from a clone at all.
include(joinpath(ROOT, "src", "prerun_config.jl")); const CFG = PrerunConfig
# make_aem_opt, needed to re-assemble Stage 1 draws when pooling several
# folders. GADI_SRC is what aem_setup.jl expects; Stage 1 sets the same.
const GADI_SRC = joinpath(ROOT, "src")
include(joinpath(ROOT, "src", "aem_setup.jl"))

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
        # Seed PER ENSEMBLE from config or the command line, and seed every
        # worker from seed + its id. Previously this was hash((tag, e)), which
        # was deterministic but not settable, so a run could not be repeated
        # with a stated seed. SEED_BASE + e keeps ensembles independent.
        let sd = SEED_BASE + e
            Random.seed!(sd)
            @sync for p in workers()
                @async remotecall_fetch(() -> (Random.seed!(sd + myid()); nothing), p)
            end
        end
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

"""
Per-cell median of the post-burn-in T=1 draws POOLED across several Stage 1
folders.

aem_posterior.jld2 keeps only 30 subsampled draws, so the full set has to be
re-assembled from each folder's chain files. Pooling ONE folder reproduces that
folder's sigma50 exactly, because Stage 1 computes sigma50 as the same quantile
over the same draws with the same burn-in fraction -- which is the test.

burninfrac is 0.5 to match tools/aem_HYP.jl's own assembleTat1 call literally;
CFG.BURNINFRAC is the same value, but matching the call is what guarantees the
reproduction is exact.
"""
function pooled_sigma(folders)
    mats = Matrix{Float64}[]; zall = nothing; perfolder = Int[]
    for d in folders
        isdir(d) || error("pool folder not found: $d")
        S = load(joinpath(d, "aem_setup.jld2"))
        znall = S["znall_aem"]
        zall === nothing && (zall = collect(Float64, S["zall_aem"]))
        opt = make_aem_opt(znall, "hypsynth_aem_"; nmax = get(S, "nmax_k", 40))
        opt.xall[:] .= znall
        old = pwd(); cd(d)
        F = try
            transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
        finally
            cd(old)
        end
        m = 10 .^ (-hcat([f[:] for f in F]...))      # S/m per AEM cell
        push!(mats, m); push!(perfolder, size(m, 2))
        @printf("    pooled %-28s %d T=1 draws\n", basename(d), size(m, 2))
    end
    σ = hcat(mats...)
    ([quantile(view(σ, i, :), 0.50) for i in axes(σ, 1)], zall, size(σ, 2), perfolder)
end

# Moved back next to inputsig on 18 Nov. The pooled_sigma docstring added
# in a39a3e9 was inserted BETWEEN this docstring and the function it
# documents, leaving one docstring documenting another. That parses but
# fails at load with "cannot document the following expression", so Stage 2
# could not run at all. It went unnoticed because m4_run.sh skips Stage 2
# and the pooling test exercises pooled_sigma without loading this file.
"""
Everything upstream whose change would change this condition's posterior,
collapsed to one UInt64 and stored in the output as `input_sig`. It is what
lets a re-run tell "already done" apart from "done from something else",
which this stage could not do: it skipped on the mere existence of the file.

Covers the SNMR data actually used (V, sigma-hat), the kernel conductivity and
both grids, the instrument geometry, the sampler's size, and for C4 the AEM
draws. `hash` is stable within a Julia version and not across them; install.jl
pins Julia 1.12, and a version change surfaces as a mismatch, which errors with
instructions rather than silently reusing the file.
"""
inputsig(cond, islog, V, σ̂, σref, zA, zsmr, R, ϕ, Be, qgrid, kernelfile, draws) =
    hash((cond, TIER, TAG, LABEL, islog, NSAMP, NENS, NCHAIN, TRUEKERNEL,
          cond == "C4" ? NDRAWS : 0, basename(kernelfile),
          V, σ̂, σref, zA, zsmr, R, ϕ, Be, qgrid, draws))

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
                 endswith(TAG, "T") ? joinpath(OUT, "aem$(TTAG)$(LABEL)", "aem_posterior.jld2") :
                 TTAG == "2" ? joinpath(OUT, "aem_occam.jld2") :
                               joinpath(OUT, "aem_occam$(TTAG).jld2")
    # --pool supersedes this single folder entirely, so the file need not exist.
    # Added with --pool and missed: the check stayed unconditional, so pooling
    # four m4 ensembles failed on the absence of an unlabelled data/aem<TAG>/
    # that nothing was going to read. Caught by the launcher's smoke, 18 Nov.
    POOLED = haskey(ARG, "pool") && !isempty(ARG["pool"])
    !(POOLED && TRUEKERNEL) || error(
        "--pool and --truekernel are alternatives: one builds the inversion " *
        "kernel from the pooled Stage 1 median, the other from the true " *
        "conductivity. Pass one.")
    POOLED || TRUEKERNEL || isfile(kernelfile) || error("""
    kernel source not found: $kernelfile

      If Stage 1 was run with a label, pass the same one:  --label <L>
      aem_HYP.jl's 4th argument writes data/aem<TAG><LABEL>/, and this stage
      has to look in the same place.

      Or pass --pool <dir1>,<dir2>,... to build the kernel from the median of
      several Stage 1 folders' pooled draws, in which case this file is not read.
    """)
    A = (POOLED || TRUEKERNEL) ? nothing : load(kernelfile)
    # C4 rebuilds the kernel per AEM draw and reads A["σ_draws"] from the single
    # folder, so it cannot run pooled. Fail here rather than on `nothing`.
    !(POOLED && "C4" in CONDS) || error(
        "--pool cannot be combined with condition C4: C4 needs the per-draw " *
        "sigma of ONE Stage 1 folder, which pooling discards.")
    global QGRID = D["qgrid"]

    zsmr = D["zsmr"]; R = D["R"]; Be = D["Be"]; ϕ = deg2rad(D["inc_deg"])
    V  = TIER == "poor" ? D["poor_V"] : D["typical_V"]
    σ̂  = TIER == "poor" ? D["poor_σ̂"] : D["typical_σ̂"]
    # --pool takes a comma-separated list of Stage 1 folders; the kernel is then
    # the median of the POOLED post-burn-in T=1 draws. Without it, behaviour is
    # exactly as before: the single folder's own sigma50.
    POOL = POOLED ?
           [isabspath(x) ? x : joinpath(OUT, x) for x in split(ARG["pool"], ",")] :
           String[]
    local zA, σref
    global POOLED_N = -1
    if TRUEKERNEL
        haskey(D, "sigma_gen") || error("""
            --truekernel needs sigma_gen in $(basename(dpath)), and it is absent.
            That file predates the 18 Nov fix (tag snmr-data-truth-20261118).
            Regenerate it: julia --project=. tools/smr_data.jl --tag $(TTAG)""")
        zt = collect(Float64, D["z_gen"]); st = collect(Float64, D["sigma_gen"])
        keep = findall(z -> z > 0, zt); zA = zt[keep]; σref = st[keep]
        @printf("  CONTROL: inversion kernel on the TRUE conductivity (sigma_gen)\n")
    elseif isempty(POOL)
        zA = A["zall_aem"]; keep = findall(z -> z > 0, zA); zA = zA[keep]
        σref = haskey(A, "σ50") ? A["σ50"][keep] : A["σ_occam"][keep]
    else
        @printf("  POOLING %d Stage 1 folders for the kernel\n", length(POOL))
        σp, zp, ntot, perf = pooled_sigma(POOL)
        keep = findall(z -> z > 0, zp); zA = zp[keep]; σref = σp[keep]
        global POOLED_N = ntot
        @printf("  pooled %d draws across %d folders (%s)\n",
                ntot, length(POOL), join(perf, "+"))
    end
    # Name what the kernel was ACTUALLY built from. When pooling, kernelfile is
    # never read, and printing it said the kernel came from a single folder that
    # had not been opened.
    kernelsrc = TRUEKERNEL ? "the TRUE conductivity (sigma_gen)" :
        POOLED ?
        @sprintf("the pooled median of %d folders (%s)", length(POOL),
                 join(basename.(POOL), ", ")) :
        (basename(dirname(kernelfile)) == "outputs" || isempty(basename(dirname(kernelfile))) ?
         basename(kernelfile) : joinpath(basename(dirname(kernelfile)), basename(kernelfile)))
    @printf("SMR kernel built from %s  (σ %.4f–%.4f S/m)\n", kernelsrc, extrema(σref)...)

    pids = workers()
    @printf("RUN FOLDER: %s\n", RUNDIR)
    @printf("run_conditions — tier %s, conds %s, %d samples × %d chains × %d ensembles\n",
            TIER, join(CONDS, ","), NSAMP, NCHAIN, NENS)
    @printf("workers %s; σ̂ median %.2f nV\n\n", string(pids), 1e9*median(σ̂))

    F50 = kernel_for(σref, zA, zsmr, R, ϕ, Be)

    todo = [c for (i, c) in enumerate(CONDS) if c == "C4" || mod(i - 1, NSHARD) == SHARD - 1]
    for cond in todo
        f = joinpath(RUNDIR, "cond_$(cond)_$(TIER)$(TAG)$(TRUEKERNEL ? "_truekernel" : "").jld2")
        cond == "C4" && (f = joinpath(OUT, "cond_C4_$(TIER)_shard$(SHARD).jld2"))
        sig = inputsig(cond, cond == "C2", V, σ̂, σref, zA, zsmr, R, ϕ, Be,
                       QGRID, kernelfile,
                       # haskey, not a bare index: C4 against the Occam column
                       # (TAG "2") has no draws, and the branch below warns and
                       # skips it gracefully. Indexing here would throw first.
                       cond == "C4" && haskey(A, "σ_draws") ? A["σ_draws"] : nothing)
        if isfile(f) && !FORCE
            stored = try
                get(load(f), "input_sig", nothing)
            catch e
                @warn "$cond output exists but will not load; recomputing" file=f err=e
                nothing
            end
            if stored !== nothing && UInt64(stored) == sig
                @printf("%s already done from these exact inputs, skipping\n", cond)
                continue
            elseif stored !== nothing
                error("""
                $cond output exists but was built from DIFFERENT inputs.

                  file: $f

                  Its input signature is $(UInt64(stored)); the inputs on disk now
                  give $(sig). Something upstream changed: the SNMR data, the AEM
                  kernel, the tier, the sample count, the ensemble or chain count,
                  or (for C4) the number of AEM draws.

                  This stage used to skip here without looking, so a re-run after
                  rebuilding Stage 1 kept the old posterior and reported it as
                  current. Choose deliberately:

                    --force          recompute and overwrite
                    rm $f            then re-run
                """)
            else
                @warn """
                $cond output exists but carries NO input signature, so it cannot be
                checked against the inputs on disk. It predates this check.

                  file: $f

                  Skipping it, as this stage always did -- but if anything upstream
                  has been rebuilt since it was written, this is a stale answer.
                  Pass --force to recompute, or delete the file to be certain.
                """
                continue
            end
        end
        t0 = time()
        if cond in ("C1", "C2")
            islog = cond == "C2"
            # TAG belongs in the chain-file prefix: without it every run wrote
            # points_/models_/misfits_C2_typical_e1_* into the shared smr_runs
            # directory, so two tags running at once overwrote each other's
            # chains. Harmless sequentially, fatal for a batch.
            @printf("── %s  prior=%s  kernel=AEM p50\n", cond, islog ? "log (corrected)" : "uniform")
            # TRUEKERNEL belongs in the chain-file prefix for exactly the reason TAG
            # does, noted above: without it the control run wrote its points_,
            # models_, misfits_ and values_nuisance_ files over the pooled run's, and
            # the pooled run's nuisance history -- the only place the inferred noise
            # factor s lives -- was lost. Found 18 Nov, when s came out identical for
            # the pooled and control runs of the same column.
            ens, _ = run_condition("$(cond)_$(TIER)$(TAG)$(TRUEKERNEL ? "_truekernel" : "")",
                                   V, σ̂, F50, islog, zsmr, pids)
            jldsave(f; ens = ens, cond = cond, tier = TIER, islog = islog,
                    kernel = (POOLED || TRUEKERNEL) ? kernelsrc : basename(kernelfile),
                    truekernel = TRUEKERNEL,
                    # The conductivity the INVERSION kernel was built on, and its
                    # grid. Stored 18 Nov: without it, anything that wants to
                    # forward a model through the same kernel afterwards has to
                    # re-derive the pooled median by re-assembling every chain,
                    # and the kernel error between the data (true sigma) and the
                    # inversion (pooled median) cannot be measured at all.
                    sigma_kernel_used = collect(Float64, σref),
                    z_kernel_used = collect(Float64, zA),
                    nsamples = NSAMP, zsmr = zsmr, input_sig = sig,
                    pooled_folders = isempty(POOL) ? [kernelfile] : POOL,
                    pooled_ndraws = isempty(POOL) ? -1 : POOLED_N,
                    rundir = RUNDIR, stamp = STAMP, nens = NENS,
                    seed_base = SEED_BASE,
                    seeds_ensemble = [SEED_BASE + e for e in 1:NENS])
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
                        σ_draws = draws, z_aem = zA, shard = SHARD,
                        nshard = NSHARD, input_sig = sig)
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
