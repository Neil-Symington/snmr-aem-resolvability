# Ported from udf_realdata/ into the repository. The include paths and the
# project are the only things changed from the working version: COMMON is
# fieldrun/common/ and the project is the repository root, not ../ek_prerun.
# run_smr_transd.jl — trans-dimensional SNMR inversion of the four chosen UDF
# soundings, with intact posterior draws retained.
#
# Mirrors ek_realdata/run_smr_transd.jl. The site/job assembly lives in
# udf_smr.jl because UDF has no EK-style inversion spreadsheet and no supplied
# resistivity profile, so `--supplied` is not offered.
#
#   julia --project=. -p 4 run_smr_transd.jl --sites A,B,C
#
# Three deliberate differences from ek_prerun/run_prerun_mcmc.jl:
#
#  1. THE CONDUCTIVE EARTH COMES FROM THE AEM INVERSION at the nearest fiducial,
#     not from the operator-supplied resistivity profile. The SMR kernel is
#     conductivity-sensitive, so this is the step that actually makes the two
#     datasets joint rather than merely co-located. The supplied profile is
#     still loaded and its conductance recorded, for comparison.
#
#  2. s PRIOR WIDTH tau = 2.0, not the pre-run's 0.5. Measured (ek_prerun
#     noise_floor.jl): the processed GMR records carry no recoverable noise, so
#     sigma-hat from est_E0_T2 is decay-model mismatch and supplies a per-moment
#     SHAPE but never a level. At tau = 0.5 the EK posterior s sat at chi2/N
#     2.90 against SSC's 1.02 -- the prior was too narrow by ~1.7x in sigma.
#
#  3. INTACT DRAWS ARE SAVED, not just q05/q50/q95. The report's detection
#     criteria are integrals over the posterior and its water-quality numbers
#     are propagated through a non-linear rock-physics transform; percentile
#     profiles are not posterior draws and cannot be used for either.

using Distributed, Printf

# The control file first: support modules read it at load time.
include(joinpath(@__DIR__, "..", "common", "stage_prelude.jl"))

# --flags are parsed in common/smr_setup.jl, included below.

# PRERUN is gone: prerun_data.jl and prerun_lib.jl are in common/

# THE CHAINS NEED PROCESSES, and this stage never asked for any. run_batch.sh
# launches `julia --project=.` with no -p, so nworkers() was 1; the old
# `pids = workers()[1:min(NCHAINS, nworkers())]` then quietly ran ONE chain
# instead of the config's nchains, and HiQGA rejected it one layer down:
#
#   ArgumentError: range(0.0, stop=0.3979400086720376, length=1): endpoints differ
#
# a log-spaced ladder of one chain between T=1 and T=2.5 has no endpoints.
# Both sites failed that way in the first full batch while the stage still
# exited 0, so the runner printed "snmr_transd ok" and the three sweep stages
# went on to skip every site for missing inputs. Ask for the processes the
# control file specifies and assert, as run_aem_transd.jl does -- a chain
# count that silently disagrees with [mcmc] is not a run worth saving.
addprocs(NCHAINS_CFG - (nworkers() == 1 ? 0 : nworkers()))
@assert nworkers() == NCHAINS_CFG "need $(NCHAINS_CFG) workers, got $(nworkers())"

# THE WORKERS NEED THE CONFIG TOO, and this is where the port bit back.
# `@everywhere` runs on every process, and COMMON is defined by the prelude on
# the MASTER only -- so this failed with `UndefVarError: COMMON not defined`
# inside a nested task, on both sites. And prerun_data.jl now reads
# FieldConfig.CFG[] at load time for its East Kimberley paths, so a worker that
# had COMMON would still have failed one line later with no config loaded.
#
# Both values are interpolated from the master with $, so each worker gets the
# literal strings and bootstraps itself: fieldconfig, then the config, then
# prerun_data. Lifting paths into a control file means every process that loads
# a module needs that file, not just the one that parsed the arguments.
@everywhere begin
    const COMMON = $(COMMON)
    include(joinpath(COMMON, "fieldconfig.jl"))
    using .FieldConfig
    FieldConfig.load_config($(CFGPATH))
    include(joinpath(COMMON, "prerun_data.jl"))
end
include(joinpath(COMMON, "prerun_lib.jl"))
include(joinpath(COMMON, "udf_smr.jl"))
using .PrerunLib
using HiQGA.transD_GP, JLD2, Statistics, Random, SNMRForward, DelimitedFiles

const HERE     = @__DIR__
# OUT comes from the prelude, which resolves [paths] outputs. Each stage and
# aem_lib.jl used to define its own `joinpath(HERE, "outputs")`, and because
# the support modules are included AFTER the prelude theirs won -- so the
# first real run wrote to fieldrun/common/outputs/ instead of where the
# control file said. Four separate definitions of one path, which is what
# the control file exists to prevent.
# RUNDIR, NSAMPLES, NENS, NCHAINS, the [noise] constants, FLAT/OCCAMK/LABEL,
# make_opts and hier live in common/smr_setup.jl, so report_smr.jl rebuilds
# exactly the options and sounding this stage sampled with.
include(joinpath(COMMON, "smr_setup.jl"))

mkpath(RUNDIR)


function clear_tag(dir, prefix)
    for f in readdir(dir)
        occursin(prefix, f) && rm(joinpath(dir, f), force = true)
    end
end

"""
Swap in the AEM conductivity model, keeping everything else about the sounding.

`ens` selects WHICH conductivity the kernel is built on. The paper's Stage 2
evaluates the kernel "on a conductivity draw from Stage 1 for each SNMR draw",
so that kernel uncertainty is propagated. transD_GP fixes the operator for a
whole chain, so one kernel per sample is not reachable without re-architecting
the sampler. What is reachable, and is done here, is one conductivity draw PER
ENSEMBLE: the pooled posterior then marginalises over `NENS` kernel
realisations instead of conditioning on one. That is a coarse version of the
paper's Stage 2, and it is labelled as such -- `n_kernel_draws` is stored
beside the result so the approximation is never silent.

Falls back to the smooth Occam model when no trans-D ensemble exists.
"""
function with_aem_conductivity(b, job, field_id; ens::Int = 0)
    td = joinpath(OUT, "aem_td_$(field_id).jld2")
    f = joinpath(OUT, "aem_occam_$(field_id).jld2")
    local σ, ztop, src, nk
    if isfile(td) && !OCCAMK
        A = load(td); D = A["sigma_draws"]
        nk = size(D, 2)
        σ = D[:, mod1(max(ens, 1), nk)]; ztop = A["zbnd"]
        src = "trans-D AEM draw $(mod1(max(ens,1),nk)) of $nk"
    else
        isfile(f) || error("no AEM inversion for $field_id — run run_aem_occam.jl first")
        A = load(f); σ = A["sigma"]; ztop = A["zbnd"]; nk = 1
        src = "AEM Occam smooth model (no trans-D ensemble present)"
    end
    # ConductivityModel takes n layer conductivities and n-1 thicknesses, the
    # last layer being a halfspace -- the same convention as the supplied
    # resistivity profiles (prerun_data.jl load_conductivity).
    thick = diff(ztop)
    @assert length(thick) == length(σ) - 1
    cm = SNMRForward.ConductivityModel(σ, thick)
    F = SMRPInversion.MRSForward_square(job.loop_side / 2, b.zbnd, b.q,
            deg2rad(job.inc_deg), 0.0, b.Be, cm)
    merge(b, (F = F, cond = cm, cond_source = src,
              aem_sigma = σ, aem_ztop = ztop, n_kernel_draws = nk))
end


# ── convergence gates and summaries ─────────────────────────────────────────
# These three live inside ek_prerun/run_prerun_mcmc.jl rather than in
# PrerunLib, so they are reproduced here verbatim (the library functions they
# call -- rhat_and_ess, stationary_trend, doi_from_width, detection_flag,
# kl_against_uniform -- are imported, not copied).

gate_idx(zbnd) = [argmin(abs.(zbnd .- z)) for z in CFG.GATE_DEPTHS]

function apply_gates(zbnd, ens)
    gi = gate_idx(zbnd)
    dz = diff(vcat(0.0, collect(Float64, zbnd)))
    i100 = findlast(zbnd .<= 100.0)
    scal = Dict{String,Vector{Vector{Float64}}}()
    for (k, idx) in zip(["theta_$(Int(z))m" for z in CFG.GATE_DEPTHS], gi)
        scal[k] = [collect(Float64, view(e.fstar, idx, :)) for e in ens if size(e.fstar, 2) > 0]
    end
    scal["totalwater_100m"] = [vec(sum(e.fstar[1:i100, :] .* dz[1:i100], dims = 1))
                               for e in ens if size(e.fstar, 2) > 0]
    scal["log_s"]  = [e.logs for e in ens if !isempty(e.logs)]
    scal["misfit"] = [e.U    for e in ens if !isempty(e.U)]
    scal["k"]      = [Float64.(e.nodes) for e in ens if !isempty(e.nodes)]
    pw = CFG.FBOUNDS[2] - CFG.FBOUNDS[1]
    res = Dict{String,NamedTuple}()
    for (k, ch) in scal
        pwk = startswith(k, "theta_") ? pw : nothing
        res[k] = length(ch) >= 2 ? rhat_and_ess(ch; priorwidth = pwk) :
                 (rhat = NaN, ess_bulk = NaN, ess_tail = NaN, determined = false)
    end
    rh = filter(isfinite, [v.rhat for v in values(res)])
    eb = filter(isfinite, [v.ess_bulk for v in values(res)])
    et = filter(isfinite, [v.ess_tail for v in values(res)])
    st = [stationary_trend(e.U; nsig = CFG.TREND_SIGMA) for e in ens if !isempty(e.U)]
    pooled = hcat([e.fstar for e in ens if size(e.fstar, 2) > 0]...)
    shmax = NaN
    if size(pooled, 2) >= 8
        h = size(pooled, 2) ÷ 2
        shmax = maximum(abs(median(pooled[i, 1:h]) - median(pooled[i, h+1:2h])) for i in gi)
    end
    (per_scalar = res, scalars = scal,
     ndetermined = count(v -> get(v, :determined, false), values(res)),
     rhat_max = isempty(rh) ? NaN : maximum(rh),
     ess_bulk_min = isempty(eb) ? NaN : minimum(eb),
     ess_tail_min = isempty(et) ? NaN : minimum(et),
     stationary = isempty(st) ? false : all(s -> s.stationary, st),
     misfit_t = isempty(st) ? NaN : maximum(abs(s.t) for s in st),
     splithalf_max = shmax, pooled = pooled)
end

function summarise(b, pooled)
    n = size(pooled, 2)
    n == 0 && return nothing
    qs(p) = [quantile(view(pooled, i, :), p) for i in axes(pooled, 1)]
    q05, q50, q95 = qs(0.05), qs(0.50), qs(0.95)
    pw = CFG.FBOUNDS[2] - CFG.FBOUNDS[1]
    doi  = doi_from_width(b.zbnd, q05, q95, pw; widthfrac = CFG.DOI_WIDTHFRAC)
    doiA = doi_from_width(b.zbnd, q05, q95, pw; widthfrac = CFG.DOI_WIDTHFRAC_ALT[1])
    doiB = doi_from_width(b.zbnd, q05, q95, pw; widthfrac = CFG.DOI_WIDTHFRAC_ALT[2])
    (q05 = q05, q50 = q50, q95 = q95, ndraws = n,
     doi = doi, doi_05 = doiA, doi_07 = doiB)
end

function run_ensemble(b, sounding_id, e, pids)
    prefix = "$(sounding_id)$(LABEL)_e$(e)_"
    clear_tag(pwd(), prefix)
    sd = seed_for(sounding_id * LABEL, e)
    @sync for p in pids
        @async remotecall_fetch(() -> (Random.seed!(sd + myid()); nothing), p)
    end
    Random.seed!(sd)
    S = hier(b)
    opt, optn = make_opts(prefix, b.zbnd)
    t = @elapsed transD_GP.main(opt, optn, S, pids;
                                nsamples = NSAMPLES, nchainsatone = 1,
                                Tmax = CFG.TMAX)
    set_run_filenames!(opt, optn)
    inv = assert_one_chain_at_T1(opt.costs_filename, length(pids))
    opt.xall[:] .= b.zbnd
    F = transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1,
                               burninfrac = CFG.BURNINFRAC)
    fstar = isempty(F) ? zeros(length(b.zbnd), 0) : hcat([f[:] for f in F]...)
    nu = try
        transD_GP.assemblenuisancesatT(optn; temperaturenum = 1,
                                       burninfrac = CFG.BURNINFRAC)
    catch e2
        @warn "nuisance assembly failed for $prefix" exception = e2
        zeros(0, 1)
    end
    logs = size(nu, 1) == 0 ? Float64[] : vec(nu[:, end])
    h = nodes_at_T1(opt.costs_filename; burninfrac = CFG.BURNINFRAC)
    @printf("      ensemble %d: %d draws, %.1f s, T=1 %s\n",
            e, size(fstar, 2), t, inv.ok ? "ok" : inv.why)
    (prefix = prefix, seed = sd, secs = t, fstar = fstar, logs = logs,
     nodes = h.nodes, U = h.U, t1_ok = inv.ok)
end

function run_site(field_id, filename, pids)
    r = udf_meta(field_id)
    path, fname = udf_matfile(field_id)
    job = merge(r, (filename = fname, matpath = path,
                    emi = occursin("noemi", fname) ? :noemi : :raw,
                    sounding_id = "$(field_id)__$(replace(fname, ".mat" => ""))"))

    b0 = build_udf_sounding(job)
    b0.ok || error("$field_id: $(b0.why)")
    # ens = 1 only for the printed summary; each ensemble rebuilds its own below
    b = udf_kernel(b0, job, field_id, OUT; ens = 1, occam = OCCAMK)

    @printf("\n%s  (%s, %.0f ms, %.0f m loop, inc %.1f deg)\n", field_id,
            job.filename, job.pulse_ms, job.loop_side, job.inc_deg)
    @printf("    nq %d, Be %.1f nT (f %.0f Hz), sigma-hat median %.2f nV\n",
            b.nq, b.Be * 1e9, b.freq, median(b.σ̂) * 1e9)
    @printf("    conductivity: %s\n", b.cond_source)

    dir = joinpath(RUNDIR, field_id * LABEL); mkpath(dir)
    old = pwd(); cd(dir)
    ens = Vector{Any}(undef, NENS)
    try
        for e in 1:NENS
            # Stage 2: this ensemble's kernel is built on its own Stage-1 draw.
            be = udf_kernel(b0, job, field_id, OUT; ens = e, occam = OCCAMK)
            ens[e] = run_ensemble(be, job.sounding_id, e, pids)
        end
    finally
        cd(old)
    end

    G = apply_gates(b.zbnd, ens)
    S = summarise(b, G.pooled)
    @printf("    pooled %d draws; rhat_max %.3f, ess_bulk_min %.0f, gates %s\n",
            S.ndraws, G.rhat_max, G.ess_bulk_min,
            (G.rhat_max <= CFG.GATE_RHAT && G.ess_bulk_min >= CFG.GATE_MIN_ESS &&
             G.stationary) ? "PASS" : "FAIL")
    s_draws = vcat([e.logs for e in ens if !isempty(e.logs)]...)
    @printf("    noise scale s: p05 %.2f p50 %.2f p95 %.2f  (x sigma-hat)\n",
            quantile(exp.(s_draws), .05), quantile(exp.(s_draws), .5),
            quantile(exp.(s_draws), .95))

    jldsave(joinpath(OUT, "smr_transd_$(field_id)$(LABEL).jld2");
        field_id = field_id, sounding_id = job.sounding_id,
        filename = job.filename, emi = String(string(job.emi)),
        pulse_ms = job.pulse_ms, loop_side_m = job.loop_side,
        inc_deg = job.inc_deg, nq = b.nq, q = b.q, V0 = b.V0, phi = b.ϕ,
        sigma_hat = b.σ̂, Be = b.Be, freq = b.freq,
        zbnd = collect(Float64, b.zbnd), znall = collect(Float64, b.znall),
        theta = G.pooled,              # intact draws, ncells x ndraws
        log_s = s_draws, nodes = vcat([Float64.(e.nodes) for e in ens]...),
        U = vcat([e.U for e in ens]...),
        q05 = S.q05, q50 = S.q50, q95 = S.q95, ndraws = S.ndraws,
        doi_width = S.doi, doi_05 = S.doi_05, doi_07 = S.doi_07,
        rhat_max = G.rhat_max, ess_bulk_min = G.ess_bulk_min,
        ess_tail_min = G.ess_tail_min, stationary = G.stationary,
        splithalf_max = G.splithalf_max,
        gates_pass = (G.rhat_max <= CFG.GATE_RHAT &&
                      G.ess_bulk_min >= CFG.GATE_MIN_ESS && G.stationary),
        cond_source = b.cond_source, n_kernel_draws = b.n_kernel_draws,
        nsamples = NSAMPLES, nensembles = NENS, nchains = NCHAINS, s_tau = S_TAU,
        flat_sigmahat = FLAT, supplied_conductivity = SUPPLIED,
        seeds = [e.seed for e in ens], secs = sum(e.secs for e in ens),
        fbounds = CFG.FBOUNDS)
    println("    wrote smr_transd_$(field_id)$(LABEL).jld2")
end

function main()
    mkpath(OUT)
    # --site (one site) or --sites (comma-separated). ONLY comes from the
    # prelude and is what run_batch.sh passes, as --site, once per site; this
    # stage read only --sites and so IGNORED it, running every site on every
    # invocation. With npar=2 that meant two processes each doing both sites
    # and racing to write the same two output files.
    sites = if !isempty(ONLY)
        [ONLY]
    elseif haskey(ARG, "sites")
        String.(split(ARG["sites"], ","))
    else
        String[]
    end
    # site -> chosen SNMR file, from the selection step
    sel = Dict{String,String}()
    raw, hdr = let r = readdlm(joinpath(FieldConfig.cfgpath(FieldConfig.CFG[], "aem_sounds"), "manifest.csv"),
                               ',', String; header = true); (r[1], vec(r[2])); end
    ci = Dict(strip(h) => i for (i, h) in enumerate(hdr))
    for i in axes(raw, 1)
        sel[raw[i, ci["field_id"]]] = ""
    end
    isempty(sites) && (sites = collect(keys(sel)))
    pids = workers()   # exactly NCHAINS, asserted at load
    @printf("trans-D SNMR: %d sites, %d ensembles x %d chains, %d samples, s tau %.1f\n",
            length(sites), NENS, length(pids), NSAMPLES, S_TAU)
    failed = String[]
    for s in sites
        try
            run_site(s, sel[s], pids)
        catch e
            @warn "site $s failed" exception = (e, catch_backtrace())
            push!(failed, s)
        end
    end
    # A caught exception is still a failure. This loop used to swallow every
    # site and return normally, so the stage exited 0 having written nothing
    # and run_batch.sh reported "ok" for a run in which both sites threw.
    if !isempty(failed)
        @printf("FAILED: %d of %d sites (%s)\n", length(failed), length(sites),
                join(failed, ", "))
        exit(1)
    end
    @printf("all %d sites ok\n", length(sites))
end

main()
