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
#   julia --project=. fieldrun/stages/run_smr_transd.jl --config <toml> [--site S]
#
# The stage adds its own chain workers ([mcmc] nchains_snmr), so no -p.
# With [mcmc] replicates_snmr = R > 1, each kernel draw gets R independent
# ladders, and one (draw, replicate) pair is a "unit". Units can run as
# separate processes, several to a node, and be combined afterwards:
#
#   ... --site S --ens 1,2 --rep 3     run only these units; save them, no gate
#   ... --site S --combine             gate and pool the units already on disk
#
# With neither flag the stage runs every unit and then combines them.
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
# The SNMR ladder can differ from the AEM one ([mcmc] nchains_snmr).
# --combine only gates units already on disk, so it needs no chains.
const COMBINE_ONLY = "--combine" in ARGS
COMBINE_ONLY || addprocs(NCHAINS_SNMR - (nworkers() == 1 ? 0 : nworkers()))
COMBINE_ONLY || @assert nworkers() == NCHAINS_SNMR "need $(NCHAINS_SNMR) workers, got $(nworkers())"

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

"""
The gate scalars of one ordered sample: theta at the gate depths and total
water to 100 m, from an ncells x ndraws `fstar`.
"""
function gate_scalars(fstar, zbnd)
    gi = gate_idx(zbnd)
    dz = diff(vcat(0.0, collect(Float64, zbnd)))
    i100 = findlast(zbnd .<= 100.0)
    d = Dict{String,Vector{Float64}}()
    for (z, idx) in zip(CFG.GATE_DEPTHS, gi)
        d["theta_$(Int(z))m"] = collect(Float64, view(fstar, idx, :))
    end
    d["totalwater_100m"] = size(fstar, 2) == 0 ? Float64[] :
                           vec(sum(fstar[1:i100, :] .* dz[1:i100], dims = 1))
    d
end

"""
THE GATE.

Three comparisons, of which the first available is gated on:

  replicate    ([mcmc] replicates_snmr >= 2) R independent ladders on each
               kernel draw. Their T=1 samples, each in iteration order, are
               the chains of a rank-normalised split-R-hat, per ensemble, on
               theta at the gate depths, total water to 100 m, log s, misfit
               and layer count. Independent ladders share nothing but the
               kernel, so disagreement is mixing failure and nothing else.
               R-hat is the worst over ensembles and scalars; ESS adds across
               ensembles, each contributing its worst scalar.
  within-kernel  (one replicate, >= 2 T=1 chain-id streams) the T=1 draws of
               one ladder split by chain id. Kept as a diagnostic. Once swaps
               are frequent the streams are short, state-selected visits (a
               chain swaps down when its misfit is low), and their R-hat is
               inflated: the 17 Nov Gadi run read 5.5-7.5 here while the T=1
               slot was stable to 1.05-1.33.
  across-ensemble  each ensemble has its own kernel draw, so this mixes
               sampler mixing with kernel spread. Always reported, as kernel
               spread; gated on only when nothing else is available.
"""
function apply_gates(zbnd, units)
    pw = CFG.FBOUNDS[2] - CFG.FBOUNDS[1]
    pwof(k) = startswith(k, "theta_") ? pw : nothing
    enss = sort(unique(u.e for u in units))
    byens = [filter(u -> u.e == e, units) for e in enss]
    nrep = minimum(length, byens)

    # ACROSS ensembles (kernel spread), on each ensemble's replicates pooled
    efs = [hcat([u.fstar for u in g]...) for g in byens]
    scal = Dict{String,Vector{Vector{Float64}}}()
    for f in efs
        size(f, 2) > 0 || continue
        for (k, v) in gate_scalars(f, zbnd)
            push!(get!(scal, k, Vector{Float64}[]), v)
        end
    end
    scal["log_s"]  = [vcat([u.logs for u in g]...) for g in byens if any(u -> !isempty(u.logs), g)]
    scal["misfit"] = [vcat([u.U for u in g]...) for g in byens if any(u -> !isempty(u.U), g)]
    scal["k"]      = [vcat([Float64.(u.nodes) for u in g]...) for g in byens if any(u -> !isempty(u.nodes), g)]
    res = Dict{String,NamedTuple}()
    for (k, ch) in scal
        res[k] = length(ch) >= 2 ? rhat_and_ess(ch; priorwidth = pwof(k)) :
                 (rhat = NaN, ess_bulk = NaN, ess_tail = NaN, determined = false)
    end
    rh = filter(isfinite, [v.rhat for v in values(res)])
    eb = filter(isfinite, [v.ess_bulk for v in values(res)])
    et = filter(isfinite, [v.ess_tail for v in values(res)])

    # WITHIN a kernel by chain id (precomputed per unit; diagnostic)
    cid_rh = filter(isfinite, [u.chainid_rhat for u in units])
    cid_ess = sum((u.chainid_ess for u in units if isfinite(u.chainid_ess)); init = 0.0)
    ncid = count(u -> isfinite(u.chainid_rhat), units)

    # BETWEEN replicate ladders on one kernel
    rep_by_scalar = Dict{String,Float64}()
    rep_rh = Float64[]; rep_ess = 0.0
    if nrep >= 2
        for g in byens
            ks = intersect([keys(u.slot) for u in g]...)
            r = Float64[]; q = Float64[]
            for k in ks
                v = rhat_and_ess([u.slot[k] for u in g]; priorwidth = pwof(k))
                isfinite(v.rhat) && (push!(r, v.rhat);
                    rep_by_scalar[k] = max(get(rep_by_scalar, k, -Inf), v.rhat))
                isfinite(v.ess_bulk) && push!(q, v.ess_bulk)
            end
            isempty(r) || push!(rep_rh, maximum(r))
            isempty(q) || (rep_ess += minimum(q))
        end
    end

    trips = [minimum(u.round_trips) for u in units]
    st = [stationary_trend(u.U; nsig = CFG.TREND_SIGMA) for u in units if !isempty(u.U)]
    pooled = hcat(efs...)
    gi = gate_idx(zbnd)
    shmax = NaN
    if size(pooled, 2) >= 8
        h = size(pooled, 2) ÷ 2
        shmax = maximum(abs(median(pooled[i, 1:h]) - median(pooled[i, h+1:2h])) for i in gi)
    end
    rk  = isempty(rh) ? NaN : maximum(rh)
    rc  = isempty(cid_rh) ? NaN : maximum(cid_rh)
    rr  = isempty(rep_rh) ? NaN : maximum(rep_rh)
    stat_ok = isempty(st) ? false : all(s -> s.stationary, st)
    ladder_ok = stat_ok && minimum(trips) >= 1
    basis, rm, me = if nrep >= 2 && length(rep_rh) == length(byens)
        "replicate", rr, rep_ess
    elseif ncid == length(units)
        "within-kernel", rc, cid_ess
    else
        "across-ensemble", rk, (isempty(eb) ? NaN : minimum(eb))
    end
    pass = rm <= CFG.GATE_RHAT && me >= CFG.GATE_MIN_ESS &&
           (basis == "across-ensemble" ? stat_ok : ladder_ok)
    (per_scalar = res, ndetermined = count(v -> get(v, :determined, false), values(res)),
     rhat_max = rm, rhat_mixing_max = rm, ess_mixing = me,
     rhat_replicate_max = rr, ess_replicate = nrep >= 2 ? rep_ess : NaN,
     rhat_replicate_by_scalar = rep_by_scalar, nreplicates = nrep,
     rhat_chainid_max = rc, ess_chainid = ncid == 0 ? NaN : cid_ess,
     rhat_kernel_max = rk, gate_basis = basis, gates_pass = pass,
     min_round_trips = trips,
     swap_rate = [u.swap_rate for u in units],
     frac_at_T1 = [u.frac_at_T1 for u in units],
     ess_bulk_min = me,
     ess_tail_min = isempty(et) ? NaN : minimum(et),
     stationary = stat_ok,
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

"""
The T=1 draws of ONE ladder, two ways, plus the ladder's health.

  streams  split by chain id. `transD_GP.assembleTat1` pools every T=1 row
           across chains and drops which chain made it; the split recovers it.
           Used only for the chain-id diagnostic (see apply_gates).
  slot     the T=1 slot in ITERATION order (t1_slot_order), as gate scalars.
           This is the sample a replicate comparison treats as one chain.
"""
function t1_streams(opt, optn, zbnd)
    rows = t1_rows_by_chain(opt.costs_filename; burninfrac = CFG.BURNINFRAC)
    nrows = t1_rows_by_chain(optn.costs_filename; burninfrac = CFG.BURNINFRAC, Tcol = 4)
    vals = readdlm(optn.vals_filename, ' ', String)
    vcid = parse.(Int, vals[:, 1]); vlogs = parse.(Float64, vals[:, end])
    cr = [split(l) for l in readlines(opt.costs_filename) if !isempty(strip(l))]
    ccid = [parse(Int, r[1]) for r in cr]
    Fs = Vector{Any}(undef, length(rows))
    st = map(rows) do r
        # The models file holds a chain's model ONLY while it is at T=1 (as
        # transD_GP.assembleTat1 indexes it), so the burn-in mask is taken over
        # that chain's T=1 rows, not over all its rows.
        F = transD_GP.history(opt; stat = :fstar, chain_idx = r.chain)
        Fs[r.chain] = F
        at1 = isapprox.(r.T, 1.0; atol = 1e-12)
        @assert length(F) == count(at1) "chain $(r.chain): $(length(F)) stored models, $(count(at1)) T=1 rows"
        keep = r.mask[at1]
        fs = any(keep) ? hcat([f[:] for f in F[keep]]...) : zeros(length(opt.xall), 0)
        ln = vlogs[vcid .== r.chain]
        nm = nrows[r.chain].mask
        (chain = r.chain, fstar = fs,
         logs = length(ln) == length(nm) ? ln[nm] : Float64[])
    end
    # the T=1 slot in iteration order
    ord = t1_slot_order([r.T for r in rows]; burninfrac = CFG.BURNINFRAC)
    crc = [cr[ccid .== c] for c in eachindex(rows)]
    lnc = [vlogs[vcid .== c] for c in eachindex(rows)]
    fslot = isempty(ord) ? zeros(length(opt.xall), 0) :
            hcat([Fs[c][k][:] for (j, c, k) in ord]...)
    slot = gate_scalars(fslot, zbnd)
    slot["misfit"] = [parse(Float64, crc[c][j][9]) for (j, c, k) in ord]
    slot["k"]      = [parse(Float64, crc[c][j][8]) for (j, c, k) in ord]
    if all(c -> length(lnc[c]) == length(rows[c].T), eachindex(rows))
        slot["log_s"] = [lnc[c][j] for (j, c, k) in ord]
    end
    st, slot, tempering_stats([r.T for r in rows])
end

"Chain-id R-hat of one unit: the worst scalar, and its ESS (diagnostic)."
function chainid_check(streams, zbnd)
    pw = CFG.FBOUNDS[2] - CFG.FBOUNDS[1]
    S = [x for x in streams if size(x.fstar, 2) >= 8]
    length(S) >= 2 || return (rhat = NaN, ess = NaN)
    gs = [gate_scalars(x.fstar, zbnd) for x in S]
    vs = [rhat_and_ess([g[k] for g in gs]; priorwidth = startswith(k, "theta_") ? pw : nothing)
          for k in keys(gs[1])]
    L = [x.logs for x in S if length(x.logs) >= 8]
    length(L) >= 2 && push!(vs, rhat_and_ess(L))
    r = filter(isfinite, [v.rhat for v in vs]); q = filter(isfinite, [v.ess_bulk for v in vs])
    (rhat = isempty(r) ? NaN : maximum(r), ess = isempty(q) ? NaN : minimum(q))
end

"""
One unit: ensemble `e` (kernel draw e), replicate `r` (seed). Samples, then
keeps what the posterior and the gates need -- not the chains, which stay in
the run directory.
"""
function run_unit(b, sounding_id, e, r, pids)
    prefix = unit_prefix(sounding_id, e, r)
    clear_tag(pwd(), prefix)
    sd = unit_seed(sounding_id, e, r)
    @sync for p in pids
        @async remotecall_fetch(() -> (Random.seed!(sd + myid()); nothing), p)
    end
    Random.seed!(sd)
    S = hier(b)
    opt, optn = make_opts(prefix, b.zbnd)
    t = @elapsed transD_GP.main(opt, optn, S, pids;
                                nsamples = NSAMPLES, nchainsatone = NCHAINSATONE,
                                Tmax = TMAX)
    set_run_filenames!(opt, optn)
    inv = assert_chains_at_T1(opt.costs_filename, length(pids), NCHAINSATONE)
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
    st, slot, tp = t1_streams(opt, optn, b.zbnd)
    ci = chainid_check(st, b.zbnd)
    @printf("      ensemble %d replicate %d: %d draws, %.1f s, T=1 %s; T=1 draws by chain %s; swaps %.2f, round trips %s\n",
            e, r, size(fstar, 2), t, inv.ok ? "ok" : inv.why,
            join([size(x.fstar, 2) for x in st], "/"), tp.swap_rate, join(tp.round_trips, "/"))
    (e = e, r = r, prefix = prefix, seed = sd, secs = t, fstar = fstar, logs = logs,
     nodes = Float64.(h.nodes), U = h.U, t1_ok = inv.ok, slot = slot,
     chainid_rhat = ci.rhat, chainid_ess = ci.ess,
     swap_rate = tp.swap_rate, round_trips = collect(Int, tp.round_trips),
     frac_at_T1 = collect(Float64, tp.frac_at_T1))
end

unit_file(dir, e, r) = joinpath(dir, "unit_e$(e)r$(r).jld2")
save_unit(dir, u) = jldsave(unit_file(dir, u.e, u.r); pairs(u)...)
function load_unit(dir, e, r)
    d = load(unit_file(dir, e, r))
    (; (Symbol(k) => v for (k, v) in d)...)
end

"""
One site. With --ens and/or --rep only those units are sampled and saved,
and the site is not combined: that is how a PBS job runs several units of one
site at once, as separate processes. --combine samples nothing and gates
whatever units are on disk. With neither, every unit is sampled in turn and
then combined, as before.
"""
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
    want = [(e, rr) for e in 1:NENS for rr in 1:REPS]
    partial = haskey(ARG, "ens") || haskey(ARG, "rep")
    haskey(ARG, "ens") && (want = filter(u -> u[1] in parse.(Int, split(ARG["ens"], ",")), want))
    haskey(ARG, "rep") && (want = filter(u -> u[2] in parse.(Int, split(ARG["rep"], ",")), want))
    if !COMBINE_ONLY
        isempty(want) && error("--ens/--rep select no unit (ensembles 1:$NENS, replicates 1:$REPS)")
        old = pwd(); cd(dir)
        try
            for (e, rr) in want
                # Stage 2: this ensemble's kernel is built on its own Stage-1 draw.
                be = udf_kernel(b0, job, field_id, OUT; ens = e, occam = OCCAMK)
                save_unit(dir, run_unit(be, job.sounding_id, e, rr, pids))
            end
        finally
            cd(old)
        end
        if partial
            println("    wrote $(length(want)) unit(s) to $dir; combine with --combine")
            return
        end
    end

    allu = [(e, rr) for e in 1:NENS for rr in 1:REPS]
    missing_u = filter(u -> !isfile(unit_file(dir, u...)), allu)
    isempty(missing_u) || error("$field_id: missing units " *
        join(("e$(u[1])r$(u[2])" for u in missing_u), ", ") * " in $dir")
    units = [load_unit(dir, u...) for u in allu]

    G = apply_gates(b.zbnd, units)
    S = summarise(b, G.pooled)
    @printf("    pooled %d draws over %d ensembles x %d replicates\n", S.ndraws, NENS, REPS)
    @printf("    mixing: replicate rhat %.3f (ESS %.0f); chain-id rhat %.3f (diagnostic); kernel spread rhat %.3f\n",
            G.rhat_replicate_max, G.ess_replicate, G.rhat_chainid_max, G.rhat_kernel_max)
    if !isempty(G.rhat_replicate_by_scalar)
        println("    replicate rhat by scalar (worst ensemble): " *
                join((@sprintf("%s %.3f", k, G.rhat_replicate_by_scalar[k])
                      for k in sort(collect(keys(G.rhat_replicate_by_scalar)))), ", "))
    end
    @printf("    ladder: %d chains (%d at T=1), Tmax %.3g; min round trips per unit %s; stationary %s; gates %s (%s)\n",
            NCHAINS, NCHAINSATONE, TMAX, join(G.min_round_trips, "/"), G.stationary,
            G.gates_pass ? "PASS" : "FAIL", G.gate_basis)
    s_draws = vcat([u.logs for u in units if !isempty(u.logs)]...)
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
        log_s = s_draws, nodes = vcat([u.nodes for u in units]...),
        U = vcat([u.U for u in units]...),
        q05 = S.q05, q50 = S.q50, q95 = S.q95, ndraws = S.ndraws,
        doi_width = S.doi, doi_05 = S.doi_05, doi_07 = S.doi_07,
        rhat_max = G.rhat_max, ess_bulk_min = G.ess_bulk_min,
        ess_tail_min = G.ess_tail_min, stationary = G.stationary,
        splithalf_max = G.splithalf_max,
        gates_pass = G.gates_pass, gate_basis = G.gate_basis,
        rhat_mixing_max = G.rhat_mixing_max, ess_mixing = G.ess_mixing,
        rhat_replicate_max = G.rhat_replicate_max, ess_replicate = G.ess_replicate,
        rhat_replicate_by_scalar = G.rhat_replicate_by_scalar,
        rhat_chainid_max = G.rhat_chainid_max, ess_chainid = G.ess_chainid,
        rhat_kernel_max = G.rhat_kernel_max, min_round_trips = G.min_round_trips,
        swap_rate = G.swap_rate, frac_at_T1 = G.frac_at_T1,
        nchainsatone = NCHAINSATONE, tmax = TMAX, nreplicates = REPS,
        cond_source = b.cond_source, n_kernel_draws = b.n_kernel_draws,
        nsamples = NSAMPLES, nensembles = NENS, nchains = NCHAINS, s_tau = S_TAU,
        flat_sigmahat = FLAT, supplied_conductivity = SUPPLIED,
        seeds = [u.seed for u in units], secs = sum(u.secs for u in units),
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
    (haskey(ARG, "ens") || haskey(ARG, "rep")) && length(sites) != 1 &&
        error("--ens/--rep run units of ONE site; give --site")
    @printf("trans-D SNMR: %d sites, %d ensembles x %d replicates x %d chains (%d at T=1, Tmax %.3g), %d samples, s tau %.1f%s\n",
            length(sites), NENS, REPS, NCHAINS, NCHAINSATONE, TMAX, NSAMPLES, S_TAU,
            COMBINE_ONLY ? "  [combine only]" : "")
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
