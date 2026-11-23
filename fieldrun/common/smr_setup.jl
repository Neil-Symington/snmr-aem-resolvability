# smr_setup.jl -- the SNMR trans-D settings, declared once.
#
# Included by run_smr_transd.jl (which samples) and report_smr.jl (which
# re-reads the chains to plot them). Both must build the SAME OptionsStat /
# OptionsNuisance and the SAME hierarchical sounding, or the report would
# describe a different run from the one on disk; so these were moved here from
# run_smr_transd.jl verbatim rather than copied.
#
# Needs, already loaded: stage_prelude.jl (OUT, FCFG, NSAMPLES_SNMR,
# NCHAINS_CFG), prerun_data.jl (CFG), HiQGA.transD_GP, SMRPInversion, Statistics,
# and smr_marginal.jl (module SMRMarginal) on every process.

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        if startswith(ARGS[i], "--")
            if i < length(ARGS) && !startswith(ARGS[i+1], "--")
                ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2
            else
                ARG[ARGS[i][3:end]] = "true"; i += 1
            end
        else; i += 1; end
    end
end

const RUNDIR   = joinpath(OUT, "smr_mcmc")
const NSAMPLES = NSAMPLES_SNMR  # [mcmc] nsamples_snmr, or --nsamples
# NENS is the ensemble count AND the kernel-draw count -- one kernel per
# ensemble, the pooled posterior marginalising over them (see run_ensemble).
# It was hardcoded 4, so [kernel] n_draws and [mcmc] nensembles were both dead
# keys and the smoke overlay's kernel_draws = 1 never took effect: "smoke" ran
# four ensembles, not one. [kernel] n_draws governs it now, being the key the
# overlay already maps onto and the one documented with the reason. In full
# mode it is 4, as every run so far has used, so no full-mode result moves.
const NENS     = haskey(ARG, "ensembles") ? parse(Int, ARG["ensembles"]) :
                 Int(FCFG["kernel"]["n_draws"])
# THE SNMR LADDER: [mcmc] nchains_snmr / nchainsatone_snmr / tmax_snmr, each
# falling back to the shared nchains / 1 / tmax (stage_prelude.jl). TMAX was
# CFG.TMAX, the 2.5 hardcoded in src/prerun_config.jl -- [mcmc] tmax was
# printed by the prelude and never reached this stage. A test holds the two
# equal, so the shared value was only ever right by agreement; the ladder is
# now read from the control file.
const NCHAINS      = NCHAINS_SNMR
# R independent ladders per kernel draw ([mcmc] replicates_snmr, or
# --replicates). Ensemble e, replicate r is one "unit": its own seed, its own
# chain files, and with --ens/--rep its own process, so a node can run several.
const REPS         = haskey(ARG, "replicates") ? parse(Int, ARG["replicates"]) :
                     REPLICATES_SNMR
@assert REPS >= 1 "[mcmc] replicates_snmr must be >= 1, got $REPS"
const NCHAINSATONE = NCHAINSATONE_SNMR
const TMAX         = TMAX_SNMR
@assert 1 <= NCHAINSATONE < NCHAINS "SNMR ladder needs 1 <= nchainsatone ($NCHAINSATONE) < nchains ($NCHAINS)"
# THE [noise] SECTION WAS DECORATIVE. fieldconfig.jl validated all five keys
# and the stage then used none of them: s_prior_tau, s_sdev, s_bounds_log10 and
# sigma_hat_branch came from hardcoded constants here or from prerun_config.jl,
# so editing the control file changed nothing. Two of them also DISAGREED with
# the code, which is worse than unused:
#
#   s_prior_tau        control file 0.5    this stage 2.0
#   s_bounds_log10     control file 1e2    this stage 1e4 (upper)
#
# 2.0 is the right tau and the header above says why: at 0.5 the East Kimberley
# posterior s sat at chi2/N 2.90 against SSC's 1.02, the prior being too narrow
# by about 1.7x in sigma. The 1e4 upper bound is what the working version ran
# with, and it is not academic -- the smoke posterior puts s in the tens, so
# an upper bound of 1e2 would pin it where 1e4 does not. Both control-file
# values have been corrected to what the code actually used, so wiring them
# through changes no result; from here the control file governs.
const S_TAU    = haskey(ARG, "stau") ? parse(Float64, ARG["stau"]) :
                 Float64(FCFG["noise"]["s_prior_tau"])
const SBOUNDS  = let b = Float64.(FCFG["noise"]["s_bounds_log10"])
    # the control file states the bounds as log10, transD_GP wants natural log
    reshape([log(10.0^b[1]), log(10.0^b[2])], 1, 2)
end
const S_SDEV   = Float64.(FCFG["noise"]["s_sdev"])
const S_PRIOR_CENTER = Float64(FCFG["noise"]["s_prior_center"])
# HOW s IS HANDLED ([noise] s_treatment). "sampled" (the default): HiQGA samples
# u = log s as a nuisance, every run before 2026-11-23. "marginal": s is
# integrated out of the likelihood (common/smr_marginal.jl), the same posterior
# on w with one coordinate fewer to sample; s is drawn after the run. The prior
# (centre, tau) and bounds are the same in both.
const S_TREATMENT = String(get(FCFG["noise"], "s_treatment", "sampled"))
const MARGINAL    = S_TREATMENT == "marginal"
# --supplied: use the operator's resistivity profile instead of the AEM model,
# as a kernel-sensitivity control.
const SUPPLIED = false   # no operator profile exists for UDF
# --flat: sigma-hat_j = median(sigma-hat_j) for every j. The LEVEL is preserved
# and the per-moment SHAPE removed. The pre-run measured that shape rising 2.5x
# with pulse moment and tracking signal amplitude, which is decay-model
# mismatch rather than the additive noise the physics implies; it is also the
# spec's own fallback. It matters here because the inferred noise level sets
# the constraint depth, which sets what the detection criteria can see.
const FLAT     = haskey(ARG, "flat")
# --occam: force the smooth Occam conductivity for the kernel even when a
# trans-D ensemble exists. The control for whether a Stage-1 conductivity
# change is what moved the Stage-2 water.
const OCCAMK   = haskey(ARG, "occam")
const LABEL    = (SUPPLIED ? "_supplied" : "") * (FLAT ? "_flat" : "") *
                 (OCCAMK ? "_occamk" : "")

function make_opts(prefix, zbnd)
    nz = length(zbnd)
    opt = transD_GP.OptionsStat(
        nmin = CFG.NMIN, nmax = CFG.NMAX,
        xbounds = [1.0 Float64(nz)], fbounds = CFG.FBOUNDS, fdataname = prefix,
        xall = permutedims(collect(1:1.0:nz)), λ = CFG.LAMBDA_GP, δ = CFG.DELTA,
        sdev_prop = 0.07 * diff(CFG.FBOUNDS, dims = 2)[:],
        sdev_pos  = [0.05 * nz],
        sdev_dc   = 0.008 * diff(CFG.FBOUNDS, dims = 2)[:],
        save_freq = CFG.SAVE_FREQ, quasimultid = false, K = CFG.GP_KERNEL)
    # with s integrated out there is no nuisance to sample: optn is nothing
    MARGINAL && return opt, nothing
    optn = transD_GP.OptionsNuisance(opt;
        sdev = copy(S_SDEV), bounds = copy(SBOUNDS),
        nsamples = NSAMPLES, updatenuisances = true)
    optn.fdataname = prefix
    opt, optn
end

"""
The SNMR sounding: hierarchical (s sampled) or marginal (s integrated out),
with the same noise shape, prior and bounds either way.
"""
function hier(b)
    σ̂ = FLAT ? fill(median(b.σ̂), length(b.σ̂)) : b.σ̂
    MARGINAL && return SMRMarginal.newSMRSoundingMarginal(b.V0, b.ϕ, b.F;
        σ̂_V0 = σ̂, s_prior_center = S_PRIOR_CENTER, s_prior_logwidth = S_TAU,
        s_bounds_ln = (SBOUNDS[1], SBOUNDS[2]), linearsat = true)
    SMRPInversion.newSMRSoundingHierarchical(b.V0, b.ϕ, b.F;
        σ̂_V0 = σ̂, s_prior_center = S_PRIOR_CENTER,
        s_prior_logwidth = S_TAU, s_init = CFG.S_INIT,
        linearsat = true, amponly = true)
end

"""
The chain-file prefix of ensemble `e`, replicate `r`.

With one replicate this is the prefix every earlier run used, so old run
directories still read. With more it is `_e<e>r<r>_`: never a substring of
another unit's prefix, which matters because each unit clears its own prefix
before sampling, and several units share the directory and may run at once.
"""
unit_prefix(sounding_id, e, r = 1) =
    REPS == 1 ? "$(sounding_id)$(LABEL)_e$(e)_" : "$(sounding_id)$(LABEL)_e$(e)r$(r)_"

"Seed of one unit: replicate 1 keeps the single-ladder seed."
unit_seed(sounding_id, e, r = 1) =
    seed_for(sounding_id * LABEL * (r == 1 ? "" : "_r$r"), e)
