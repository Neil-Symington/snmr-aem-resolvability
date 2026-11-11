# smr_setup.jl -- the SNMR trans-D settings, declared once.
#
# Included by run_smr_transd.jl (which samples) and report_smr.jl (which
# re-reads the chains to plot them). Both must build the SAME OptionsStat /
# OptionsNuisance and the SAME hierarchical sounding, or the report would
# describe a different run from the one on disk; so these were moved here from
# run_smr_transd.jl verbatim rather than copied.
#
# Needs, already loaded: stage_prelude.jl (OUT, FCFG, NSAMPLES_SNMR,
# NCHAINS_CFG), prerun_data.jl (CFG), HiQGA.transD_GP, SMRPInversion, Statistics.

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
const NCHAINS  = NCHAINS_CFG    # from [mcmc], not prerun_config.jl
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
    optn = transD_GP.OptionsNuisance(opt;
        sdev = copy(S_SDEV), bounds = copy(SBOUNDS),
        nsamples = NSAMPLES, updatenuisances = true)
    optn.fdataname = prefix
    opt, optn
end

"Hierarchical sounding with the widened s prior."
function hier(b)
    σ̂ = FLAT ? fill(median(b.σ̂), length(b.σ̂)) : b.σ̂
    SMRPInversion.newSMRSoundingHierarchical(b.V0, b.ϕ, b.F;
        σ̂_V0 = σ̂, s_prior_center = S_PRIOR_CENTER,
        s_prior_logwidth = S_TAU, s_init = CFG.S_INIT,
        linearsat = true, amponly = true)
end
