# aem_lib.jl — loading and inversion helpers for the AEM soundings extracted at
# the SMR sites. Written for UDF (SkyTEM); [aem] system now also selects the
# Birrindudu XCite. Mirrors ek_realdata/aem_lib.jl; the differences are the
# system module, the manifest columns that this delivery carries, and the
# column name `d_pV` (Ord-Keep's exporter wrote `data_pV`).

using DelimitedFiles, Statistics, Printf, LinearAlgebra
using HiQGA.transD_GP

# THE AEM SYSTEM, chosen by [aem] system. Each system module parses its own .des
# and cross-checks it, which is stage work, so both live in stages/. They export
# the same API -- build_grid, build_operator(g, height; ...), predicted(aem),
# data_blocks(g), gate_times(), describe(), MU0, PV, rTx -- and nothing below or
# in the stages names a system. A config without the key is a SkyTEM config:
# that is what every config written before this switch meant.
const AEM_SYSTEM = String(get(FieldConfig.CFG[]["aem"], "system", "skytem"))
if AEM_SYSTEM == "skytem"
    include(joinpath(@__DIR__, "..", "stages", "udf_system.jl"));  using .UDFSystem
elseif AEM_SYSTEM == "xcite"
    include(joinpath(@__DIR__, "..", "stages", "xcite_system.jl")); using .XCiteSystem
else
    error("[aem] system = \"$AEM_SYSTEM\": known systems are \"skytem\" and \"xcite\"")
end
include(joinpath(@__DIR__, "occam_udf.jl")); using .OccamUDF

const HERE = @__DIR__
# From [paths] aem_sounds, not beside this file. The stages were written
# in the directory that held the soundings; in the repository they are
# shipped under fieldrun/example/aem/soundings and the control file says so.
const SND  = FieldConfig.cfgpath(FieldConfig.CFG[], "aem_sounds")
# OUT comes from the prelude, which resolves [paths] outputs. Each stage and
# aem_lib.jl used to define its own `joinpath(HERE, "outputs")`, and because
# the support modules are included AFTER the prelude theirs won -- so the
# first real run wrote to fieldrun/common/outputs/ instead of where the
# control file said. Four separate definitions of one path, which is what
# the control file exists to prevent.
"Read a CSV with a header row into (header => column) pairs of Strings."
function readcsv_dict(path)
    raw, hdr = readdlm(path, ',', String; header = true)
    Dict(strip(h) => raw[:, i] for (i, h) in enumerate(vec(hdr)))
end

"The site manifest as a vector of named tuples."
function manifest()
    d = readcsv_dict(joinpath(SND, "manifest.csv"))
    n = length(d["field_id"])
    [(field_id = d["field_id"][i],
      line = parse(Int, d["line"][i]),
      fiducial = parse(Float64, d["fiducial"][i]),
      rec = parse(Int, d["rec"][i]),
      sep_m = parse(Float64, d["sep_m"][i]),
      easting = parse(Float64, d["easting"][i]),
      northing = parse(Float64, d["northing"][i]),
      height_m = parse(Float64, d["height_m"][i]),
      dtm_m = parse(Float64, d["dtm_m"][i]),
      plm = parse(Float64, d["plm"][i]),
      bore = d["bore"][i],
      bore_sep_m = parse(Float64, d["bore_sep_m"][i]),
      bore_td_m = parse(Float64, d["bore_td_m"][i])) for i in 1:n]
end

"""
Z-component data and sigma for one site, in delivered pV units.

Every system gets `d`, `σ`, `gate`, `moment` in the data order its operator
uses (SkyTEM: LM then HM). SkyTEM also keeps the per-moment split
(`d_lm`, `σ_lm`, ... ) that its operator is built from.
"""
function load_gates(field_id)
    d = readcsv_dict(joinpath(SND, "$(field_id).csv"))
    mom = String.(strip.(d["moment"]))
    dat = parse.(Float64, d["d_pV"])
    sig = parse.(Float64, d["sigma_pV"])
    gat = parse.(Int, d["gate"])
    if AEM_SYSTEM == "skytem"
        il, ih = mom .== "LM", mom .== "HM"
        o = [findall(il); findall(ih)]
        (d_lm = dat[il], σ_lm = sig[il], g_lm = gat[il],
         d_hm = dat[ih], σ_hm = sig[ih], g_hm = gat[ih],
         d = dat[o], σ = sig[o], gate = gat[o], moment = mom[o])
    else
        all(==("Z"), mom) || error("$field_id.csv: an XCite sounding has only moment \"Z\", found $(unique(mom))")
        issorted(gat) || error("$field_id.csv: gates are not in order")
        # halt_pV is the MEASURED high-altitude noise; sigma_pV is it combined
        # with the nominal_unc convention (prep_aem_xcite.jl). The two-term
        # calibration rescales only the convention, so it needs the floor.
        halt = haskey(d, "halt_pV") ? parse.(Float64, d["halt_pV"]) : Float64[]
        (d = dat, σ = sig, gate = gat, moment = mom, halt = halt)
    end
end

"""
    fit_unconstrained(g, height; grid, σ_scale, starts)

Smooth (Occam) fit with no constraint, from several uniform starts. Returns
`(aem, m, chi2, ndata)` with `m` in log10 conductivity on `grid.zall`.
"""
function fit_unconstrained(g, height; grid = build_grid(), σ_scale = 1.0,
                           starts = (-2.6, -2.0, -1.4), bounds = Bounds(),
                           m0val = -2.0)
    aem = build_operator(g, height;
                         grid = grid, calcjacobian = true, σ_scale = σ_scale)
    n = length(aem.ρ) - grid.nfixed
    ndata = length(aem.res)
    R = makeR1(n); m0 = fill(m0val, n)
    tr = Dict{Symbol,Any}()
    m, χ² = occam_multi(aem, [fill(s, n) for s in starts], m0, R, ndata;
                        mode = :minimise, bounds = bounds, trace = tr)
    (aem = aem, m = m, chi2 = χ², ndata = ndata, R = R, m0 = m0, n = n,
     l2 = get(tr, :l2, NaN), bounds = bounds)
end

"""
    effective_parameters(fit)

p_eff = tr(Jw (Jw'Jw + l2 (R'R + beta2 I))^-1 Jw'), Jw = W J, at the fitted model
and the lambda^2 it was solved at: the trace of the data-resolution matrix, so
the expected chi2 of a fit that absorbs only noise is ndata - p_eff, not ndata.
"""
function effective_parameters(fit)
    isfinite(fit.l2) || error("no Occam step was accepted, so lambda^2 is undefined and p_eff with it")
    chi2(fit.aem, fit.m; computeJ = true)
    Jw = Matrix(fit.aem.W * fit.aem.J)
    H = Jw' * Jw .+ fit.l2 .* (Matrix(fit.R' * fit.R) + OccamUDF.Β2 * I)
    p = tr(Jw * (H \ Jw'))
    chi2(fit.aem, fit.m)          # leave aem.res at the accepted model
    p
end

"Cells of a log10-conductivity model sitting on a bound."
on_bound(m, b::Bounds) = count(x -> x <= b.lo + 1e-9 || x >= b.hi - 1e-9, m)

"""
    calibrate_noise_twoterm(g, height; grid, tol, maxit, fmin, bounds)

TWO-SIDED, and only where the delivered noise is a convention. XCite sigma is
sqrt((u d)^2 + halt^2): u = nominal_unc is a processing convention, halt is the
MEASURED high-altitude noise, and on the late gates halt is most of sigma. So
this scales the convention alone,

    sigma_i(f) = sqrt((f u d_i)^2 + halt_i^2),

and never takes a gate below its measured noise, as a whole-sigma scale would.
f is solved for

    chi2_min(f) = ndata - p_eff(f),

the maximum-likelihood level with the degrees of freedom the fit absorbs:
calibrating to chi2 = ndata instead ignores p_eff and overtightens. p_eff
depends on the lambda^2 the fit lands on, which depends on f, so this is a
fixed point; p_eff is recomputed every pass and every pass is reported.

REFUSES, rather than returning a number, when (1) the passes do not settle in
`maxit`, (2) f would fall below `fmin`, or (3) the fixed point puts any cell on a
model bound. (3) is the important one: a rougher model resolves more parameters
and fits better, so targeting a misfit can drift toward an overfit, and the
synthetics' lambda^2 = 0.01 models pin 3-28 of 65 cells on the bounds where the
smoothest fitting models pin none. A noise level calibrated to a bound-pinned
model is calibrated to an overfit.
"""
function calibrate_noise_twoterm(g, height; grid = build_grid(), tol = 0.02, maxit = 8,
                                 fmin = 0.01, bounds = Bounds(), verbose = true)
    AEM_SYSTEM == "xcite" || error("the two-term noise calibration is written for XCite only")
    length(g.halt) == length(g.σ) ||
        error("no halt_pV column in the sounding file; the two-term calibration needs the measured floor")
    conv = sqrt.(max.(g.σ .^ 2 .- g.halt .^ 2, 0.0))       # the u*d term, as delivered
    all(>(0), conv) || error("sigma_pV does not exceed halt_pV at every gate; cannot separate the terms")
    σof(f) = sqrt.((f .* conv) .^ 2 .+ g.halt .^ 2)
    hist = NamedTuple[]
    f = 1.0
    verbose && @printf("    %4s %8s %9s %7s %8s %7s %9s %10s %6s\n",
                       "pass", "f", "chi2", "p_eff", "N-p_eff", "ratio", "lambda2", "roughness", "bound")
    for it in 1:maxit
        fit = fit_unconstrained(merge(g, (σ = σof(f),)), height; grid = grid, bounds = bounds)
        p = effective_parameters(fit)
        target = fit.ndata - p
        r = fit.chi2 / target
        h = (pass = it, f = f, chi2 = fit.chi2, p_eff = p, target = target, ratio = r,
             l2 = fit.l2, roughness = norm(fit.R * fit.m)^2, n_bound = on_bound(fit.m, bounds))
        push!(hist, h)
        verbose && @printf("    %4d %8.4f %9.4f %7.2f %8.2f %7.4f %9.3g %10.4g %6d\n",
                           h.pass, h.f, h.chi2, h.p_eff, h.target, h.ratio, h.l2, h.roughness, h.n_bound)
        if abs(r - 1) <= tol
            h.n_bound == 0 || error("the calibrated fixed point (f = $(round(f, sigdigits=4))) puts ",
                "$(h.n_bound) cells on a model bound: the noise is calibrated to an overfit and is not usable")
            return (f = f, σ = σof(f), fit = fit, p_eff = p, history = hist)
        end
        # Secant in log f on log(ratio) once two passes exist: chi2 responds to f
        # more weakly than 1/f^2 because halt does not scale. First pass: the
        # pure-convention guess. Steps are held to a factor of 3.
        lf = if it == 1
            log(f) + 0.5 * log(r)
        else
            a, b = hist[end-1], hist[end]
            sl = (log(b.ratio) - log(a.ratio)) / (log(b.f) - log(a.f))
            isfinite(sl) && sl < 0 ? log(b.f) - log(b.ratio) / sl : log(b.f) + 0.5 * log(b.ratio)
        end
        f = exp(clamp(lf, log(f) - log(3), log(f) + log(3)))
        f < fmin && error("f would fall below $fmin (to $(round(f, sigdigits=3))): the measured floor ",
                          "alone overstates the misfit, or the fit absorbs it; not calibrating")
    end
    error("the noise calibration did not settle in $maxit passes; not calibrating (passes above)")
end

"""
    calibrated_gates(g, A)

The gates and sigma scale to build the operator with, from what run_aem_occam.jl
stored in `A`: the two-term calibrated per-gate sigma, or the delivered sigma and
the scalar (inflate-only) scale. Every stage that rebuilds the operator goes
through here, so none can apply the scalar to a two-term calibration.
"""
function calibrated_gates(g, A)
    cal = get(A, "noise_calibration", "inflate")
    if cal == "two_term"
        σc = A["sigma_cal"]
        length(σc) == length(g.σ) || error("sigma_cal has $(length(σc)) gates, the sounding $(length(g.σ))")
        return merge(g, (σ = σc,)), 1.0
    end
    cal == "inflate" || error("unknown noise_calibration \"$cal\" in the Occam output")
    g, A["noise_scale"]
end

"""
    calibrate_noise(g, height; grid, tol, maxit)

The admissibility level the report uses is chi2 <= ndata + 3.84, which presumes
a noise model under which a good model reaches chi2/ndata ~ 1. The contractor's
3%-plus-relative-uncertainty model is a processing convention, not a measured
error bar on the forward operator we are using here, so on real data it lands
wherever it lands. This rescales sigma by a single scalar until the best smooth
fit sits at chi2/ndata = 1, and reports the factor.

This is the two-pass calibration already used for the SNMR soundings in
ek_inversions (sigma_fit -> sigma_cal there), applied to AEM: the maximum
likelihood noise scale for this model class. It makes chi2 comparable between
sites and makes "ndata + 3.84" mean what it means in the report.

INFLATE ONLY. The scale is floored at 1, so the contractor's error bars are
never made SMALLER. Measured here, the smooth fit reaches chi2/ndata below 1 at
every site under the delivered noise model, so a two-sided calibration would
shrink sigma and hand back a TIGHTER admissible band -- a stronger claim
resting on a weaker justification. The admissibility level is referenced to
ndata, not to the achieved minimum, so a minimum below ndata needs no
correction; only a minimum far above it does, and that is model error, which
inflation is the honest way to absorb.
"""
function calibrate_noise(g, height; grid = build_grid(), tol = 0.02, maxit = 6,
                         bounds = Bounds(), verbose = true)
    scale = 1.0
    fit = nothing
    for it in 1:maxit
        fit = fit_unconstrained(g, height; grid = grid, σ_scale = scale,
                                bounds = bounds)
        r = fit.chi2 / fit.ndata
        verbose && @printf("    pass %d: sigma x %.3f -> chi2/N %.4f\n", it, scale, r)
        (r <= 1 + tol) && break
        scale *= sqrt(r)
    end
    (scale = scale, fit = fit)
end

"Geometric mean conductivity over a depth window of a log10-conductivity model."
function unit_geomean(m, zall, za, zb)
    k = findall(z -> za <= z <= zb, zall)
    isempty(k) && error("window $za-$zb m contains no cells")
    10^mean(m[k])
end

"Arithmetic mean (conductance-weighted) conductivity over the same window."
function unit_arithmean(m, zall, za, zb)
    k = findall(z -> za <= z <= zb, zall)
    mean(10 .^ m[k])
end
