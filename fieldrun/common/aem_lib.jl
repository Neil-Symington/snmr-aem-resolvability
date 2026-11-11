# aem_lib.jl — loading and inversion helpers for the UDF AEM soundings
# extracted at the four chosen SMR sites. Mirrors ek_realdata/aem_lib.jl; the
# differences are the system module, the manifest columns that this delivery
# carries, and the column name `d_pV` (Ord-Keep's exporter wrote `data_pV`).

using DelimitedFiles, Statistics, Printf, LinearAlgebra
using HiQGA.transD_GP

# udf_system.jl is a STAGE, not common: it parses the SkyTEM .des and
# cross-checks it, which is stage work. This include pointed beside itself and
# only loading the file found that -- it parses either way.
include(joinpath(@__DIR__, "..", "stages", "udf_system.jl"));  using .UDFSystem
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

"Z-component data and sigma for one site, split by moment, in delivered pV units."
function load_gates(field_id)
    d = readcsv_dict(joinpath(SND, "$(field_id).csv"))
    mom = strip.(d["moment"])
    dat = parse.(Float64, d["d_pV"])
    sig = parse.(Float64, d["sigma_pV"])
    gat = parse.(Int, d["gate"])
    il, ih = mom .== "LM", mom .== "HM"
    (d_lm = dat[il], σ_lm = sig[il], g_lm = gat[il],
     d_hm = dat[ih], σ_hm = sig[ih], g_hm = gat[ih])
end

"""
    fit_unconstrained(g, height; grid, σ_scale, starts)

Smooth (Occam) fit with no constraint, from several uniform starts. Returns
`(aem, m, chi2, ndata)` with `m` in log10 conductivity on `grid.zall`.
"""
function fit_unconstrained(g, height; grid = build_grid(), σ_scale = 1.0,
                           starts = (-2.6, -2.0, -1.4), bounds = Bounds(),
                           m0val = -2.0)
    aem = build_operator(g.d_lm, g.σ_lm, g.d_hm, g.σ_hm, height;
                         grid = grid, calcjacobian = true, σ_scale = σ_scale)
    n = length(aem.ρ) - grid.nfixed
    ndata = length(aem.res)
    R = makeR1(n); m0 = fill(m0val, n)
    m, χ² = occam_multi(aem, [fill(s, n) for s in starts], m0, R, ndata;
                        mode = :minimise, bounds = bounds)
    (aem = aem, m = m, chi2 = χ², ndata = ndata, R = R, m0 = m0, n = n)
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
