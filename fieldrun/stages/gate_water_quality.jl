# Ported from udf_realdata/ into the repository. The include paths and the
# project are the only things changed from the working version: COMMON is
# fieldrun/common/ and the project is the repository root, not ../ek_prerun.
# The control file first: support modules read it at load time.
include(joinpath(@__DIR__, "..", "common", "stage_prelude.jl"))
# gate_water_quality.jl — "what is the water quality?", two ways.
#
#   julia --project=. gate_water_quality.jl
#
# (1) THE GATE, prior-free. synthetics/gate_threshold.jl, with the truth's water
#     content and clay term replaced by what is actually available at a real
#     site. The gate certifies the water fresher than X when EVERY conductivity
#     the AEM data admit implies pore water below X:
#
#         sigma_crit(X) = geomean over the interval of [X*theta_total^m + surf]
#         X_cert        = min { X : sigma_crit(X) >= sigma_hi }
#
#     with sigma_hi the upper edge of the 95% admissible band from
#     run_aem_profile.jl. At a real site neither theta_total nor the clay term
#     is known, so the HEADLINE uses the conservative assumption-free limit:
#
#         clay = 0  =>  theta_total = theta_SMR  and  surf = 0
#
#     Both substitutions push X_cert UP, i.e. they weaken the certification:
#     setting surf = 0 attributes all of sigma_bulk to pore water, and
#     clay = 0 makes theta_total as small as it can be, so theta_total^m is
#     small and X_cert = sigma_hi / geomean(theta_SMR^m) is as large as it can
#     be. A site that certifies under this limit certifies under any clay model.
#     The clay-prior version is reported beside it and is always tighter.
#
#     theta_SMR comes from the SNMR posterior DRAWS and m from its U(1.5, 2.2)
#     prior, so X_cert is a distribution; the conservative quote is its p95.
#
# (2) THE WATER-QUALITY DISTRIBUTION, through the rock physics. The SNMR theta
#     draws are paired with sigma_bulk profiles drawn from the AEM PROFILE
#     LIKELIHOOD -- the constrained models of run_aem_profile.jl, resampled with
#     weight exp(-(chi2 - chi2_min)/2) -- and pushed through rockphys.jl.
#
#     THIS IS NOT A POSTERIOR AND IS NOT LABELLED AS ONE. It is a
#     likelihood-weighted ensemble over the one direction the gate constrains,
#     the interval's mean conductivity, with the smooth-model shape attached. A
#     trans-D AEM posterior would add the shape freedom this misses. What it
#     does do is propagate intact PROFILES rather than a percentile profile,
#     which is the standing requirement here: pushing a p50 sigma_bulk through a
#     non-linear transform moved published numbers by a factor of 4 when it was
#     caught, so it is not done.

using JLD2, Printf, Statistics, Random, DelimitedFiles

include(joinpath(COMMON, "rockphys.jl")); using .RockPhys
using .RockPhys: unit_sigma_w, unit_theta_total

const OUT = joinpath(@__DIR__, "outputs")
const NTHETA = 4000       # SNMR draws subsampled for the transform
const NPER   = 3          # clay/m draws per SNMR draw
const NSB    = 1200       # sigma_bulk profiles resampled from the profile likelihood
const FRESH  = 0.15       # S/m, the report's fresh line (~1,000 mg/L at 6400 mg/L per S/m)
const LEVELS = ((0.078, "500 mg/L  good drinking"),
                (0.156, "1,000 mg/L  drinking limit"),
                (0.234, "1,500 mg/L"),
                (0.469, "3,000 mg/L  dairy cattle"),
                (0.781, "5,000 mg/L"),
                (1.094, "7,000 mg/L  sheep"))

"Piecewise-constant resample of v(z_src) onto z_dst (aem_setup.jl convention)."
function piecewise(z_dst, z_src, v)
    out = similar(collect(Float64, z_dst))
    for (i, zt) in enumerate(z_dst)
        j = clamp(searchsortedlast(z_src, zt), 1, length(v))
        out[i] = v[j]
    end
    out
end

"""
sigma_bulk profiles from the AEM profile likelihood.

The stored constrained models span the admissible range of the interval's mean
conductivity. Resampling them with weight exp(-Dchi2/2) gives an ensemble whose
spread in that direction is the likelihood's, with no prior and no sampler.

INTERPOLATED, not resampled from the 16 grid points directly. The profile is
steep -- tens of chi2 across the admissible band -- so exp(-Dchi2/2) puts almost
all the weight on two or three grid points and a direct resample has an
effective sample size of about 3, which would make sigma_bulk a three-valued
discrete variable. The misfit curve is smooth in log sigma, so chi2 is
interpolated onto a fine grid, and each drawn model is the linear interpolation
(in log10 sigma) of its two bracketing constrained models, shifted so the
window's mean log10 conductivity equals the drawn target exactly. That last
shift is what keeps the ensemble consistent with the constraint the profile
was computed under.
"""
function sigma_bulk_ensemble(P, z_dst, nsb; seed = 20261026, nfine = 600)
    χ = P["chi2_prof"]; M = P["models"]; zall = P["zall"]
    fix = P["fixidx"]
    lt = log10.(P["targets"])
    ok = findall(isfinite, χ)
    length(ok) < 2 && error("profile has fewer than two finite points")
    ltok, χok = lt[ok], χ[ok]
    lf = collect(range(first(ltok), last(ltok), length = nfine))
    # linear interpolation of chi2 onto the fine grid
    χf = similar(lf)
    for (i, x) in enumerate(lf)
        j = clamp(searchsortedlast(ltok, x), 1, length(ltok) - 1)
        t = (x - ltok[j]) / (ltok[j+1] - ltok[j])
        χf[i] = (1 - t) * χok[j] + t * χok[j+1]
    end
    w = exp.(-(χf .- minimum(χf)) ./ 2)
    w ./= sum(w)
    rng = MersenneTwister(seed)
    cw = cumsum(w)
    SB = Matrix{Float64}(undef, length(z_dst), nsb)
    drawn = Vector{Float64}(undef, nsb)
    for k in 1:nsb
        i = clamp(searchsortedfirst(cw, rand(rng)), 1, nfine)
        x = lf[i]; drawn[k] = x
        j = clamp(searchsortedlast(ltok, x), 1, length(ltok) - 1)
        t = (x - ltok[j]) / (ltok[j+1] - ltok[j])
        m = (1 - t) .* M[:, ok[j]] .+ t .* M[:, ok[j+1]]
        m .+= x - mean(m[fix])          # put the window exactly on the target
        SB[:, k] = piecewise(z_dst, zall, 10 .^ m)
    end
    (SB = SB, weights = w, drawn = drawn, neff = 1 / sum(abs2, w),
     nfine = nfine, lf = lf, chi2f = χf)
end

"X_cert by bisection: the lowest pore-water conductivity the band can certify."
function xcert(θt_geo_pow, σ_hi)
    # sigma_crit(X) = X * geomean(theta_total^m) when the surface term is zero,
    # so this inverts analytically; kept as a function so the clay-prior version
    # can bisect the same way.
    σ_hi / θt_geo_pow
end

function run_site(field_id)
    pf = joinpath(OUT, "aem_profile_$(field_id).jld2")
    sf = joinpath(OUT, "smr_transd_$(field_id).jld2")
    df = joinpath(OUT, "detect_$(field_id).jld2")
    all(isfile, (pf, sf, df)) || (println("  $field_id: missing inputs, skipping"); return nothing)
    P = load(pf); S = load(sf); Dt = load(df)

    za, zb = P["win"]
    σ_hi = P["edge_hi"]; σ_lo = P["edge_lo"]
    Θall = S["theta"]; z = S["zbnd"]
    idx = findall(zi -> za <= zi <= zb, z)
    nd = size(Θall, 2)
    rng = MersenneTwister(4242)
    sel = nd > NTHETA ? sort(randperm(rng, nd)[1:NTHETA]) : collect(1:nd)
    Θ = Θall[:, sel]

    println("="^100)
    @printf("%s — interval %.1f-%.1f m (%s), %d SNMR draws used of %d\n",
            field_id, za, zb, Dt["win_source"], size(Θ, 2), nd)
    if !isfinite(σ_hi)
        println("  NO admissible conductivity band — the forced profile never came")
        println("  within ndata + 3.84 of the level. The gate cannot be evaluated here;")
        println("  this is a statement about the misfit floor, not about the water.")
        return nothing
    end
    @printf("  AEM 95%% admissible band: %.5f - %.5f S/m (%.0f - %.0f ohm.m), %.2f decades\n",
            σ_lo, σ_hi, 1 / σ_hi, 1 / σ_lo, log10(σ_hi / σ_lo))
    if P["edge_hi_open"]
        println("  WARNING: the band is OPEN at the conductive end — it runs off the")
        println("  forced grid, so sigma_hi is a grid limit, not a crossing, and X_cert")
        println("  below is a LOWER bound on the certifiable threshold.")
    end
    P["edge_lo_open"] && println("  NOTE: the band is open at the resistive end (expected; the")
    P["edge_lo_open"] && println("  gate is one-sided — certifying fresh needs only the upper edge).")

    # ── (1) the gate, assumption-free limit ────────────────────────────────
    dz = Dt["dz"]
    θ_int = vec(sum(Θ[idx, :] .* dz[idx], dims = 1) ./ sum(dz[idx]))   # interval mean theta
    ms = 1.5 .+ 0.7 .* rand(rng, length(θ_int))
    # geomean over the interval of theta^m, per draw
    θpow = [10^mean(ms[j] .* log10.(max.(Θ[idx, j], 1e-6))) for j in axes(Θ, 2)]
    Xc = [xcert(θpow[j], σ_hi) for j in axes(Θ, 2)]
    @printf("\n  GATE (clay = 0, the conservative limit)\n")
    @printf("    interval mean theta_SNMR: p05 %.4f  p50 %.4f  p95 %.4f v/v\n",
            quantile(θ_int, .05), quantile(θ_int, .5), quantile(θ_int, .95))
    @printf("    X_cert: p05 %.4f  p50 %.4f  p95 %.4f S/m\n",
            quantile(Xc, .05), quantile(Xc, .5), quantile(Xc, .95))
    @printf("    conservative quote (p95): %.4f S/m = %.0f mg/L\n",
            quantile(Xc, .95), quantile(Xc, .95) * MGL)
    verdict = quantile(Xc, .95) <= FRESH ? "CERTIFIED FRESH" : "NOT CERTIFIED at the 1,000 mg/L line"
    @printf("    verdict at %.2f S/m (~1,000 mg/L): %s\n", FRESH, verdict)
    println("\n    certification against management thresholds:")
    for (X, lbl) in LEVELS
        @printf("      %-28s %s\n", lbl, quantile(Xc, .95) <= X ? "certified" : "—")
    end

    # ── (2) the water-quality distribution ─────────────────────────────────
    E = sigma_bulk_ensemble(P, z, NSB)
    @printf("\n  sigma_bulk ensemble from the profile likelihood: %d profiles drawn\n", NSB)
    @printf("    over %d interpolated points; effective sample size %.0f\n", E.nfine, E.neff)
    sbu = vec(10 .^ mean(log10.(E.SB[idx, :]), dims = 1))
    @printf("    interval sigma_bulk: p05 %.4f  p50 %.4f  p95 %.4f S/m\n",
            quantile(sbu, .05), quantile(sbu, .5), quantile(sbu, .95))

    # pair each theta draw with an independently drawn sigma_bulk profile: the
    # two datasets are independent, so independent pairing is the correct joint
    SBpair = E.SB[:, rand(rng, 1:NSB, size(Θ, 2))]

    results = Dict{String,Any}()
    for (tag, cp) in (("clay prior 0.10-2.0 S/m", ClayPrior(:sigclay, 0.10, 2.0)),
                      ("product prior 0.002-0.8 S/m", ClayPrior(:product, 0.002, 0.8)),
                      ("no clay credit", ClayPrior(:none, 0.0, 0.0)),
                      ("clay zero (gate limit)", ClayPrior(:zero, 0.0, 0.0)))
        En = RockPhys.ensemble(Θ, SBpair, z, cp; nper = NPER, window = (za, zb),
                               seed = 20261026)
        if En.n == 0
            @printf("\n    %-30s NO admissible draws (accept %.3f)\n", tag, En.accept)
            continue
        end
        σw = unit_sigma_w(En, idx)
        θt = unit_theta_total(En, idx)
        pfresh = mean(σw .< FRESH)
        results[tag] = (σw = σw, pfresh = pfresh, accept = En.accept, n = En.n,
                        θt = θt, clay = vec(mean(En.clay[idx, :], dims = 1)))
        @printf("\n    %-30s accept %.3f, %d draws\n", tag, En.accept, En.n)
        @printf("      sigma_w  p05 %.4f  p50 %.4f  p95 %.4f S/m   (%.0f / %.0f / %.0f mg/L)\n",
                quantile(σw, .05), quantile(σw, .5), quantile(σw, .95),
                quantile(σw, .05) * MGL, quantile(σw, .5) * MGL, quantile(σw, .95) * MGL)
        @printf("      P(fresh, sigma_w < %.2f S/m) = %.3f\n", FRESH, pfresh)
        @printf("      implied clay fraction over the interval: p50 %.3f\n",
                quantile(results[tag].clay, .5))
    end

    main = get(results, "product prior 0.002-0.8 S/m", nothing)
    jldsave(joinpath(OUT, "gate_$(field_id).jld2");
        field_id = field_id, win = [za, zb], win_source = Dt["win_source"],
        sigma_lo = σ_lo, sigma_hi = σ_hi,
        theta_interval = θ_int, m_draws = ms, theta_pow = θpow,
        Xcert = Xc, Xcert_p95 = quantile(Xc, .95), Xcert_p50 = quantile(Xc, .5),
        mgl = MGL, fresh_line = FRESH, verdict = verdict,
        sigma_bulk_unit = sbu, sb_weights = E.weights, sb_neff = E.neff,
        sb_drawn = E.drawn, sb_lf = E.lf, sb_chi2f = E.chi2f,
        tags = collect(keys(results)),
        sigma_w = Dict(k => v.σw for (k, v) in results),
        pfresh = Dict(k => v.pfresh for (k, v) in results),
        accept = Dict(k => v.accept for (k, v) in results),
        clay_interval = Dict(k => v.clay for (k, v) in results),
        ntheta = size(Θ, 2), nper = NPER, nsb = NSB)
    (field_id = field_id, win = (za, zb), σ_lo = σ_lo, σ_hi = σ_hi,
     Xc95 = quantile(Xc, .95), Xc50 = quantile(Xc, .5), verdict = verdict,
     pfresh = main === nothing ? NaN : main.pfresh,
     σw50 = main === nothing ? NaN : quantile(main.σw, .5))
end

function main()
    d = readdlm(joinpath(FieldConfig.cfgpath(FieldConfig.CFG[], "aem_sounds"), "manifest.csv"), ',', String; header = true)
    ci = Dict(strip(h) => i for (i, h) in enumerate(vec(d[2])))
    rows = []
    for s in d[1][:, ci["field_id"]]
        r = try run_site(s) catch e
            @warn "gate failed for $s" exception = (e, catch_backtrace()); nothing
        end
        r !== nothing && push!(rows, r)
    end
    println("\n", "="^108)
    @printf("%-14s %-15s %10s %10s %11s %9s %9s  %s\n", "site", "interval m",
            "sigma_lo", "sigma_hi", "X_cert p95", "mg/L", "P(fresh)", "verdict")
    for r in rows
        @printf("%-14s %6.1f-%-8.1f %10.5f %10.5f %11.4f %9.0f %9.3f  %s\n",
                r.field_id, r.win[1], r.win[2], r.σ_lo, r.σ_hi, r.Xc95,
                r.Xc95 * MGL, r.pfresh, r.verdict)
    end
    println("="^108)
    println("X_cert p95 is the conservative, clay-free certification: every conductivity")
    println("the AEM data admit implies pore water fresher than this, with no clay credit.")
end

main()
