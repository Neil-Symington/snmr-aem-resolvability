# section6_bore_value.jl -- Section 6. Narrowing the Stage 0 priors with
# borehole logs. Test two only.
#
#   julia --project=. section6_bore_value.jl
#
# THIS SCRIPT READS TRUTH ON PURPOSE AND IS NOT PART OF THE STAGE 0-6 CHAIN.
# It simulates borehole logs, and a simulated log has to be centred on
# something. Nothing it computes feeds a stage. The firewall applies to the
# inversion, not to an information-value experiment about it.
#
# IT DOES NOT RE-RUN THE MONTE CARLO. Stage 5 already stored its 40,000 draws
# -- sigma_bulk, clay_product, theta_draw, m_draw -- so the four states are the
# SAME draws with the clay prior swapped and/or an induction likelihood applied
# as a weight. That is what makes the survey-only row reproduce 5.2 exactly
# rather than to within Monte Carlo error: it IS 5.2, recomputed from its own
# arrays. The reproduction is asserted below and the script stops if it fails.
#
# THE MODEL, as everywhere else:
#     sigma_bulk = sigma_w * theta^m + cp          cp = sigma_clay * C   [S/m]
#
# THE FOUR STATES
#   survey only  Stage 0 priors unchanged.
#   gamma        the clay prior replaced by a log-normal centred on the true
#                window clay term, 95% interval x/div 2.50.
#   induction    each Stage 1 draw WEIGHTED by a log-normal likelihood on its
#                window conductance, centred on the truth, x/div 1.30. A
#                likelihood, not a replacement: the draws are kept and
#                reweighted, so the AEM's own spread still shows through.
#   both         gamma and induction together.
#
# CLAY FRACTION vs CLAY TERM. The brief says "clay fraction prior". Stage 0 has
# no prior on clay fraction -- it has one log-uniform on the PRODUCT
# sigma_clay*C, because that product is the only combination entering
# sigma_bulk. The gamma log is therefore applied to the product. A gamma log
# constrains the fraction and not sigma_clay, so x/div 2.50 on the product is
# the optimistic reading of what a gamma log buys.
#
# CONDUCTIVITY vs CONDUCTANCE. Every Stage 3 window on both columns is exactly
# 16.00 m thick (Stage 3 slides a fixed-length window), so a likelihood on the
# window's mean conductivity and one on its conductance are the same function.
# It is applied to sigma_bulk, which is the quantity the truth is read on.
#
# NEGATIVE sigma_w IS A RESULT, NOT AN ERROR: the drawn clay term alone exceeds
# the measured conductance, so the data are consistent with arbitrarily fresh
# water. Counted as fresh in P(fresh), excluded from the quantiles, and the
# ACCEPTANCE RATE is the fraction above zero.
using JLD2, Statistics, Printf, Random
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "config.jl")); using .SiteConfig
const Z95 = 1.959963984540054

"""Weighted quantile, the StatsBase/Hyndman-Fan type 7 generalisation. With
equal weights it reduces EXACTLY to `Statistics.quantile`, which is what lets
the survey-only state reproduce the stored Stage 5 numbers to the last digit --
a nearest-order-statistic weighted quantile does not, and reported p5 816 as
815.87."""
function wquantile(v::AbstractVector, w::AbstractVector, q::Real)
    isempty(v) && return NaN
    pp = sortperm(v); vs = v[pp]; ws = w[pp]
    W = sum(ws); W <= 0 && return NaN
    S = cumsum(ws)
    h = q * (W - ws[1]) + ws[1]
    k = searchsortedfirst(S, h)
    k <= 1 && return vs[1]
    k > length(vs) && return vs[end]
    Sk1 = S[k-1]
    den = S[k] - Sk1
    den <= 0 ? vs[k] : vs[k-1] + (h - Sk1)/den * (vs[k] - vs[k-1])
end
ess(w) = (s = sum(w); s <= 0 ? 0.0 : s^2 / sum(w .^ 2))

"""log-normal with a 95% interval of x/div `fac` about `mu`."""
lognorm_draw(rng, mu, fac, n) = exp.(log(mu) .+ (log(fac)/Z95) .* randn(rng, n))
"""Likelihood kernel in `x` for a log-normal measurement of `x` at `mu`.
A likelihood, so no 1/x density factor."""
loglik(x, mu, fac) = (s = log(fac)/Z95; x <= 0 ? -Inf : -0.5*((log(x)-log(mu))/s)^2)

"""The window truths a simulated log would read, on the representative window."""
function window_truths(cfg)
    T = load(joinpath(ROOT, cfg["data"]["truth"]))
    G = load(joinpath(ROOT, cfg["output"]["dir"], "stage4_$(cfg["name"]).jld2"))
    w = (G["window"].top, G["window"].base); L = w[2] - w[1]
    z, dz = T["z"], T["dz"]
    k = findall(i -> z[i] > w[1] && z[i] <= w[2], eachindex(z))
    (window = w,
     sigma_bulk = sum(T["σ_bulk"][k]) * dz / L,
     cp         = sum(T["σ_s"][k] .* T["clay"][k]) * dz / L,
     sigma_w    = mean(T["σw"][k]))
end

"""One prior state on the stored Stage 5 draws."""
function runstate(S5, bt, MGL, Xlim; gamma = false, induction = false,
                  cp_fixed = nothing, rng = nothing)
    SB = S5["sigma_bulk"]; TH = S5["theta_draw"]; M = S5["m_draw"]
    n  = length(SB)
    CP = cp_fixed !== nothing ? fill(cp_fixed, n) :
         gamma ? lognorm_draw(rng, bt.cp, 2.50, n) : copy(S5["clay_product"])
    lw = induction ? [loglik(SB[i], bt.sigma_bulk, 1.30) for i in 1:n] : zeros(n)
    wt = exp.(lw .- maximum(lw))
    SWv = (SB .- CP) ./ TH .^ M
    TDS = SWv .* MGL
    tot = sum(wt)
    pos = findall(>(0), SWv)
    acc = sum(wt[pos]) / tot
    Pfresh = sum(wt[i] for i in 1:n if TDS[i] < Xlim; init = 0.0) / tot
    (p05 = wquantile(TDS[pos], wt[pos], 0.05), p50 = wquantile(TDS[pos], wt[pos], 0.50),
     p95 = wquantile(TDS[pos], wt[pos], 0.95),
     cp05 = wquantile(CP[pos], wt[pos], 0.05), cp50 = wquantile(CP[pos], wt[pos], 0.50),
     cp95 = wquantile(CP[pos], wt[pos], 0.95),
     accept = acc, Pfresh = Pfresh, ess = ess(wt),
     TDS = TDS[pos], wt = wt[pos])
end

verdict(r) = r.Pfresh >= 0.95 ? "PASS" : (r.accept < 0.10 ? "FRAGILE" : "UNDETERMINED")

const STATES = (("S0","survey only",false,false), ("G","gamma",true,false),
                ("I","induction",false,true),     ("GI","both",true,true))

function main()
    cfg2P = SiteConfig.resolve(joinpath(ROOT, "configs", "synth_2P.toml"))
    cfg2S = SiteConfig.resolve(joinpath(ROOT, "configs", "synth_2S.toml"))
    pr = cfg2P["priors"]; hy = cfg2P["hypothesis"]
    MGL = Float64(pr["mg_per_sm"]); Xlim = Float64(hy["salinity"]["value"])
    S5 = load(joinpath(ROOT, cfg2P["output"]["dir"], "stage5_synth_2P.jld2"))
    bt = window_truths(cfg2P)
    rng = MersenneTwister(20261103)

    println("SECTION 6 -- narrowing the Stage 0 priors with borehole logs.  test two (2P)")
    @printf("  Stage 5 draws reused: %s\n", S5["source_note"])
    @printf("  representative window %.2f-%.2f m\n", bt.window...)
    @printf("  truth on that window: sigma_bulk %.5f S/m, clay term %.5f S/m (%.0f%% of bulk),\n",
            bt.sigma_bulk, bt.cp, 100*bt.cp/bt.sigma_bulk)
    @printf("                        sigma_w %.5f S/m = %.0f mg/L\n", bt.sigma_w, bt.sigma_w*MGL)

    # ---- the check the brief demands, before anything else ----
    s0 = runstate(S5, bt, MGL, Xlim)
    ok = isapprox(s0.p50, S5["TDS_p50"]; atol = 1e-6) &&
         isapprox(s0.p05, S5["TDS_p5"];  atol = 1e-6) &&
         isapprox(s0.p95, S5["TDS_p95"]; atol = 1e-6) &&
         isapprox(s0.accept, S5["acceptance"]; atol = 1e-9)
    println("\nCHECK: survey-only against the Stage 5 result quoted in 5.2")
    @printf("    p50      %14.6f   stored %14.6f\n", s0.p50, S5["TDS_p50"])
    @printf("    p5       %14.6f   stored %14.6f\n", s0.p05, S5["TDS_p5"])
    @printf("    p95      %14.6f   stored %14.6f\n", s0.p95, S5["TDS_p95"])
    @printf("    accept   %14.6f   stored %14.6f\n", s0.accept, S5["acceptance"])
    ok || error("survey-only does NOT reproduce the stored Stage 5 result -- stopping")
    println("    reproduced exactly (same draws, same priors)")

    res = Dict{String,Any}(); order = String[]
    println("\n  " * "-"^92)
    @printf("  %-14s %22s %9s %8s %11s %8s  %s\n", "prior state", "TDS p50 (p5-p95) mg/L",
            "P(fresh)", "accept", "clay p50", "clay x", "verdict")
    println("  " * "-"^92)
    for (key, name, g, i) in STATES
        r = key == "S0" ? s0 : runstate(S5, bt, MGL, Xlim; gamma = g, induction = i, rng = rng)
        res[key] = r; push!(order, key)
        @printf("  %-14s %7.0f (%5.0f-%6.0f) %9.3f %7.1f%% %11.5f %8.1f  %s\n",
                name, r.p50, r.p05, r.p95, r.Pfresh, 100r.accept, r.cp50, r.cp95/r.cp05,
                verdict(r))
    end
    println("  " * "-"^92)
    @printf("  truth %.0f mg/L; fresh line %.0f mg/L; true clay term %.5f S/m\n",
            bt.sigma_w*MGL, Xlim, bt.cp)
    for (k,_,_,i) in STATES
        i && @printf("  ESS, %-12s %8.0f of %d draws (induction reweighting)\n",
                     k, res[k].ess, length(S5["sigma_bulk"]))
    end

    # ------------------------------------------------ 6.2 transfer test
    println("\nTRANSFER TEST (6.2). The clay term fixed at its TRUE value for the")
    println("window, used as the clay prior at a survey-only sounding. No induction.")
    println("Each column takes its own true clay term over its own window; the two")
    println("columns share one clay profile and differ only in their window.")
    xfer = Dict{String,Any}()
    @printf("\n  %-44s %9s %9s %8s\n", "case", "TDS p50", "P(fresh)", "truth")
    println("  " * "-"^74)
    for (tag, cfgN, nm) in (("2P", cfg2P, "fresh neighbour"), ("2S", cfg2S, "saline neighbour"))
        S5N = load(joinpath(ROOT, cfgN["output"]["dir"], "stage5_synth_$tag.jld2"))
        btN = window_truths(cfgN)
        for (ck, clab, cf) in (("prior", "Stage 0 clay prior", nothing),
                               ("fixed", "clay term fixed at truth", btN.cp))
            r = runstate(S5N, btN, MGL, Xlim; cp_fixed = cf)
            xfer["$(tag)_$(ck)"] = (r = r, truth = btN.sigma_w*MGL, label = "$nm, $clab",
                                    cp_true = btN.cp)
            @printf("  %-44s %9.0f %9.3f %8.0f\n", "$nm, $clab", r.p50, r.Pfresh,
                    btN.sigma_w*MGL)
        end
    end
    println("  " * "-"^74)

    odir = joinpath(ROOT, "outputs", "section6"); mkpath(odir)
    jldopen(joinpath(odir, "section6_2P.jld2"), "w") do f
        f["order"] = order
        f["names"] = Dict(k => n for (k,n,_,_) in STATES)
        f["truth_TDS"] = bt.sigma_w * MGL
        f["Xlim"] = Xlim
        f["window_truths"] = (sigma_bulk = bt.sigma_bulk, cp = bt.cp,
                              sigma_w = bt.sigma_w, window = bt.window)
        for k in order
            r = res[k]
            f["state_$k"] = (p05=r.p05, p50=r.p50, p95=r.p95, cp05=r.cp05, cp50=r.cp50,
                             cp95=r.cp95, accept=r.accept, Pfresh=r.Pfresh, ess=r.ess,
                             verdict=verdict(r))
        end
        for (k, v) in xfer
            f["xfer_$k"] = (p50=v.r.p50, p05=v.r.p05, p95=v.r.p95, Pfresh=v.r.Pfresh,
                            truth=v.truth, label=v.label, cp_true=v.cp_true,
                            TDS=v.r.TDS, wt=v.r.wt)
        end
    end
    println("\nwrote outputs/section6/section6_2P.jld2")
end
main()
