# stage4_post.jl -- Stage 4's two outputs, from the saved sweep.
#
#   julia --project=. stage4_post.jl --config ../configs/synth_3P.toml
#
# Deliberately separate from stage4_gate_v2.jl. The sweep costs hours (one
# Jacobian is 52x a forward model and each node needs tens of them); the band,
# p(C), the salinity conversion and Delta cost seconds. Keeping them apart means
# a mistake in the cheap half does not cost another run of the expensive half --
# which has already happened in this project more than once.
#
# (a) PROFILE BAND. chi2(C) <= N + 3.84, referenced to the TARGET misfit N
#     rather than to the best fit found. Edges by bisection on the interpolated
#     curve. Each edge carries an OPEN flag: an edge at the limit of the grid is
#     not a measurement, and on 19 field soundings every lower edge ran to the
#     sweep floor while being read as physics.
#     The upper edge converts to a CLAY-FREE salinity -- the DATA-ONLY BOUND.
#
# (b) BAYESIAN CURVE. p(C) from the Laplace log L under a flat-in-log-C prior.
#     P(TDS on the wrong side of X) integrates theta over the SNMR draws and m
#     over its prior. Delta is the MAXIMUM change in that P across lambda-hat/10,
#     lambda-hat, lambda-hat x10 and across dropping one node at each grid end.
#
# WHY THE CONVERSION IS CLAY-FREE. sigma_bulk = sigma_w*theta^m + sigma_s*clay,
# and the AEM cannot separate the two terms. Dropping the clay term attributes
# ALL the conductance to the pore water, so the recovered sigma_w is an
# OVER-estimate and P(too salty) is an over-estimate too. That is the
# conservative direction, and it is the only reason the number is reportable
# without knowing the clay. Archie with no surface term read as an estimate
# rather than a bound FAILS on clay-rich ground -- it reads clay as salt.
using JLD2, Statistics, Printf, LinearAlgebra, Random
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "config.jl")); using .SiteConfig

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) ? (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (i += 1)
    end
end
const CFGPATH = get(ARG, "config", "")
isempty(CFGPATH) && error("usage: stage4_post.jl --config <path>")

cellsof(z) = (mid = (z[1:end-1] .+ z[2:end]) ./ 2;
              top = vcat(0.0, mid);
              bot = vcat(mid, z[end] + (z[end] - z[end-1]) / 2);
              (top = top, bot = bot, dz = bot .- top))

"""Overlap-truncated cell weights for the window [a,b] on a grid of centres."""
function winweights(z, a, b)
    g = cellsof(z)
    w = [max(0.0, min(g.bot[i], b) - max(g.top[i], a)) for i in eachindex(z)]
    idx = findall(>(1e-12), w)
    (idx = idx, w = w[idx])
end

"""Linear interpolation of y(x) at xq, NaN outside the finite support."""
function interp(x, y, xq)
    k = findall(isfinite, y); isempty(k) && return NaN
    (xq < x[k[1]] || xq > x[k[end]]) && return NaN
    j = searchsortedlast(x, xq); j = clamp(j, 1, length(x) - 1)
    (isfinite(y[j]) && isfinite(y[j+1])) || return NaN
    t = (xq - x[j]) / (x[j+1] - x[j])
    (1 - t) * y[j] + t * y[j+1]
end

"""Bisect for f(x) = target between bracketing grid points. Returns the edge
and whether it is OPEN (the crossing was never found inside the grid)."""
function edge(x, y, target, from_left)
    k = findall(isfinite, y)
    length(k) < 2 && return (NaN, true)
    xs, ys = x[k], y[k]
    below = ys .<= target
    any(below) || return (NaN, true)
    if from_left
        i = findfirst(below)
        i == 1 && return (xs[1], true)                 # open: already inside at the floor
        lo, hi = xs[i-1], xs[i]
    else
        i = findlast(below)
        i == length(xs) && return (xs[end], true)      # open: still inside at the ceiling
        lo, hi = xs[i], xs[i+1]
    end
    for _ in 1:60
        mid = (lo + hi) / 2
        v = interp(xs, ys, mid)
        isfinite(v) || break
        (v <= target) == from_left ? (hi = mid) : (lo = mid)
        abs(hi - lo) < 1e-6 && break
    end
    ((lo + hi) / 2, false)
end

function main()
    cfg  = SiteConfig.resolve(CFGPATH)
    info = SiteConfig.assertvalid(cfg; root = ROOT)
    d, g, pr, hy, out = cfg["data"], cfg["gate"], cfg["priors"], cfg["hypothesis"], cfg["output"]
    odir = joinpath(ROOT, out["dir"])
    G  = load(joinpath(odir, "stage4_$(cfg["name"]).jld2"))
    S3 = load(joinpath(odir, "stage3_$(cfg["name"]).jld2"))
    A  = load(joinpath(odir, "stage0_$(cfg["name"]).jld2"))

    @printf("Stage 4 outputs -- site %s  (%s)\n", cfg["name"], cfg["label"])
    @printf("  lineage: column %s, parent %s, role %s\n",
            d["column"], isempty(info.parent) ? "-" : info.parent, info.role)

    lt, Cg, thick = G["lt"], G["Cg"], G["thick"]
    lhat = G["lambda_hat"]; ndata = G["ndata"]
    SWK = G["sweeps"]; keysf = sort(parse.(Float64, collect(keys(SWK))))
    wtop, wbot = G["window"].top, G["window"].base
    @printf("  window %.2f-%.2f m (%.3f m), lambda-hat %.4g, N = %d\n", wtop, wbot, thick, lhat, ndata)

    # ---- theta^m integral over the window, from the SNMR draws ---------
    ep = A["snmr_ensemble"]
    epath = isfile(joinpath(ROOT, ep)) ? joinpath(ROOT, ep) : ep
    E = load(epath); Θ = Matrix{Float64}(reduce(hcat, E["ens"])); zs = E["zsmr"]
    ww = winweights(zs, wtop, wbot)
    @printf("  SNMR: %d draws; window support on the SNMR grid %.4f m (window %.4f m)\n",
            size(Θ, 2), sum(ww.w), wbot - wtop)
    mlo, mhi = info.m_prior
    MGL = Float64(pr["mg_per_sm"])
    Xlim = Float64(hy["salinity"]["value"])
    wrongside_gt = hy["salinity"]["direction"] == "less_than"   # "too salty" = TDS > X

    rng = MersenneTwister(4242)
    NS = 20000
    jd = rand(rng, 1:size(Θ, 2), NS)
    ms = mlo .+ (mhi - mlo) .* rand(rng, NS)
    # K_s = int theta^m dz over the window, one per sample
    K = [sum((@view Θ[ww.idx, jd[s]]) .^ ms[s] .* ww.w) for s in 1:NS]
    @printf("  int theta^m dz over the window: p50 %.4f m, p5-p95 %.4f-%.4f  (m ~ U(%.2f,%.2f))\n",
            median(K), quantile(K, .05), quantile(K, .95), mlo, mhi)

    # ---- (a) the profile band -----------------------------------------
    # The misfit profile: the LOWEST lambda2 swept, i.e. the least prior-shaped
    # curve available under the same constraint. Stated, not assumed.
    l2prof = minimum(keysf)
    SP = SWK[string(l2prof)]
    chi = copy(SP.chi2); chi[.!SP.good] .= NaN
    target = ndata + Float64(g["band"]["delta_chi2"])

    # HOW MANY NODES SURVIVED. The grid runs to the model bounds by design and
    # nodes near them cannot satisfy the constraint, so a sizeable exclusion is
    # expected -- but if NOTHING survives there is no band and no p(C), and that
    # must be reported as an UNDETERMINED result rather than crashing on an
    # empty reduction. A 3-node grid over 4 decades excludes all three.
    println("\nADMISSIBILITY (nodes excluded, never normalised over)")
    for l2 in keysf
        SWv = SWK[string(l2)]
        @printf("    lambda2 %-10.4g %3d of %3d admissible, %3d on a bound, %3d off-constraint\n",
                l2, count(SWv.good), length(SWv.good), count(SWv.nbound .> 0),
                count(SWv.conres .>= 1e-8))
    end
    nadm = count(SP.good)
    if nadm < 2
        println("""
  NO ADMISSIBLE BAND. Fewer than two nodes on the grid satisfied the
  constraint with an interior optimum, so chi2(C) is not defined anywhere the
  Laplace approximation is valid. This is an UNDETERMINED result, not a wide
  one: the grid's reachable region was not sampled. Either the grid is too
  coarse for this column (nodes 0.1 dex apart are needed near the optimum) or
  the column's conductance lies outside the model bounds.""")
        jldsave(joinpath(odir, "stage4post_$(cfg["name"]).jld2");
            column = d["column"], parent = get(d,"parent",""), role = info.role,
            window = (top = wtop, base = wbot), thick = thick, lambda_hat = lhat,
            band_C = (lo = NaN, hi = NaN), band_dex = NaN,
            band_open = (lo = true, hi = true), chi2_target = target,
            profile_lambda2 = l2prof, tds_data_only_p95 = NaN, tds_best_case = NaN,
            P_wrong = NaN, P_by_lambda2 = Dict{Float64,Float64}(), P_drop_ends = NaN,
            delta_lambda2 = NaN, delta_drop = NaN, delta = NaN,
            delta_tol = Float64(g["bayes"]["delta_tol"]),
            verdict_label = "UNDETERMINED (no admissible node)",
            K_theta_m = (p50 = median(K), p5 = quantile(K,.05), p95 = quantile(K,.95)),
            m_prior = collect(info.m_prior), X_limit = Xlim, mg_per_sm = MGL,
            n_admissible = nadm)
        println("\nwrote ", relpath(joinpath(odir, "stage4post_$(cfg["name"]).jld2"), ROOT),
                "  (UNDETERMINED)")
        return
    end
    # Prefer the edges Stage 4 located by TRUE BISECTION (solving at new C),
    # falling back to interpolation between grid nodes for runs made before the
    # bisection existed. With a 15-node fine grid the interpolated edge is too
    # coarse to quote, so this is not a cosmetic preference.
    if haskey(G, "edge_lo_dex") && isfinite(G["edge_lo_dex"])
        elo, open_lo = G["edge_lo_dex"], G["edge_lo_open"]
        ehi, open_hi = G["edge_hi_dex"], G["edge_hi_open"]
        @printf("    edges from Stage 4 bisection (%d forward solves)\n", get(G,"edge_solves",0))
    else
        (elo, open_lo) = edge(lt, chi, target, true)
        (ehi, open_hi) = edge(lt, chi, target, false)
        println("    edges interpolated between grid nodes (no bisection in this run)")
    end
    Clo = isfinite(elo) ? 10^elo * thick : NaN
    Chi = isfinite(ehi) ? 10^ehi * thick : NaN
    println("\n(a) PROFILE BAND   chi2 <= N + ", g["band"]["delta_chi2"], " = ", round(target, digits=2),
            "   [profile at lambda2 = ", round(l2prof, sigdigits=3), "]")
    fchi = filter(isfinite, chi)
    @printf("    min chi2 on the grid %.2f (/N %.3f) at sigma_eff %.4g S/m\n",
            minimum(fchi), minimum(fchi)/ndata, 10^lt[argmin(replace(chi, NaN=>Inf))])
    @printf("    band  C  %.4g - %.4g S   (sigma_eff %.4g - %.4g S/m)\n",
            Clo, Chi, isfinite(elo) ? 10^elo : NaN, isfinite(ehi) ? 10^ehi : NaN)
    @printf("    width %.3f dex   lower edge %s   upper edge %s\n",
            (isfinite(elo) && isfinite(ehi)) ? ehi - elo : NaN,
            open_lo ? "OPEN (at the grid floor -- NOT a measurement)" : "closed",
            open_hi ? "OPEN (at the grid ceiling -- NOT a measurement)" : "closed")
    # Under the two-pass grid an OPEN edge carries an extra caveat: the fine
    # grid was bracketed by the coarse pass, so a reachable region can lie up to
    # one COARSE step beyond the fine grid's end and never be sampled. The edge
    # is then not merely open, it is open by at least that much.
    cstep = get(G, "coarse_step_dex", NaN)
    if get(G, "grid_mode", "") == "two_pass" && (open_lo || open_hi) && isfinite(cstep)
        @printf("    two-pass caveat: an open edge may lie up to one COARSE step (%.3f dex)\n", cstep)
        println("      beyond the fine grid, which was never sampled. Widen n_nodes_coarse")
        println("      to close it, or report the edge as open by at least that margin.")
    end
    # the data-only bound: the upper edge, clay-free, through theta and m
    tds_bound = NaN
    if isfinite(Chi)
        sw = Chi ./ K
        tds_bound = quantile(sw .* MGL, 0.95)
        @printf("    DATA-ONLY BOUND: upper edge clay-free -> TDS p95 = %.0f mg/L\n", tds_bound)
        @printf("       (p50 %.0f, p5 %.0f; clay-free, so this is an UPPER bound on salinity)\n",
                quantile(sw .* MGL, .5), quantile(sw .* MGL, .05))
    else
        println("    DATA-ONLY BOUND: not defined -- the upper edge is open")
    end

    # ---- (b) the Bayesian curve ---------------------------------------
    function pofC(SWv)
        lL = -SWv.phi ./ 2 .- SWv.ldH ./ 2 .+ SWv.ldP ./ 2
        lL[.!SWv.good] .= NaN
        k = findall(isfinite, lL); isempty(k) && return (nothing, nothing)
        lL = lL .- maximum(lL[k])
        w = zeros(length(lL)); w[k] = exp.(lL[k])      # flat in log C
        s = sum(w[k]) * (lt[2] - lt[1]); w ./= s
        (w, k)
    end
    """P(TDS on the wrong side of Xlim), integrating C over p(C) and theta,m over draws."""
    function Pwrong(SWv; drop_ends = 0)
        w, k = pofC(SWv); w === nothing && return NaN
        kk = collect(k)
        if drop_ends > 0 && length(kk) > 2*drop_ends
            kk = kk[1+drop_ends:end-drop_ends]
        end
        isempty(kk) && return NaN
        ww2 = w[kk] ./ sum(w[kk])
        # sample C from the discrete p(C), pair with a theta,m sample
        cdf = cumsum(ww2)
        rngl = MersenneTwister(99)
        tot = 0.0
        for s in 1:NS
            u = rand(rngl); i = kk[searchsortedfirst(cdf, u)]
            swv = Cg[i] / K[s] * MGL
            tot += wrongside_gt ? (swv > Xlim) : (swv < Xlim)
        end
        tot / NS
    end

    println("\n(b) BAYESIAN CURVE   Laplace log L, p(C) flat in log C")
    @printf("    %-12s %-10s %-12s %-10s %s\n", "lambda2", "usable", "P(TDS>X)", "C p50 (S)", "note")
    Ps = Dict{Float64,Float64}()
    for l2 in keysf
        SWv = SWK[string(l2)]
        w, k = pofC(SWv)
        P = Pwrong(SWv); Ps[l2] = P
        cm = if w === nothing; NaN else
            cdf = cumsum(w[k]) ./ sum(w[k]); Cg[k[searchsortedfirst(cdf, 0.5)]]
        end
        @printf("    %-12.4g %-10d %-12.4f %-10.4g %s\n", l2, count(SWv.good), P, cm,
                l2 == lhat ? "<- lambda-hat" : "")
    end
    Pbest = Ps[lhat]

    # Delta: lambda2 spread AND dropping one node at each grid end
    Pdrop = Pwrong(SWK[string(lhat)]; drop_ends = Int(g["bayes"]["delta_drop_end_nodes"]))
    cand = [Ps[l] for l in keysf if isfinite(Ps[l])]
    d_lam  = (isempty(cand) || !isfinite(Pbest)) ? NaN : maximum(abs.(cand .- Pbest))
    d_drop = (isfinite(Pdrop) && isfinite(Pbest)) ? abs(Pdrop - Pbest) : NaN
    dfin = filter(isfinite, [d_lam, d_drop])
    Δ = isempty(dfin) ? NaN : maximum(dfin)
    dtol = Float64(g["bayes"]["delta_tol"])
    label = !isfinite(Δ) ? "UNDETERMINED (Delta not computable)" :
            Δ <= dtol ? "DATA-CARRIED" : "PRIOR-CARRIED"
    @printf("\n    P(TDS > %.0f mg/L) at lambda-hat = %.4f\n", Xlim, Pbest)
    @printf("    Delta from lambda2 (x1/10, x10)        = %.4f\n", d_lam)
    @printf("    Delta from dropping %d node at each end = %.4f\n",
            Int(g["bayes"]["delta_drop_end_nodes"]), d_drop)
    @printf("    Delta = max = %.4f   vs tolerance %.2f   ->  %s\n", Δ, dtol, label)

    # best-case salinity for Section 7.5: median C, no clay
    w, k = pofC(SWK[string(lhat)])
    Cmed = if w === nothing; NaN else
        cdf = cumsum(w[k]) ./ sum(w[k]); Cg[k[searchsortedfirst(cdf, 0.5)]]
    end
    tds_best = isfinite(Cmed) ? median(Cmed ./ K .* MGL) : NaN
    @printf("\n    BEST-CASE SALINITY (median C, no clay) = %.0f mg/L   [Section 7.5]\n", tds_best)

    jldsave(joinpath(odir, "stage4post_$(cfg["name"]).jld2");
        column = d["column"], parent = get(d,"parent",""), role = info.role,
        window = (top = wtop, base = wbot), thick = thick, lambda_hat = lhat,
        band_C = (lo = Clo, hi = Chi), band_dex = (isfinite(elo)&&isfinite(ehi)) ? ehi-elo : NaN,
        band_open = (lo = open_lo, hi = open_hi), chi2_target = target,
        profile_lambda2 = l2prof,
        tds_data_only_p95 = tds_bound, tds_best_case = tds_best,
        P_wrong = Pbest, P_by_lambda2 = Ps, P_drop_ends = Pdrop,
        delta_lambda2 = d_lam, delta_drop = d_drop, delta = Δ,
        delta_tol = dtol, verdict_label = label,
        K_theta_m = (p50 = median(K), p5 = quantile(K,.05), p95 = quantile(K,.95)),
        m_prior = collect(info.m_prior), X_limit = Xlim, mg_per_sm = MGL)
    println("\nwrote ", relpath(joinpath(odir, "stage4post_$(cfg["name"]).jld2"), ROOT))
end
main()
