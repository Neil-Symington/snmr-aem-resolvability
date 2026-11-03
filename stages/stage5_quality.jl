# stage5_quality.jl -- Stage 5. Water quality under the prior.
#
#   julia --project=. stage5_quality.jl --config ../configs/synth_3P.toml
#   julia --project=. stage5_quality.jl --config ... --source gate   # comparison only
#
# BULK CONDUCTIVITY COMES FROM THE STAGE 1 TRANS-D DRAWS, inside each draw's
# Stage 3 window. The gate's Laplace machinery enters Stage 4 and nowhere else.
#
# WHY, on the record: p(C) is a constrained-optimisation object built to locate
# an admissible BAND, not a posterior meant to be propagated. Using it put
# column 2P's Stage 5 median 2.31x above the truth (C p50 2.014 S against a
# true 0.8727 S) and pushed the accepted 90% interval off the truth entirely.
# The `--source gate` path is retained only long enough to document that
# comparison once, and is then removed.
#
#   sigma_bulk  = (1/L) * integral of sigma_j dz over the window   [Stage 1 draw j]
#   sigma_bulk  = sigma_w * theta_total^m  +  (sigma_clay * C)
#   => sigma_w  = (sigma_bulk - clay_product) / theta_total^m
#
# Each Monte Carlo sample pairs one Stage 1 conductivity draw with one Stage 3
# WINDOW draw. The two are independent posteriors given the data -- the windows
# come from the SNMR ensemble, the conductivities from the AEM ensemble -- so
# they are sampled independently, and the window's own thickness is used, not a
# fixed one.
#
# THE CLAY PRIOR IS ONE LOG-UNIFORM ON THE PRODUCT sigma_clay*C, in S/m: the
# only combination entering sigma_bulk and the only one the data could
# constrain. It replaced a logit-normal clay fraction times a log-uniform
# sigma_clay, which on 3P put 98.1% of its mass above the TOTAL measured
# conductance.
#
# NEGATIVE sigma_w IS A RESULT, NOT AN ERROR: the drawn clay term alone exceeds
# the observed conductance, so the data are consistent with arbitrarily fresh
# water. Counted, never clamped. The ACCEPTANCE RATE goes in the verdict table
# and below 10% the estimate is FRAGILE -- a statement about the prior's tail.
#
# NOTE ON THE INPUT FILE: Stage 1 draws are read from aem_draws_<tag>.jld2 and
# ONLY the keys `sigma` and `zall_aem` are loaded. aem_posterior.jld2 also
# stores `sigma_aem_truth`, and no stage after Stage 0 may be in a position to
# read a truth array even by accident.
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
isempty(CFGPATH) && error("usage: stage5_quality.jl --config <path> [--source stage1|gate]")
const SOURCE = get(ARG, "source", "stage1")
SOURCE in ("stage1", "gate") || error("--source must be stage1 or gate")

"""Cumulative integral of a piecewise-constant profile, evaluated at arbitrary
depths. `zb` are cell BOTTOMS; the first cell starts at 0."""
function antideriv(V::AbstractMatrix, zb::AbstractVector, d::AbstractVector)
    nz, nd = size(V)
    dz = diff(vcat(0.0, collect(Float64, zb)))
    cum = zeros(nz + 1, nd)
    @inbounds for i in 1:nz, j in 1:nd
        cum[i+1, j] = cum[i, j] + V[i, j] * dz[i]
    end
    F = zeros(length(d), nd)
    @inbounds for (k, dk) in enumerate(d)
        i = clamp(searchsortedlast(zb, dk), 0, nz)
        if i == nz
            for j in 1:nd; F[k, j] = cum[nz+1, j]; end
        else
            frac = dk - (i == 0 ? 0.0 : zb[i])
            for j in 1:nd; F[k, j] = cum[i+1, j] + V[i+1, j] * frac; end
        end
    end
    F
end

function main()
    cfg  = SiteConfig.resolve(CFGPATH)
    info = SiteConfig.assertvalid(cfg; root = ROOT)
    d, pr, hy, out = cfg["data"], cfg["priors"], cfg["hypothesis"], cfg["output"]
    odir = joinpath(ROOT, out["dir"])
    S3 = load(joinpath(odir, "stage3_$(cfg["name"]).jld2"))
    P4 = load(joinpath(odir, "stage4post_$(cfg["name"]).jld2"))

    @printf("Stage 5 -- water quality under the prior.  site %s   [source: %s]\n",
            cfg["name"], SOURCE == "stage1" ? "STAGE 1 trans-D draws" : "gate p(C)  (COMPARISON ONLY)")
    @printf("  lineage: column %s, parent %s, role %s\n",
            d["column"], isempty(info.parent) ? "-" : info.parent, info.role)

    MGL = Float64(pr["mg_per_sm"]); Xlim = Float64(hy["salinity"]["value"])
    fresh_is_less = hy["salinity"]["direction"] == "less_than"
    mlo, mhi = info.m_prior
    plo, phi_ = Float64.(pr["porosity_range"])
    cplo, cphi = Float64.(pr["clay"]["product"])

    println("\nPRIORS USED (printed beside the answer)")
    @printf("  m             ~ U(%.2f, %.2f)          [lithology %s]\n", mlo, mhi,
            get(pr, "m_lithology", "explicit"))
    @printf("  theta_total   ~ U(%.3f, %.3f)        [formation range, declared]\n", plo, phi_)
    @printf("  clay product  ~ log-uniform(%.1e, %.1e) S/m   [sigma_clay * C]\n", cplo, cphi)

    rng = MersenneTwister(20261101)
    NS = 40000

    # ---- sigma_bulk, per sample ---------------------------------------
    local SB, thks, srcnote
    if SOURCE == "stage1"
        sp = d["stage1_draws"]
        spath = isfile(joinpath(ROOT, sp)) ? joinpath(ROOT, sp) : sp
        Σ = Matrix{Float64}(load(spath, "sigma"))        # only these two keys
        zb = collect(Float64, load(spath, "zall_aem"))
        naem = size(Σ, 2)
        tops, bases = S3["tops"], S3["bases"]
        nwin = length(tops)
        @printf("  Stage 1 draws: %s, %d draws over %d cells to %.1f m\n",
                basename(spath), naem, size(Σ,1), last(zb))
        @printf("  Stage 3 windows: %d draws, top p50 %.2f m, base p50 %.2f m\n",
                nwin, median(tops), median(bases))
        # pair one conductivity draw with one window draw per sample
        ja = rand(rng, 1:naem, NS); jw = rand(rng, 1:nwin, NS)
        # evaluate the antiderivative once on the union of window edges
        edges = sort(unique(vcat(tops, bases)))
        F = antideriv(Σ, zb, edges)
        pos_of = Dict(e => i for (i, e) in enumerate(edges))
        SB = Vector{Float64}(undef, NS); thks = Vector{Float64}(undef, NS)
        @inbounds for s in 1:NS
            w = jw[s]; a, b = tops[w], bases[w]
            ia, ib = pos_of[a], pos_of[b]
            L = b - a
            SB[s] = (F[ib, ja[s]] - F[ia, ja[s]]) / L
            thks[s] = L
        end
        srcnote = @sprintf("Stage 1 draws x Stage 3 windows (%d x %d)", naem, nwin)
        @printf("  sigma_bulk over the window: p50 %.5f S/m, p5-p95 %.5f-%.5f S/m\n",
                median(SB), quantile(SB,.05), quantile(SB,.95))
    else
        G4 = load(joinpath(odir, "stage4_$(cfg["name"]).jld2"))
        thick = G4["thick"]; Cg = G4["Cg"]; lhat = G4["lambda_hat"]
        SW = G4["sweeps"][string(lhat)]
        lL = -SW.phi ./ 2 .- SW.ldH ./ 2 .+ SW.ldP ./ 2
        lL[.!SW.good] .= NaN
        k = findall(isfinite, lL)
        isempty(k) && error("no admissible node at lambda-hat; the gate path has no p(C)")
        lL = lL .- maximum(lL[k]); w = zeros(length(lL)); w[k] = exp.(lL[k])
        wk = w[k] ./ sum(w[k]); cdfC = cumsum(wk)
        drawC(u) = Cg[k[clamp(searchsortedfirst(cdfC, u), 1, length(k))]]
        SB = [drawC(rand(rng)) / thick for _ in 1:NS]
        thks = fill(thick, NS)
        srcnote = @sprintf("gate p(C) at lambda-hat = %.4g (COMPARISON ONLY)", lhat)
        @printf("  sigma_bulk from p(C): p50 %.5f S/m, p5-p95 %.5f-%.5f S/m\n",
                median(SB), quantile(SB,.05), quantile(SB,.95))
    end

    # ---- the three priors ---------------------------------------------
    cps = exp.(log(cplo) .+ (log(cphi) - log(cplo)) .* rand(rng, NS))
    ms  = mlo .+ (mhi - mlo) .* rand(rng, NS)
    ths = plo .+ (phi_ - plo) .* rand(rng, NS)
    @printf("  clay product drawn: p50 %.4g, p10-p90 %.4g-%.4g S/m\n",
            median(cps), quantile(cps,.1), quantile(cps,.9))

    sigw(sb, cp, th, m) = (sb - cp) / th^m
    SW_ = sigw.(SB, cps, ths, ms)
    neg = count(<=(0), SW_); accept = 1 - neg/NS
    TDS = SW_ .* MGL
    pos = findall(>(0), SW_)

    println("\nSIGMA_W POSTERIOR (clay term subtracted)")
    @printf("  ACCEPTANCE RATE (sigma_w > 0): %.1f%%  %s\n", 100accept,
            accept < 0.10 ? "<-- BELOW 10%: Stage 5 estimate is FRAGILE" : "")
    if neg > 0
        @printf("  %d of %d draws negative: the drawn clay term alone exceeds the measured\n", neg, NS)
        println("  conductance -- those draws are consistent with ARBITRARILY FRESH water.")
    end
    if !isempty(pos)
        @printf("  sigma_w (accepted) p50 %.4f S/m, p5-p95 %.4f-%.4f S/m\n",
                median(SW_[pos]), quantile(SW_[pos],.05), quantile(SW_[pos],.95))
        @printf("  TDS     (accepted) p50 %.0f mg/L, p5-p95 %.0f-%.0f mg/L\n",
                median(TDS[pos]), quantile(TDS[pos],.05), quantile(TDS[pos],.95))
    end
    Pcrit = fresh_is_less ? count(t -> t < Xlim, TDS)/NS : count(t -> t > Xlim, TDS)/NS
    @printf("  P(TDS %s %.0f mg/L) = %.4f   (MC SE %.4f, over ALL draws incl. negatives)\n",
            fresh_is_less ? "<" : ">", Xlim, Pcrit, sqrt(Pcrit*(1-Pcrit)/NS))
    @printf("  Stage 4's clay-free bound, for comparison: TDS p95 = %.0f mg/L, P(too salty) = %.4f\n",
            P4["tds_data_only_p95"], P4["P_wrong"])

    # ---- SWING TABLE, baseline = median of the ACCEPTED draws ----------
    # The all-inputs median was not usable as a baseline: on 3P the median clay
    # product (0.0102 S/m) exceeded the median sigma_bulk, so the baseline was
    # itself a REJECTED draw and the m and porosity rows read "negative" at both
    # ends -- uninformative by construction. The baseline is now the median of
    # the accepted draws, and the p10/p90 ranges are taken over those too.
    println("\nSWING TABLE   (one input at p10 / p90, all others at the ACCEPTED-draw median)")
    if isempty(pos)
        println("  no accepted draws: no baseline exists, so no swing table")
        swing = Dict{String,Any}(); flips = String[]
    else
        bSB, bcp, bth, bm = median(SB[pos]), median(cps[pos]), median(ths[pos]), median(ms[pos])
        @printf("  baseline: sigma_bulk %.5f S/m, clay %.5f S/m, theta %.3f, m %.3f -> TDS %.0f mg/L\n",
                bSB, bcp, bth, bm, sigw(bSB, bcp, bth, bm) * MGL)
        tds1(sb, cp, th, m) = sigw(sb, cp, th, m) * MGL
        side(t) = fresh_is_less ? ((t > 0 && t < Xlim) ? 1 : 0) : ((t > Xlim) ? 1 : 0)
        q10(v) = quantile(v, .10); q90(v) = quantile(v, .90)
        inputs = Any[
          ("sigma_bulk",   q10(SB[pos]),  q90(SB[pos]),  (x) -> tds1(x, bcp, bth, bm)),
          ("clay product", q10(cps[pos]), q90(cps[pos]), (x) -> tds1(bSB, x, bth, bm)),
          ("m",            q10(ms[pos]),  q90(ms[pos]),  (x) -> tds1(bSB, bcp, bth, x)),
          ("porosity",     q10(ths[pos]), q90(ths[pos]), (x) -> tds1(bSB, bcp, x, bm)),
        ]
        @printf("  %-14s %-11s %-11s %-14s %-14s %s\n",
                "input", "p10", "p90", "TDS at p10", "TDS at p90", "P range")
        swing = Dict{String,Any}(); flips = String[]
        for (nm, a, b, f) in inputs
            ta, tb = f(a), f(b); pa, pb = side(ta), side(tb)
            fmt(t) = t <= 0 ? "neg (clay>sb)" : @sprintf("%.0f mg/L", t)
            @printf("  %-14s %-11.4g %-11.4g %-14s %-14s %d -> %d%s\n",
                    nm, a, b, fmt(ta), fmt(tb), pa, pb,
                    pa == pb ? "" : "   <-- FLIPS the verdict alone")
            swing[nm] = Dict("p10"=>a, "p90"=>b, "tds_p10"=>ta, "tds_p90"=>tb,
                             "P_p10"=>pa, "P_p90"=>pb)
            pa == pb || push!(flips, nm)
        end
        println(isempty(flips) ? "  no single input flips the verdict on its own" :
                "  inputs that flip the verdict alone: " * join(flips, ", "))
    end

    tag = SOURCE == "stage1" ? "" : "_gatesrc"
    jldsave(joinpath(odir, "stage5_$(cfg["name"])$(tag).jld2");
        column = d["column"], parent = get(d,"parent",""), role = info.role,
        source = SOURCE, source_note = srcnote,
        window_p50 = (top = median(S3["tops"]), base = median(S3["bases"])),
        thickness_p50 = median(thks),
        sigma_bulk = SB, sigma_w = SW_, TDS = TDS, frac_negative = neg/NS,
        # the clay product per draw, needed for the (sigma_clay, TDS) trade-off
        clay_product = cps, m_draw = ms, theta_draw = ths,
        sigma_bulk_p50 = median(SB),
        sigma_w_p50 = isempty(pos) ? NaN : median(SW_[pos]),
        TDS_p50 = isempty(pos) ? NaN : median(TDS[pos]),
        TDS_p5  = isempty(pos) ? NaN : quantile(TDS[pos],.05),
        TDS_p95 = isempty(pos) ? NaN : quantile(TDS[pos],.95),
        P_criterion = Pcrit, X_limit = Xlim,
        priors = (m = [mlo,mhi], porosity = [plo,phi_], clay_product = [cplo,cphi]),
        acceptance = accept, fragile = accept < 0.10,
        swing = swing, flips = flips, undetermined = false)
    println("\nwrote ", relpath(joinpath(odir, "stage5_$(cfg["name"])$(tag).jld2"), ROOT))
end
main()
