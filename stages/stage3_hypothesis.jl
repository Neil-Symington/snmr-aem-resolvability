# stage3_hypothesis.jl -- Stage 3 of the Section 3 workflow.
#
#   julia --project=. stage3_hypothesis.jl --config ../configs/synth_3P.toml
#
# Reads the hypothesis from the settings file as three optional parts --
# interval, water criterion, salinity criterion -- sweeps every window the
# interval spec allows, and reports the posterior on the window and the
# probability the water criterion is met.
#
# NO TRUTH IS READ HERE. The theta draws and their depth grid come from the
# SNMR ensemble; the hypothesis comes from the config. Nothing else is opened.
#
# TWO THINGS THIS STAGE GETS RIGHT THAT EARLIER CODE DID NOT
#
# 1. WINDOWS ARE NOT QUANTISED TO CELLS. The depth grid is log-growing, so a
#    cell whose centre lies inside a window does not lie inside it, and
#    summing whole cells makes a nominal 10.0 m window carry 8.8 m of
#    support. That error has already been found once in this project, in the
#    Stage 4 truth statistic. Here the integral is computed from a piecewise-
#    linear antiderivative evaluated at the window's actual edges, so
#    W(a,b) = F(b) - F(a) is exact for piecewise-constant theta and any real
#    a, b. The window edges are free, not grid points.
#
# 2. THE INTEGRAL, NOT A THRESHOLD. W = int theta dz over the window. The
#    integral is 8-10x better determined than theta at a point on these
#    synthetics and a median 2.2x on field data, and a pointwise threshold
#    pair moves P(aquifer) from 0.04 to 0.78 on the SAME draws -- which is a
#    statement about the threshold, not the aquifer.
#
# STAGE 3 NO LONGER PICKS A BEST WINDOW. For every candidate window it reports
# p_pass, the fraction of draws with W >= w, and then the INTERVALS those
# windows form: every window with p_pass >= 0.5, overlapping ones merged into
# one interval reported by its extent and its maximum p_pass, plus any secondary
# interval reaching p_pass >= 0.02 for the record. A single "best window" was
# never identified by the data anyway -- 3,093 of 5,709 windows tied at p_pass =
# 1.0000 on 3P, so the pick came from a tie-break, not from the measurement.
#
# The tie-break window is STILL COMPUTED and still written to the file, because
# Stages 4 and 5 have already run on it and changing it silently would make
# their outputs inconsistent with this one. It is reported as "the window Stages
# 4 and 5 used", and this stage now checks it lies inside a reported interval.
#
# AN AMBIGUITY IN THE EARLIER BRIEF, RESOLVED BY REPORTING BOTH. "Sweep all windows
# satisfying the interval spec on every draw; keep the window maximising
# P(criterion)" admits two readings, because a single draw has a W, not a P:
#
#   PER-DRAW  (reported as the primary result) for each draw, take the window
#     maximising that draw's W. This is the only reading that yields a
#     POSTERIOR on top, base and length, which is what the stage's output
#     list asks for. P(criterion) is then the fraction of draws whose own
#     best window clears the threshold.
#   FIXED-WINDOW  the single window maximising P(W >= w) across all draws,
#     with that P. One window, no posterior on its position.
#
# Both are computed and printed. They answer different questions: the
# per-draw posterior says where the water is if it is anywhere, the
# fixed-window P says how well one screen placed in advance would do.
#
# TWO DEGENERACIES IN THAT RULE, FOUND BY RUNNING IT, HANDLED EXPLICITLY
#
# A. MAXIMISING W IS SOLVED BY L = L_max, ALWAYS. theta >= 0, so W is
#    monotone non-decreasing in window length and the per-draw best window
#    is the longest one allowed, on every draw, on every column. The
#    resulting "posterior" on LENGTH is the bound L_max and nothing else --
#    measured at 100% of draws on all four tier-P columns. The top and base
#    posteriors ARE informative; the length is not, and is reported as
#    spec-determined rather than dressed up as an inference.
#
# B. A SATURATED CRITERION DOES NOT IDENTIFY A WINDOW. If P(W >= w) = 1 for
#    many windows, "the window maximising P" is whichever the loop happened
#    to see first -- scan order, not data. So ties are broken by a stated
#    rule (most robust, then shortest, then shallowest) and saturation is
#    detected and announced, because a criterion that cannot fail tells you
#    nothing about the column and must not be reported as a pass earned.
using JLD2, Statistics, Printf, TOML
const HERE = @__DIR__
include(joinpath(HERE, "..", "src", "config.jl")); using .SiteConfig

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) ? (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (i += 1)
    end
end
const CFGPATH = get(ARG, "config", "")
isempty(CFGPATH) && error("usage: stage3_hypothesis.jl --config <path/to/site.toml>")
const ROOT = normpath(joinpath(HERE, ".."))

# Cell edges from cell CENTRES, the convention the rest of this project uses
# (stage1_integral.jl `cellsof`): edges at the midpoints, first top at 0, last
# bottom extrapolated by the final half-spacing.
cellsof(z) = (mid = (z[1:end-1] .+ z[2:end]) ./ 2;
              top = vcat(0.0, mid);
              bot = vcat(mid, z[end] + (z[end] - z[end-1]) / 2);
              (top = top, bot = bot, dz = bot .- top))

"""Antiderivative of piecewise-constant theta, evaluated at arbitrary depths.
Returns F with F[k, j] = int_0^{d[k]} theta_j dz, exact for any d."""
function antideriv(Θ::AbstractMatrix, g, d::AbstractVector)
    nz, nd = size(Θ)
    cum = zeros(nz + 1, nd)                      # cum[i+1] = integral to bot[i]
    @inbounds for i in 1:nz, j in 1:nd
        cum[i+1, j] = cum[i, j] + Θ[i, j] * g.dz[i]
    end
    F = zeros(length(d), nd)
    @inbounds for (k, dk) in enumerate(d)
        i = searchsortedlast(g.bot, dk)          # cells fully above dk
        i = clamp(i, 0, nz)
        if i == nz
            for j in 1:nd; F[k, j] = cum[nz+1, j]; end
        else
            frac = (dk - (i == 0 ? 0.0 : g.bot[i]))   # part-way into cell i+1
            for j in 1:nd; F[k, j] = cum[i+1, j] + Θ[i+1, j] * frac; end
        end
    end
    F
end

function main()
    cfg = SiteConfig.resolve(CFGPATH)
    info = SiteConfig.assertvalid(cfg; root = ROOT)
    d, hy, out = cfg["data"], cfg["hypothesis"], cfg["output"]

    @printf("Stage 3 -- hypothesis.  site %s  (%s)\n", cfg["name"], cfg["label"])
    @printf("  lineage: column %s, parent %s, role %s\n",
            d["column"], isempty(info.parent) ? "-" : info.parent, info.role)
    @printf("  config:  %s\n", cfg["_source"])

    # ---- the SNMR ensemble. No truth file is opened. -------------------
    ep = get(d, "snmr_ensemble", "")
    isempty(ep) && error("""
        [data].snmr_ensemble is empty for column $(d["column"]).
        Stage 3 needs theta draws. Run Stage 2 for this column first
        (run_step7_pair.sh $(d["column"]) builds the pair), then set the path.""")
    epath = isfile(joinpath(ROOT, ep)) ? joinpath(ROOT, ep) : ep
    E = load(epath)
    Θ = Matrix{Float64}(reduce(hcat, E["ens"]))      # theta (v/v) directly
    z = E["zsmr"]; g = cellsof(z); nd = size(Θ, 2)
    @printf("  ensemble: %s\n", basename(epath))
    @printf("            %d draws over %d cells, %.2f-%.2f m, kernel %s, tier %s\n",
            nd, length(z), g.top[1], g.bot[end], E["kernel"], E["tier"])

    # ---- the hypothesis, as asserted ----------------------------------
    haskey(hy, "interval") || error("[hypothesis.interval] is required")
    iv = hy["interval"]
    step = Float64(get(iv, "search_step", 0.25))
    Lmin, Lmax = Float64(iv["L_min"]), Float64(iv["L_max"])
    dmin, dmax = Float64(iv["depth_min"]), Float64(iv["depth_max"])
    has_water = haskey(hy, "water")
    wtype = has_water ? get(hy["water"], "type", "integral") : ""
    wval  = has_water ? Float64(hy["water"]["value"]) : NaN

    println("\nHYPOTHESIS AS ASSERTED (pre-registered; not edited after the run)")
    @printf("  interval : %s, L %.1f-%.1f m, within %.1f-%.1f m, swept every %.2f m\n",
            iv["mode"], Lmin, Lmax, dmin, dmax, step)
    if has_water
        @printf("  water    : %s, %s\n", wtype,
                wtype == "integral" ? @sprintf("W = int theta dz >= %.2f m", wval) :
                                      @sprintf("theta > %.3f throughout", wval))
    else
        println("  water    : NOT ASSERTED -- reported as such, not as passed")
    end
    if haskey(hy, "salinity")
        @printf("  salinity : TDS %s %.0f mg/L   (gated in Stage 4, not here)\n",
                hy["salinity"]["direction"] == "less_than" ? "<" : ">",
                Float64(hy["salinity"]["value"]))
    else
        println("  salinity : NOT ASSERTED")
    end
    wtype == "integral" || error("only the integral water criterion is supported here; " *
                                 "pointwise is retained in stage1_integral.jl for old figures")

    # ---- sweep every allowed window -----------------------------------
    dgrid = collect(dmin:step:dmax)
    F = antideriv(Θ, g, dgrid)
    ng = length(dgrid)
    # (i, j) with dgrid[j] - dgrid[i] in [Lmin, Lmax]
    pairs = [(i, j) for i in 1:ng-1 for j in i+1:ng
             if Lmin - 1e-9 <= dgrid[j] - dgrid[i] <= Lmax + 1e-9]
    isempty(pairs) && error("no window satisfies the interval spec")
    @printf("\n  swept %d windows (%d tops x lengths %.2f-%.2f m) on each of %d draws\n",
            length(pairs), ng, Lmin, Lmax, nd)

    # Per-draw best window, and every window's P -- kept, not just the best,
    # because the SPREAD of P over windows is what tells you whether the
    # criterion discriminates at all.
    bestW = fill(-Inf, nd); besti = zeros(Int, nd); bestj = zeros(Int, nd)
    np = length(pairs)
    Pw = zeros(np); W5 = zeros(np)          # P and the 5th percentile of W
    wbuf = zeros(nd)
    for (q_, (i, j)) in enumerate(pairs)
        @inbounds begin
            nhit = 0
            for k in 1:nd
                w = F[j, k] - F[i, k]
                wbuf[k] = w
                if w > bestW[k]; bestW[k] = w; besti[k] = i; bestj[k] = j; end
                w >= wval && (nhit += 1)
            end
            Pw[q_] = nhit / nd
            W5[q_] = quantile(wbuf, 0.05)
        end
    end

    # TIE-BREAK, stated rather than left to scan order. Among windows whose P
    # is within one Monte Carlo standard error of the maximum, prefer the most
    # ROBUST (highest 5th percentile of W), then the SHORTEST (the most
    # conservative screen that still clears the criterion), then the
    # SHALLOWEST (cheapest to drill). Every step is a tie-break, not a model
    # choice: it only ever picks among windows the data cannot separate.
    Pmax = maximum(Pw)
    tol  = max(sqrt(max(Pmax*(1-Pmax), 0.0)/nd), 1/nd)
    cand = findall(>=(Pmax - tol), Pw)
    bykey = sort(cand, by = q_ -> (-W5[q_],
                                   dgrid[pairs[q_][2]] - dgrid[pairs[q_][1]],
                                   dgrid[pairs[q_][1]]))
    bq = bykey[1]; bPi, bPj = pairs[bq]; bestP = Pw[bq]
    n_tied = length(cand)
    saturated = Pmax >= 1 - 1/nd && n_tied > 1

    tops  = dgrid[besti]; bases = dgrid[bestj]; lens = bases .- tops

    # DOES THE CRITERION DISCRIMINATE? If the WORST window the spec allows
    # already clears the threshold on most draws, the criterion cannot fail
    # and a "PASS" on it is a property of the spec, not of the column.
    Wmin_any = minimum(W5); Wmax_any = maximum(W5)
    frac_windows_pass = count(>=(1 - 1/nd), Pw) / np
    P_perdraw = count(>=(wval), bestW) / nd
    se(p) = sqrt(max(p * (1 - p), 0.0) / nd)
    q(v, p) = quantile(v, p)

    println("\n" * "-"^70)
    println("PER-DRAW RESULT (primary) -- each draw's own best window")
    @printf("  top    p50 %6.2f m   p5-p95  %6.2f - %6.2f m\n", q(tops,.5),  q(tops,.05),  q(tops,.95))
    @printf("  base   p50 %6.2f m   p5-p95  %6.2f - %6.2f m\n", q(bases,.5), q(bases,.05), q(bases,.95))
    @printf("  length p50 %6.2f m   p5-p95  %6.2f - %6.2f m\n", q(lens,.5),  q(lens,.05),  q(lens,.95))
    @printf("  W      p50 %6.3f m   p5-p95  %6.3f - %6.3f m\n", q(bestW,.5), q(bestW,.05), q(bestW,.95))
    @printf("  P(W >= %.2f m) = %.4f   (Monte Carlo SE %.4f on %d draws)\n",
            wval, P_perdraw, se(P_perdraw), nd)
    # top/base multimodality: p50 with p5-p95 hides a two-aquifer column
    srt = sort(tops); gaps = diff(srt)
    bigg = findall(>(3.0), gaps)          # a >3 m hole in the support
    if !isempty(bigg)
        @printf("  WARNING: the window-top posterior is MULTIMODAL -- %d gap(s) > 3 m in its\n", length(bigg))
        println("           support. p50 with p5-p95 misrepresents it; the modes are:")
        edges = vcat(1, bigg .+ 1, length(srt) + 1)
        for b in 1:length(edges)-1
            seg = srt[edges[b]:edges[b+1]-1]
            @printf("           mode %d: %6.2f - %6.2f m, %5.1f%% of draws\n",
                    b, first(seg), last(seg), 100*length(seg)/nd)
        end
    end
    # A length posterior pinned at L_max is the spec biting, not the data.
    frac_at_max = count(>=(Lmax - step/2), lens) / nd
    frac_at_min = count(<=(Lmin + step/2), lens) / nd
    @printf("  length at L_max in %.1f%% of draws, at L_min in %.1f%%", 100frac_at_max, 100frac_at_min)
    println(frac_at_max > 0.5 ?
        "\n  NOTE: the length posterior is pinned at L_max -- the bound is carrying this,\n" *
        "        not the data. Widen L_max or report the window as unresolved in length." : "")

    println("-"^70)
    println("DOES THE WATER CRITERION DISCRIMINATE?")
    @printf("  P(W >= %.2f m) = 1 for %.1f%% of the %d allowed windows\n",
            wval, 100*frac_windows_pass, np)
    @printf("  the WORST allowed window still gives W(p5) = %.3f m; the best, %.3f m\n",
            Wmin_any, Wmax_any)
    if Wmin_any >= wval
        @printf("""
  SATURATED. Every window the interval spec allows clears W >= %.2f m on at
  least 95%% of draws, so this criterion CANNOT FAIL on this column and a
  PASS on it is a property of the spec, not a finding about the ground.
  Reported as PASS (uninformative) in the verdict table, not as PASS.
  To make it discriminating, w would have to exceed %.2f m -- but changing it
  now, after seeing this, is no longer a pre-registered threshold.
""", wval, Wmin_any)
    end
    # ---- INTERVALS, the stage's output -------------------------------
    tops_w  = [dgrid[pairs[q][1]] for q in 1:np]
    bases_w = [dgrid[pairs[q][2]] for q in 1:np]
    """Merge overlapping windows above `thr` into intervals; report extent and max p_pass."""
    function intervals(thr)
        k = findall(>=(thr), Pw)
        isempty(k) && return NTuple{3,Float64}[]
        o = sortperm(tops_w[k]); k = k[o]
        acc = NTuple{3,Float64}[]
        ct, cb, cp = tops_w[k[1]], bases_w[k[1]], Pw[k[1]]
        for q in k[2:end]
            if tops_w[q] <= cb            # overlaps the interval being built
                cb = max(cb, bases_w[q]); cp = max(cp, Pw[q])
            else
                push!(acc, (ct, cb, cp)); ct, cb, cp = tops_w[q], bases_w[q], Pw[q]
            end
        end
        push!(acc, (ct, cb, cp)); acc
    end
    """Merge windows whose p_pass lies in [lo, hi) -- so a SECONDARY interval is
    one that reached lo but never reached the primary threshold. Merging the
    whole >= 0.02 set instead returns a superset that swallows the primaries
    and says nothing."""
    function intervals_between(lo, hi)
        k = findall(q -> Pw[q] >= lo && Pw[q] < hi, 1:np)
        isempty(k) && return NTuple{3,Float64}[]
        o = sortperm(tops_w[k]); k = k[o]
        acc = NTuple{3,Float64}[]
        ct, cb, cp = tops_w[k[1]], bases_w[k[1]], Pw[k[1]]
        for q in k[2:end]
            if tops_w[q] <= cb
                cb = max(cb, bases_w[q]); cp = max(cp, Pw[q])
            else
                push!(acc, (ct, cb, cp)); ct, cb, cp = tops_w[q], bases_w[q], Pw[q]
            end
        end
        push!(acc, (ct, cb, cp)); acc
    end
    primary = intervals(0.5)
    # A SECONDARY interval is a region that reaches p_pass >= 0.02 but lies
    # OUTSIDE every primary. Merging the 0.02-to-0.5 band instead returns a
    # span that overlaps the primaries, because at any window top some lengths
    # pass and others do not -- that reported 8.25-60.00 m on every column and
    # said nothing. Subtracting the primary extents isolates the real second
    # region, which is what the confined column's deeper aquifer needs.
    # NOTE: the accumulator is `acc`, never `out` -- `out` is cfg["output"] in
    # the enclosing scope and a nested function assigning it silently replaces
    # it, which made jldsave fail with "invalid index: \"dir\"".
    function subtract_primaries(ivs)
        acc = NTuple{3,Float64}[]
        for (a, b, pp) in ivs
            segs = [(a, b)]
            for (pa, pb, _) in primary
                nxt = Tuple{Float64,Float64}[]
                for (sa, sb) in segs
                    sb <= pa || sa >= pb ? push!(nxt, (sa, sb)) :
                        (sa < pa && push!(nxt, (sa, pa)); sb > pb && push!(nxt, (pb, sb)))
                end
                segs = nxt
            end
            for (sa, sb) in segs
                sb - sa < 1e-9 && continue
                # windows that OVERLAP the leftover segment, not ones contained
                # in it: a segment can be shorter than L_min, and requiring
                # containment then returns NaN for a real region.
                k = findall(q -> tops_w[q] < sb - 1e-9 && bases_w[q] > sa + 1e-9, 1:np)
                push!(acc, (sa, sb, isempty(k) ? NaN : maximum(Pw[k])))
            end
        end
        acc
    end
    secondary = subtract_primaries(intervals(0.02))

    println("-"^70)
    @printf("INTERVALS  (p_pass = fraction of draws with W >= %.2f m, over %d windows)\n", wval, np)
    if isempty(primary)
        println("  no interval reaches p_pass >= 0.5")
    else
        println("  primary, p_pass >= 0.5:")
        for (a,b,pp) in primary
            @printf("    %6.2f - %6.2f m   extent %5.2f m   max p_pass %.4f\n", a, b, b-a, pp)
        end
    end
    if isempty(secondary)
        println("  no secondary interval reaches p_pass >= 0.02")
    else
        println("  secondary, p_pass >= 0.02 (for the record):")
        for (a,b,pp) in secondary
            @printf("    %6.2f - %6.2f m   extent %5.2f m   max p_pass %.4f\n", a, b, b-a, pp)
        end
    end

    println("-"^70)
    println("THE WINDOW STAGES 4 AND 5 USED (tie-break, retained for consistency)")
    @printf("  window %.2f - %.2f m  (length %.2f m)\n", dgrid[bPi], dgrid[bPj], dgrid[bPj]-dgrid[bPi])
    if saturated
        @printf("  NOT IDENTIFIED BY THE DATA: %d of %d windows tie at P = %.4f.\n",
                n_tied, np, Pmax)
        println("  This one comes from the stated tie-break (most robust, then shortest,")
        println("  then shallowest), NOT from the criterion.")
    elseif n_tied > 1
        @printf("  %d windows within one MC SE of the best; tie broken by the stated rule.\n", n_tied)
    end
    Wfix = F[bPj, :] .- F[bPi, :]
    @printf("  W      p50 %6.3f m   p5-p95  %6.3f - %6.3f m\n", q(Wfix,.5), q(Wfix,.05), q(Wfix,.95))
    @printf("  P(W >= %.2f m) = %.4f   (Monte Carlo SE %.4f)\n", wval, bestP, se(bestP))
    println("-"^70)

    inside = any(a -> dgrid[bPi] >= a[1] - 1e-9 && dgrid[bPj] <= a[2] + 1e-9, primary)
    @printf("  lies inside a primary interval: %s\n", inside ? "YES" : "NO")
    inside || println("  WARNING: Stages 4 and 5 ran on a window no reported interval contains.")

    odir = joinpath(ROOT, out["dir"]); mkpath(odir)
    jldsave(joinpath(odir, "stage3_$(cfg["name"]).jld2");
        column = d["column"], parent = get(d, "parent", ""), role = info.role,
        config_source = cfg["_source"], ensemble = basename(epath),
        m_prior = info.m_prior,
        interval_spec = (Lmin = Lmin, Lmax = Lmax, dmin = dmin, dmax = dmax, step = step),
        water_value = wval, ndraws = nd,
        tops = tops, bases = bases, lengths = lens, W = bestW,
        P_perdraw = P_perdraw, P_se = se(P_perdraw),
        fixed_window = (top = dgrid[bPi], base = dgrid[bPj]), W_fixed = Wfix,
        P_fixed = bestP, P_fixed_se = se(bestP),
        win_tops = tops_w, win_bases = bases_w, win_ppass = Pw,
        intervals_primary = isempty(primary) ? zeros(0,3) :
            reduce(vcat, [[a b pp] for (a,b,pp) in primary]),
        intervals_secondary = isempty(secondary) ? zeros(0,3) :
            reduce(vcat, [[a b pp] for (a,b,pp) in secondary]),
        n_tied = n_tied, n_windows = np, saturated = saturated,
        frac_windows_pass = frac_windows_pass,
        W5_min_any_window = Wmin_any, W5_max_any_window = Wmax_any,
        length_at_Lmax = frac_at_max, length_at_Lmin = frac_at_min)
    println("\nwrote ", joinpath(out["dir"], "stage3_$(cfg["name"]).jld2"))
end

main()
