# fig07_gate_sigmin1e-3.jl -- Figure 7, two panels on one conductance axis.
#
# (a) chi2 of the forced fit against window conductance, fresh columns only.
# (b) the bar strip: the range of window conductance the data admit.
#
# The file name carries sigma_min = 1e-3 S/m because that is the run these
# edges come from. The sigma_min = 1e-4 re-run has NOT been made; the only
# 1e-4 work is the targeted open-edge walk (check_open_edge.jl), which found
# no chi2 crossing above either column's upper edge. A 1e-4 re-run can
# therefore only move these edges outward, never inward.
#
# THREE DEPARTURES FROM THE BRIEF, ALL FORCED BY THE RUN.
# 1. The brief asks for lambda2 = 1. That sweep exists on the test-one columns
#    only. On the test-two columns lambda-hat = 100 and the lowest rung swept
#    is 10, so test two is drawn at its profile lambda2 = 10 -- the same
#    lambda2 its bar already used.
# 2. The brief asks for the fine-grid nodes. The band on the strip comes from
#    the SINGLE-PASS 41-node grid for every column (stage4post was never
#    re-run against the two-pass output), so the curve is taken from the same
#    single-pass sweep. Curve and bar then agree by construction, which is the
#    point of sharing the axis. Asserted below.
# 3. Edges are the first and last admissible node, not bisected: the two-pass
#    run recorded edge_solves = 0, so no bisection fired in any saved run.
using JLD2, Statistics, Printf, PyPlot, Random, DelimitedFiles
include("figstyle.jl"); using .FigStyle
const SYN = normpath(joinpath(@__DIR__, "..", "data"))
const OUT = normpath(joinpath(@__DIR__,"..","outputs"))
const MGL = 6400.0

cellsof(z)=(mid=(z[1:end-1].+z[2:end])./2; (top=vcat(0.0,mid), bot=vcat(mid,z[end]+(z[end]-z[end-1])/2)))
function winw(z,a,b); g=cellsof(z); w=[max(0.0,min(g.bot[i],b)-max(g.top[i],a)) for i in eachindex(z)]
    ix=findall(>(1e-12),w); (ix=ix,w=w[ix]); end

"the stage4 file whose grid produced the band in stage4post (see departure 2)"
s4path(tag) = isdir(joinpath(OUT,"synth_$tag","singlepass_archive")) ?
    joinpath(OUT,"synth_$tag","singlepass_archive","stage4_synth_$tag.jld2") :
    joinpath(OUT,"synth_$tag","stage4_synth_$tag.jld2")

function row(tag, ens)
    P4 = load(joinpath(OUT,"synth_$tag","stage4post_synth_$tag.jld2"))
    G  = load(s4path(tag))
    T  = load(joinpath(SYN,"truth_murray$(replace(tag, "J" => "")).jld2"))
    E  = load(joinpath(SYN,"cond_C2_typical$ens.jld2"))
    w  = (G["window"].top, G["window"].base)
    zt, sbt, dzt = T["z"], T["σ_bulk"], T["dz"]
    it = findall(k -> zt[k] > w[1] && zt[k] <= w[2], eachindex(zt))
    ww = winw(E["zsmr"], w...); Θ = Matrix{Float64}(reduce(hcat, E["ens"])); Random.seed!(1)
    ms = 1.5 .+ 0.4 .* rand(4000); js = rand(1:size(Θ,2), 4000)
    K = median([sum((@view Θ[ww.ix,js[t]]) .^ ms[t] .* ww.w) for t in 1:4000])
    l2 = P4["profile_lambda2"]
    sw = G["sweeps"][string(l2)]
    g  = findall(sw.good)
    # curve and bar must be the same object: the admissible node range has to
    # bracket the band, or the two panels are showing different runs.
    @assert G["Cg"][first(g)] <= P4["band_C"].lo * (1 + 1e-9) "lo node/band mismatch on $tag"
    @assert G["Cg"][last(g)]  >= P4["band_C"].hi * (1 - 1e-9) ||
            P4["band_open"].hi == false "hi node/band mismatch on $tag"
    (lo = P4["band_C"].lo, hi = P4["band_C"].hi,
     lo_open = P4["band_open"].lo, hi_open = P4["band_open"].hi,
     truth = sum(sbt[it] .* dzt), ticks = [X/MGL*K for X in (1000.0,1500.0,3000.0)],
     l2 = l2, ndata = G["ndata"],
     Cn = G["Cg"][g], chi2n = sw.chi2[g])
end

"monotone-ish smooth interpolation in log C; falls back to linear"
function smoothcurve(C, chi, lo, hi; n = 400)
    x = log10.(C); xs = range(log10(lo), log10(hi), length = n)
    try
        P = PyPlot.PyCall.pyimport("scipy.interpolate").PchipInterpolator
        f = P(x, chi, extrapolate = true)
        return 10 .^ collect(xs), collect(f(collect(xs)))
    catch
        itp = [begin
                 j = clamp(searchsortedlast(x, t), 1, length(x)-1)
                 chi[j] + (t - x[j])/(x[j+1]-x[j])*(chi[j+1]-chi[j])
               end for t in xs]
        return 10 .^ collect(xs), itp
    end
end

function main()
    # (col label, data, colour, faded?, joint?). The joint bars appear as soon
    # as their Stage 4 post exists, so this script needs no edit when they land.
    foreach(assert_column, ("3P","3S","2P","2S")); provenance(("3P","3S","2P","2S"))
    rows = Any[("test one",         row("3P","3PT"), C_FRESH,  false, false),
               ("test one, saline", row("3S","3ST"), C_FRESH,  true,  false),
               ("test two",         row("2P","2PT"), C_SALINE, false, false),
               ("test two, saline", row("2S","2ST"), C_SALINE, true,  false)]
    for (t, e, lab, ctrl) in (("2PJ","2PJT","test two + sounding",false),
                              ("2SJ","2SJT","test two saline + sounding",true))
        isfile(joinpath(OUT,"synth_$t","stage4post_synth_$t.jld2")) || continue
        push!(rows, (lab, row(t, e), C_SALINE, ctrl, true))
    end
    N      = rows[1][2].ndata
    target = N + 3.84
    XLO, XHI = 0.01, 20.0
    YLO, YHI = 20.0, 56.0

    fig = figure(figsize = (MM170, 4.3))
    gs  = fig.add_gridspec(2, 1, height_ratios = [1.25, 1.0], hspace = 0.12)
    ax  = fig.add_subplot(gs[1])          # chi2 (top)
    bx  = fig.add_subplot(gs[2], sharex = ax)   # bar strip (bottom)

    # ---------------- (a) chi2 of the forced fit ----------------
    ax.axhspan(target, YHI, color = "0.5", alpha = 0.13, lw = 0, zorder = 0)
    ax.axhline(target, color = "0.2", lw = 1.0, zorder = 3)
    ax.axhline(N, color = "0.45", lw = 0.9, ls = ":", zorder = 3)
    ax.text(XLO*1.15, target + 1.0, "excluded by the data", fontsize = 7.0,
            color = "0.25", va = "bottom", ha = "left")
    ax.text(XHI/1.1, N - 0.8, @sprintf("N = %d", N), fontsize = 7.0,
            color = "0.45", va = "top", ha = "right")

    for (name, d, c, ctrl, joint) in rows
        (ctrl || joint) && continue          # controls and the joint bars are not drawn in (a)
        xs, ys = smoothcurve(d.Cn, d.chi2n, d.lo, d.hi)
        ax.plot(xs, ys, color = c, lw = 1.6, zorder = 5, solid_capstyle = "round")
        # truth marker, on the curve
        yt = ys[argmin(abs.(log10.(xs) .- log10(d.truth)))]
        ax.plot(d.truth, yt, "o", ms = 6, mfc = "white", mec = c, mew = 1.5, zorder = 7)
        # direct label, no legend box
        k   = name == "test one" ? 2 : length(xs) - 60
        lab = @sprintf("%s  (λ² = %g)", name, d.l2)
        ax.annotate(lab, (xs[k], ys[k]), textcoords = "offset points",
                    xytext = name == "test one" ? (3, 7) : (-4, 9),
                    ha = name == "test one" ? "left" : "right",
                    fontsize = 7.5, color = c, fontweight = "bold", zorder = 8)
    end
    ax.set_ylim(YLO, YHI); ax.set_ylabel("χ² of the forced fit")
    ax.tick_params(labelbottom = false)
    for s in ("top","right"); ax.spines[s].set_visible(false); end
    panel(ax, "a")

    # ---------------- (b) the bar strip ----------------
    n = length(rows)
    RECT = PyPlot.PyCall.pyimport("matplotlib.patches").Rectangle
    for (i,(name, d, c, ctrl, joint)) in enumerate(rows)
        y = n - i
        a = ctrl ? 0.45 : 1.0
        if joint
            # outlined and hatched, so the joint case reads as a different
            # object at a glance without changing colour
            bx.add_patch(RECT((d.lo, y-0.17), d.hi-d.lo, 0.34, facecolor = "none",
                              edgecolor = c, lw = 1.1, hatch = "///", alpha = a, zorder = 3))
        else
            bx.plot([d.lo, d.hi], [y, y], color = c, lw = 7, alpha = a, solid_capstyle = "butt")
        end
        d.lo_open && bx.plot(d.lo, y, marker = "<", ms = 7, color = c, alpha = a, clip_on = false)
        d.hi_open && bx.plot(d.hi, y, marker = ">", ms = 7, color = c, alpha = a, clip_on = false)
        bx.plot(d.truth, y, "o", ms = 7, mfc = "white", mec = c, mew = 1.6, zorder = 5)
        for (t, lab) in zip(d.ticks, ("1,000","1,500","3,000 mg/L"))
            bx.plot([t, t], [y-0.3, y+0.3], color = "0.2", lw = 1.0, zorder = 4)
            ctrl || bx.text(t, y+0.36, lab, ha = "center", va = "bottom", fontsize = 6.2, color = "0.2")
        end
        bx.text(0.0105, y, name, ha = "left", va = "center", fontsize = 7.5)
    end
    bx.set_xscale("log"); bx.set_xlim(XLO, XHI)
    bx.set_ylim(-0.6, n - 0.05); bx.set_yticks([])
    bx.set_xlabel("window conductance the data allow (S)")
    L2D = PyPlot.PyCall.pyimport("matplotlib.lines").Line2D
    h1 = L2D([0],[0], ls = "none", marker = "o", ms = 7, mfc = "white", mec = "0.3",
             mew = 1.6, label = "truth")
    h2 = L2D([0],[0], ls = "none", marker = ">", ms = 7, color = "0.3",
             label = "limit of the model, not the data")
    hs = [h1, h2]
    if any(r -> r[5], rows)
        PATCH = PyPlot.PyCall.pyimport("matplotlib.patches").Patch
        push!(hs, PATCH(facecolor = "none", edgecolor = "0.3", hatch = "///",
                        label = "AEM + galvanic sounding"))
    end
    bx.legend(handles = hs, loc = "upper right", frameon = false, fontsize = 6.8,
              handletextpad = 0.4, labelspacing = 0.35, borderaxespad = 0.0)
    for s in ("top","right","left"); bx.spines[s].set_visible(false); end
    panel(bx, "b")

    fig.subplots_adjust(left = 0.075, right = 0.985, top = 0.98, bottom = 0.12)
    FigStyle.finish(fig, "fig07_gate_sigmin1e-3")

    # ---------------- export: per-node chi2 ----------------
    tagof = Dict("test one"=>"3P", "test one, saline"=>"3S",
                 "test two"=>"2P", "test two, saline"=>"2S")
    open(joinpath(@__DIR__, "numbers_chi2.csv"), "w") do f
        println(f, "column,tag,C_node_S,chi2,lambda2,in_band,figure,source")
        nr = 0
        for (name, d, _, _) in rows, i in eachindex(d.Cn)
            inb = d.Cn[i] >= d.lo*(1-1e-9) && d.Cn[i] <= d.hi*(1+1e-9)
            @printf(f, "\"%s\",\"%s\",%.6g,%.6g,%g,%s,\"7a\",\"stage4_synth_%s.jld2\"\n",
                    name, tagof[name], d.Cn[i], d.chi2n[i], d.l2, inb, tagof[name])
            nr += 1
        end
        @printf("  wrote %-34s %d rows\n", "numbers_chi2.csv", nr)
    end
end
main()
