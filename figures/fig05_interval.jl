# fig05_interval.jl -- Figure 5. What Stage 3 chooses, row per column.
#
# (a, d) the per-draw PREFERRED window: for each SNMR draw, the window inside
#        the passing zone holding the most water. Stage 3 already stores these
#        as `tops` and `bases` -- its own argmax-W per draw -- and on both
#        columns 100% of them lie wholly inside a primary (p_pass >= 0.5) zone,
#        so no extra restriction is imposed. Counts on a log axis.
# (b, e) the W posterior over the representative window Stages 4 and 5 used.
#
# Panel letters run a, b, c, d. The true aquifer is read from the truth file's
# unit table rather than hardcoded -- test two's Calivil aquifer is 18-34 m,
# and an earlier version drew 21-31 m, which is that aquifer inset by 3 m and
# not its extent.
#
# EVERY DRAW'S WINDOW IS EXACTLY 16 m LONG on both columns -- L pins at L_max
# in 100% of draws, which is spec-determined and not inferred -- so (a) and (d)
# show one distribution shifted by 16 m, not two independent ones.
using JLD2, Statistics, Printf, PyPlot
include("figstyle.jl"); using .FigStyle
const SYN = normpath(joinpath(@__DIR__, "..", "data"))
const OUT = normpath(joinpath(@__DIR__,"..","outputs"))
const WMIN = 2.0

"the first unit whose label names an aquifer, as (top, base) in metres"
function aquifer(T)
    for u in T["units"]
        k = u isa Dict ? u : Dict(string(f)=>getfield(u,f) for f in propertynames(u))
        occursin("AQUIFER", uppercase(string(k["3"]))) && return (Float64(k["1"]), Float64(k["2"]))
    end
    error("no aquifer unit in the truth")
end
cellsof(z)=(mid=(z[1:end-1].+z[2:end])./2; (top=vcat(0.0,mid), bot=vcat(mid,z[end]+(z[end]-z[end-1])/2)))
function winw(z,a,b); g=cellsof(z); w=[max(0.0,min(g.bot[i],b)-max(g.top[i],a)) for i in eachindex(z)]
    ix=findall(>(1e-12),w); (ix=ix,w=w[ix]); end

function main()
    cols = (("3P","3PT", C_FRESH,  "#9dbfe0", (7.5,23.5)),
            ("2P","2PT", C_SALINE, "#f3c69b", (17.5,33.5)))
    foreach(c -> assert_column(c[1]), cols); provenance(first.(cols))
    fig, axs = FigStyle.newfig(width = MM170, h = 4.4, nrows = 2, ncols = 2)
    letters = (("a","b"), ("c","d"))
    exp_win = Any[]; exp_W = Any[]
    for (i,(tag, ens, c1, c2, rep)) in enumerate(cols)
        aq = aquifer(load(joinpath(SYN,"truth_murray$tag.jld2")))
        S3 = load(joinpath(OUT,"synth_$tag","stage3_synth_$tag.jld2"))
        tp, bs = S3["tops"], S3["bases"]
        append!(exp_win, [(tag, tp[k], bs[k]) for k in eachindex(tp)])

        # ---- (a, d) preferred window top and base ----
        ax = axs[i, 1]
        FigStyle.aquifer_band(ax, aq[1], aq[2]; horiz = true)
        # the aquifer's own boundaries, so the band has edges a reader can place
        for zb in aq
            ax.axvline(zb, color = "0.45", ls = "--", lw = 0.8, zorder = 2)
        end
        bins = collect(0.0:0.5:60.0)
        ax.hist(tp, bins = bins, color = c1, lw = 0, label = "top")
        ax.hist(bs, bins = bins, color = c2, lw = 0, label = "base")
        ax.set_yscale("log"); ax.set_xlim(0, 60)
        ax.set_xlabel("depth (m)"); ax.set_ylabel("draws")
        ax.legend(loc = "upper right", fontsize = 6.4)
        FigStyle.panel(ax, letters[i][1])
        for s in ("top","right"); ax.spines[s].set_visible(false); end

        # ---- (b, e) W over the representative window ----
        bx = axs[i, 2]
        E = load(joinpath(SYN,"cond_C2_typical$ens.jld2"))
        Θ = Matrix{Float64}(reduce(hcat, E["ens"])); ww = winw(E["zsmr"], rep...)
        Wr = [sum(Θ[ww.ix,j] .* ww.w) for j in 1:size(Θ,2)]
        append!(exp_W, [(tag, w) for w in Wr])
        bx.hist(Wr, bins = 60, color = c1, alpha = 0.85, lw = 0)
        bx.axvline(WMIN, color = "0.25", ls = "--", lw = 1.0)
        bx.set_yticks([])
        bx.set_xlabel(@sprintf("W over %.1f–%.1f m (m)", rep...))
        FigStyle.panel(bx, letters[i][2])
        for s in ("top","right","left"); bx.spines[s].set_visible(false); end
        @printf("%s  top p5-p95 %.2f-%.2f  base p5-p95 %.2f-%.2f | W p50 %.3f  P(W>=2) %.4f\n",
                tag, quantile(tp,.05), quantile(tp,.95), quantile(bs,.05), quantile(bs,.95),
                median(Wr), count(>=(WMIN),Wr)/length(Wr))
    end
    fig.subplots_adjust(left = 0.075, right = 0.985, top = 0.965, bottom = 0.105,
                        wspace = 0.16, hspace = 0.34)
    FigStyle.finish(fig, "fig05_interval")

    open(joinpath(@__DIR__,"numbers_fig05_window.csv"),"w") do f
        println(f, "column,draw,preferred_top_m,preferred_base_m,figure,source")
        for (k,(t,a,b)) in enumerate(exp_win)
            @printf(f, "\"%s\",%d,%.6g,%.6g,\"5a/5d\",\"stage3_synth_%s.jld2\"\n", t, k, a, b, t)
        end
    end
    open(joinpath(@__DIR__,"numbers_fig05_W.csv"),"w") do f
        println(f, "column,W_m,figure,source")
        for (t,w) in exp_W
            @printf(f, "\"%s\",%.6g,\"5b/5e\",\"cond_C2_typical%sT.jld2\"\n", t, w, t)
        end
    end
    @printf("  wrote numbers_fig05_window.csv (%d rows), numbers_fig05_W.csv (%d rows);\n",
            length(exp_win), length(exp_W))
    println("  numbers_fig05_ppass.csv kept from the previous build")
end
main()
