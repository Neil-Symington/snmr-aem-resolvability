# fig08b_clay_salinity_tradeoff.jl -- the clay / salinity trade-off.
#
# NOTE ON NUMBERING: `fig08_package_vs_window` already holds Figure 8 from the
# earlier brief. This is written as 08b pending a renumber; nothing downstream
# assumes either number.
#
# One panel per column. x: the clay term sigma_clay*C in S/m, log. y: TDS in
# mg/L, log. The density of ACCEPTED Stage 5 draws -- those with sigma_w > 0 --
# over white. Truth as a filled marker, the 1,000 mg/L line dashed.
#
# TWO TRADE-OFF CURVES. For a fixed bulk conductivity and a fixed theta^m, the
# rock physics fixes the trade-off exactly:
#     TDS(clay) = (sigma_bulk - clay) / theta^m  x  mg_per_sm
# drawn once at the AEM's MEDIAN window conductivity -- the line the draws
# actually follow -- and once at the TRUE window conductivity, the line they
# should follow. The gap between them is the conductance error, and on test two
# it is the 2.66x: the draws sit on the upper curve and the truth on the lower.
# Each curve meets TDS = 0 at clay = sigma_bulk, which is the clay floor.
using JLD2, Statistics, Printf, PyPlot
include("figstyle.jl"); using .FigStyle
const ROOT = normpath(joinpath(@__DIR__, ".."))
const SYN = normpath(joinpath(@__DIR__, "..", "data"))
const OUT = normpath(joinpath(@__DIR__,"..","outputs"))
const MGL = 6400.0

function main()
    # Panel (c) appears as soon as the joint Stage 5 exists, so this script
    # builds either way and nothing has to be edited when the run lands.
    cols = Any[("3P","test one — unconfined"), ("2P","test two — confined")]
    foreach(c -> assert_column(c[1]), cols); provenance(first.(cols))
    # Stage 5 is what this figure reads, and a missing one used to surface as a
    # raw ArgumentError naming a path with no pointer to what writes it. A
    # reviewer hit that on a tree that had 3P Stage 0 and 3 but not 4 and 5.
    for (tag, _) in cols
        f = joinpath(ROOT, "outputs", "synth_$tag", "stage5_synth_$tag.jld2")
        isfile(f) || error("""
    no Stage 5 output for column $tag: $(relpath(f, ROOT))

      This figure reads Stage 5. Produce it with
        julia --project=. run_site.jl configs/synth_$tag.toml
      Stage 4 must bracket first -- for 3P the smoke gate grid (4 coarse,
      3 fine) does NOT, and an 8/4 grid does. See the README's Figures section.
    """)
    end
    isfile(joinpath(OUT,"synth_2PJ","stage5_synth_2PJ.jld2")) &&
        push!(cols, ("2PJ","test two + sounding"))
    fig, axs = FigStyle.newfig(width = MM170, h = 3.4, nrows = 1, ncols = length(cols))
    xl = (1e-4, 1.0); yl = (30.0, 1e5)
    for (i,(tag,label)) in enumerate(cols)
        ax = axs[i]
        S5 = load(joinpath(OUT,"synth_$tag","stage5_synth_$tag.jld2"))
        G  = load(joinpath(OUT,"synth_$tag","stage4_synth_$tag.jld2"))
        T  = load(joinpath(SYN,"truth_murray$(replace(tag, "J" => "")).jld2"))
        w  = (G["window"].top, G["window"].base)

        pos = findall(>(0), S5["sigma_w"])
        cl  = S5["clay_product"][pos]; td = S5["TDS"][pos]

        # density of accepted draws, over white
        xe = 10 .^ range(log10(xl[1]), log10(xl[2]), length = 80)
        ye = 10 .^ range(log10(yl[1]), log10(yl[2]), length = 80)
        H = zeros(length(ye)-1, length(xe)-1)
        for (x,y) in zip(cl, td)
            bx = searchsortedlast(xe, x); by = searchsortedlast(ye, y)
            (1 <= bx <= size(H,2) && 1 <= by <= size(H,1)) && (H[by,bx] += 1)
        end
        H[H .== 0] .= NaN
        cm = PyPlot.matplotlib.pyplot.get_cmap("YlGnBu").copy(); cm.set_bad("white")
        pcm = ax.pcolormesh(xe, ye, H, cmap = cm, shading = "flat",
                            norm = PyPlot.matplotlib.colors.LogNorm(vmin = 1, vmax = maximum(filter(isfinite,H))),
                            rasterized = true)

        # the two trade-off curves, at the median theta^m of the accepted draws
        thm = median(S5["theta_draw"][pos] .^ S5["m_draw"][pos])
        sb_aem = median(S5["sigma_bulk"])
        zt, sbt, dzt = T["z"], T["σ_bulk"], T["dz"]
        it = findall(k -> zt[k] > w[1] && zt[k] <= w[2], eachindex(zt))
        sb_true = sum(sbt[it] .* dzt) / (w[2] - w[1])
        cgrid = 10 .^ range(log10(xl[1]), log10(xl[2]), length = 400)
        for (sb, ls, nm) in ((sb_aem, "-", "trade-off at the AEM median σ_bulk"),
                             (sb_true, "--", "trade-off at the true σ_bulk"))
            y = (sb .- cgrid) ./ thm .* MGL
            k = findall(v -> v > yl[1], y)
            ax.plot(cgrid[k], y[k], color = C_TRUTH, lw = 0.9, ls = ls, label = nm)
        end

        # truth
        swtrue = mean(T["σw"][it]) * MGL
        claytrue = sb_true - mean(T["σw"][it]) * mean(T["θ_total"][it])^FigStyle.archie_m(T)
        ax.plot([claytrue], [swtrue], "o", ms = 6, mfc = C_SALINE, mec = "black",
                mew = 0.8, zorder = 8, label = "truth")

        ax.axhline(1000, color = C_MGR, ls = "--", lw = 0.9, label = "1,000 mg/L")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(xl); ax.set_ylim(yl)
        ax.set_xlabel("clay term (S/m)")
        i == 1 ? ax.set_ylabel("TDS (mg/L)") : ax.set_yticklabels([])
        ax.set_title(label, fontsize = 8, pad = 3)
        FigStyle.panel(ax, "abc"[i])
        # the numbers live in the corner box, not the legend, so one shared
        # legend can serve both panels and nothing sits under the panel label
        ax.text(0.975, 0.035,
                @sprintf("acceptance %.0f%%\nAEM σ_bulk %.4f S/m\ntrue σ_bulk %.4f S/m",
                         100*S5["acceptance"], sb_aem, sb_true),
                transform = ax.transAxes, ha = "right", va = "bottom", fontsize = 6.2,
                linespacing = 1.35,
                bbox = Dict("facecolor"=>"white","edgecolor"=>"none","alpha"=>0.85,"pad"=>2.0))
        if i == length(cols)
            cb = fig.colorbar(pcm, ax = ax, pad = 0.02, fraction = 0.045)
            cb.set_label("accepted draws", fontsize = 6.5); cb.ax.tick_params(labelsize = 6)
        end
    end
    lg = axs[1].legend(loc = "upper right", fontsize = 6.0)
    lg.set_zorder(20)
    fig.subplots_adjust(wspace = 0.08, left = 0.085, right = 0.93, top = 0.925, bottom = 0.145)
    FigStyle.finish(fig, "fig06_clay_salinity_tradeoff")
end
main()
