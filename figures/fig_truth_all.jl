# fig_truth_all.jl -- the true conductivity and water content of every column.
#
# Three grounds across, two quantities down. Fresh and brackish differ only in
# pore water, so they share one water-content panel: the lower row carries one
# pair of curves per ground, not two.
#
# Light theme, as Figure 2. Depth 0-70 m. Words, not symbols.
using JLD2, Printf, PyPlot
include("figstyle.jl")

const DATA = joinpath(@__DIR__, "..", "data")
const ZMAX = 70.0
const BG, PANEL, INK, MUTE, RULE = "#ffffff", "#f4f6f8", "#16202a", "#5a6b78", "#c9d3da"
const TRUTH = "#16202a"
# the brackish colours of the gate figures, darkened for white paper
const BRACK = Dict("3S" => "#6d28d9", "2S" => "#047857", "2Sv2" => "#047857")

ld(t) = load(joinpath(DATA, "truth_murray$t.jld2"))

const GROUNDS = (
  (fresh = "3P",   brack = "3S",   win = (7.5, 23.5),  label = "Test one (unconfined)"),
  (fresh = "2P",   brack = "2S",   win = (17.5, 33.5), label = "Test two (confined), v1"),
  (fresh = "2Pv2", brack = "2Sv2", win = (17.5, 33.5), label = "Test two (confined), v2"),
)

function style!()
    rc("figure", facecolor = BG); rc("savefig", facecolor = BG, dpi = 300)
    rc("axes", facecolor = PANEL, edgecolor = RULE, labelcolor = INK, linewidth = 0.7)
    rc("xtick", color = MUTE, labelsize = 7.0); rc("ytick", color = MUTE, labelsize = 7.0)
    rc("text", color = INK); rc("font", family = "sans-serif", size = 7.5)
end

function main()
    style!()
    fig, axs = subplots(2, 3, figsize = (180/25.4, 130/25.4))
    letters = "abcdef"

    for (c, g) in enumerate(GROUNDS)
        TF, TS = ld(g.fresh), ld(g.brack)
        z = collect(Float64, TF["z"]); zi = findall(<=(ZMAX), z)
        bc = BRACK[g.brack]

        # ---- row 1: bulk conductivity ----
        ax = axs[1, c]
        ax.plot(collect(Float64, TS["σ_bulk"])[zi], z[zi], color = bc, lw = 1.3,
                label = "brackish water")
        ax.plot(collect(Float64, TF["σ_bulk"])[zi], z[zi], color = TRUTH, lw = 1.1,
                label = "fresh water")
        ax.set_xscale("log"); ax.set_xlim(1e-3, 3)
        ax.set_xticks([1e-3, 1e-2, 1e-1, 1.0])
        ax.set_xticklabels(["0.001", "0.01", "0.1", "1"])
        ax.set_xlabel("bulk conductivity (S/m)")

        # ---- row 2: water content (identical for the pair) ----
        ax = axs[2, c]
        @assert isapprox(sum(TF["θ_total"]), sum(TS["θ_total"]); rtol = 1e-12) "\
            $(g.fresh) and $(g.brack) do not share a water column; the lower row \
            may not be drawn once for the pair"
        ax.plot(collect(Float64, TF["θ_total"])[zi], z[zi], color = TRUTH, lw = 1.1,
                label = "total porosity")
        ax.plot(collect(Float64, TF["θ_SMR"])[zi], z[zi], color = TRUTH, lw = 1.1,
                ls = (0, (4, 2)), label = "SNMR water content")
        ax.set_xlim(0, 0.45); ax.set_xlabel("water content (v/v)")

        for r in 1:2
            a = axs[r, c]
            a.set_ylim(ZMAX, 0)
            for w in g.win
                a.axhline(w, color = "#3d4e5c", lw = 0.8, ls = (0, (4, 3)), zorder = 5)
            end
            c == 1 ? a.set_ylabel("depth (m)") : a.set_yticklabels([])
            for k in ("top","bottom","left","right"); a.spines[k].set_color(RULE); end
            a.text(0.0, 1.012, string("(", letters[(r-1)*3 + c], ")"),
                   transform = a.transAxes, fontsize = 7.5, fontweight = "bold",
                   color = INK, va = "bottom")
        end
        axs[1, c].text(0.0, 1.075, g.label, transform = axs[1, c].transAxes,
                       fontsize = 8.0, fontweight = "bold", color = INK, va = "bottom")
    end

    L = matplotlib.lines.Line2D
    hs = [L([], [], color = TRUTH, lw = 1.1),
          L([], [], color = BRACK["3S"], lw = 1.3),
          L([], [], color = BRACK["2S"], lw = 1.3),
          L([], [], color = TRUTH, lw = 1.1, ls = (0, (4, 2)))]
    ls = ["fresh water, and total porosity", "brackish water, test one",
          "brackish water, test two", "SNMR water content"]
    leg = fig.legend(hs, ls, loc = "lower center", ncol = 4, frameon = false,
                     bbox_to_anchor = (0.5, -0.02), fontsize = 6.8, handlelength = 1.8)
    for t in leg.get_texts(); t.set_color(INK); end

    fig.subplots_adjust(wspace = 0.12, hspace = 0.42, left = 0.085, right = 0.985,
                        top = 0.93, bottom = 0.14)
    p = joinpath(@__DIR__, "fig_truth_all.png")
    fig.savefig(p, dpi = 300, facecolor = BG, bbox_inches = "tight", pad_inches = 0.03)
    fig.savefig(replace(p, ".png" => ".pdf"), facecolor = BG, bbox_inches = "tight", pad_inches = 0.03)
    println("  wrote ", p)
end
main()
