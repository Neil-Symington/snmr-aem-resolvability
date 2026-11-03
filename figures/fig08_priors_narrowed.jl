# fig08_priors_narrowed.jl -- Figure 8. Narrowing the Stage 0 priors with
# borehole logs, test two. One panel, six rows.
#
# Rows 1-4 are the prior states at the bore itself; rows 5-6, set apart, are
# the transfer test -- the clay term fixed at its TRUE value for the window and
# carried to a survey-only neighbour. The saline neighbour has its own truth,
# 1,736 mg/L, which is why the open ring moves on that row.
#
# All bars grey: nothing on this column passes. The orange PASS rule is kept in
# the code for the field version, but it has no legend entry here because it
# never fires -- a legend key for a colour that does not appear is clutter.
using JLD2, Statistics, Printf, PyPlot
include("figstyle.jl"); using .FigStyle
const S6 = normpath(joinpath(@__DIR__,"..","outputs","section6","section6_2P.jld2"))

thou(x) = begin
    s = @sprintf("%.0f", x); o = ""
    for (i, ch) in enumerate(reverse(collect(s)))
        o = string(ch) * ((i > 1 && (i-1) % 3 == 0) ? "," : "") * o
    end
    o
end

function main()
    D = load(S6)
    order = D["order"]; names = D["names"]
    truth = D["truth_TDS"]; Xlim = D["Xlim"]
    XLO, XHI = 100.0, 1e5

    # (label, p05, p50, p95, Pfresh, truth, passes?)
    rows = Any[]
    for k in order
        s = D["state_$k"]
        push!(rows, (names[k], s.p05, s.p50, s.p95, s.Pfresh, truth, s.verdict == "PASS"))
    end
    for (k, lab) in (("2P_fixed","clay term fixed, fresh neighbour"),
                     ("2S_fixed","clay term fixed, saline neighbour"))
        x = D["xfer_$k"]
        push!(rows, (lab, x.p05, x.p50, x.p95, x.Pfresh, x.truth, false))
    end
    # four rows, a gap, then the two transfer rows
    ys = [5.0, 4.0, 3.0, 2.0, 0.7, -0.3]

    fig, ax = FigStyle.newfig(width = MM170, h = 2.9)
    ax.axvline(Xlim, color = "0.25", lw = 1.0, ls = "--", zorder = 2)
    for (i, (lab, p05, p50, p95, Pf, tr, pass)) in enumerate(rows)
        y = ys[i]
        c  = pass ? C_SALINE : "0.62"
        cm = pass ? "#a8551a" : "0.28"
        # truth ring first and unfilled, so it never hides a p50 sitting on it
        ax.plot(tr, y, "o", ms = 7, mfc = "none", mec = "black", mew = 1.2, zorder = 4)
        ax.plot([p05, p95], [y, y], color = c, lw = 7, solid_capstyle = "butt", zorder = 3)
        ax.plot(p50, y, "o", ms = 5, mfc = cm, mec = "white", mew = 0.9, zorder = 6)
        ax.text(p95*1.25, y, @sprintf("P(fresh) %.2f", Pf), ha = "left", va = "center",
                fontsize = 6.4, color = "0.3")
    end
    ax.text(Xlim*1.07, ys[1]+0.45, @sprintf("%s mg/L", thou(Xlim)),
            fontsize = 6.8, color = "0.25", ha = "left", va = "bottom")
    ax.set_xscale("log"); ax.set_xlim(XLO, XHI)
    ax.set_ylim(minimum(ys)-0.65, maximum(ys)+0.75)
    ax.set_yticks(ys); ax.set_yticklabels([r[1] for r in rows], fontsize = 7.5)
    ax.tick_params(axis = "y", length = 0)
    ax.set_xlabel("pore-water TDS (mg/L)")
    L2D = PyPlot.PyCall.pyimport("matplotlib.lines").Line2D
    ax.legend(handles = [L2D([0],[0], ls = "none", marker = "o", ms = 7, mfc = "none",
                             mec = "black", mew = 1.2, label = "truth")],
              loc = "lower right", frameon = false, fontsize = 6.8,
              handletextpad = 0.5, borderaxespad = 0.1)
    for s in ("top","right","left"); ax.spines[s].set_visible(false); end
    fig.subplots_adjust(left = 0.215, right = 0.99, top = 0.97, bottom = 0.17)
    FigStyle.finish(fig, "fig08_priors_narrowed")
end
main()
