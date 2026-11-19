# fig02_columns.jl -- Figure 2. The two test grounds, and what changes between
# fresh and brackish pore water.
#
# Dark theme, the depth figures' palette and panel styling. Five panels a row:
# the stratigraphic column in the Synthetic Section Builder's symbology, then
# water, clay, pore-water salinity and bulk conductivity. One shared key below.
#
# WORDS, NOT SYMBOLS, throughout: "water content (v/v)", "bulk conductivity
# (S/m)". Fresh is white, as the truth is in Figures 3 and 4; the brackish
# partner takes the colour its column has in the gate figures.
using JLD2, Statistics, Printf, PyPlot, TOML
include("lithology.jl"); using .Lithology

const SYN  = joinpath(@__DIR__, "..", "data")
const CFG  = joinpath(@__DIR__, "..", "synthetics", "config")
const OUT  = @__DIR__
const MGL  = 6400.0
const ZMAX = 70.0

# the depth figures' palette
const BG, PANEL, INK, MUTE, RULE = "#0b0d10", "#162023", "#e3ebea", "#8098a0", "#2a393d"
const WHITE = "#ffffff"
const BRACK = Dict("3S" => "#a78bfa", "2S" => "#34d399")   # as in the gate figures

ld(t) = load(joinpath(SYN, "truth_murray$t.jld2"))

"""
Role, not rock type. The .toml names carry formation names and must not appear;
the key gives the lithologies. Roles are assigned from the unit's position and
properties, so the two columns are labelled consistently.
"""
function roles(us, zmax)
    tops = Float64[Float64(u["top"]) for u in us]
    out = Tuple{Float64,Float64,String,String}[]      # top, base, lith, role
    naq = 0
    for (i, u) in enumerate(us)
        base = i < length(us) ? tops[i+1] : zmax
        lith = Lithology.lithof(u)
        nm = uppercase(string(get(u, "name", "")))
        role = if occursin("UNSATURATED", nm); "unsaturated"
        elseif occursin("AQUITARD", nm) || occursin("CONFINING", nm); "clay aquitard"
        elseif occursin("CAP", nm); "clay cap"
        elseif occursin("AQUIFER", nm); naq += 1; naq == 1 ? "aquifer" : "lower sand"
        elseif lith in ("clay", "sandy clay"); "clay cap"
        else; "lower sand"
        end
        push!(out, (tops[i], base, lith, role))
    end
    out
end

watertable(us) = (i = findfirst(u -> occursin("UNSATURATED", uppercase(string(get(u,"name","")))), us);
                  i === nothing ? NaN : Float64(us[min(i+1, length(us))]["top"]))

function style!()
    rc("figure", facecolor = BG); rc("savefig", facecolor = BG, dpi = 300)
    rc("axes", facecolor = PANEL, edgecolor = RULE, labelcolor = INK, linewidth = 0.7)
    rc("xtick", color = MUTE, labelsize = 7.0); rc("ytick", color = MUTE, labelsize = 7.0)
    rc("text", color = INK)
    rc("font", family = "sans-serif", size = 7.5)
end

function main()
    style!()
    cols = (("3P", "3S", "Test one (unconfined)", (7.5, 23.5), "paper_3P.toml"),
            ("2P", "2S", "Test two (confined)",   (17.5, 33.5), "confined_aquifer.toml"))
    fig, axs = subplots(2, 5, figsize = (180/25.4, 170/25.4),
                        gridspec_kw = Dict("width_ratios" => [0.80, 1, 1, 1, 1]))
    letters = "abcdefghij"; li = 1; allliths = String[]

    for (r, (ft, st, label, win, cfgname)) in enumerate(cols)
        TF, TS = ld(ft), ld(st)
        z = TF["z"]; zi = findall(<=(ZMAX), z)
        us = TOML.parsefile(joinpath(CFG, cfgname))["unit"]
        blocks = roles(us, ZMAX)
        bc = BRACK[st]

        # ---- 0. stratigraphic column ----------------------------------
        ax = axs[r, 1]
        Lithology.drawcolumn(ax, [(b[1], b[2], b[3]) for b in blocks], ZMAX;
                             dark = true, px_per_m = 9.0)
        append!(allliths, [b[3] for b in blocks if b[1] < ZMAX])
        ax.set_xlim(0, 1.45)          # the right third carries the window bracket
        seen = String[]               # a role is named once; two units share "lower sand"
        for b in blocks
            b[1] >= ZMAX && continue
            (min(b[2], ZMAX) - b[1]) < 0.055 * ZMAX && continue
            b[4] in seen && continue
            push!(seen, b[4])
            zc = (b[1] + min(b[2], ZMAX))/2
            for w in win       # keep the name off the window lines
                abs(zc - w) < 2.0 && (zc = w + (zc >= w ? 2.5 : -2.5))
            end
            ax.text(0.5, zc, b[4], ha = "center", va = "center",
                    fontsize = 5.6, color = "#0b0f14", zorder = 6,
                    bbox = Dict("facecolor" => "#e9eff5", "alpha" => 0.88,
                                "edgecolor" => "none", "pad" => 1.0))
        end
        wt = watertable(us)
        if isfinite(wt)
            ax.plot([0.12], [wt], marker = "v", ms = 6, color = WHITE,
                    markeredgecolor = "#0b0f14", markeredgewidth = 0.5, zorder = 7)
        end
        # the gate window, as a bracket on the right edge
        ax.plot([1.06, 1.14, 1.14, 1.06], [win[1], win[1], win[2], win[2]],
                color = INK, lw = 1.0, zorder = 7)
        ax.text(1.20, (win[1]+win[2])/2, "gate window", rotation = 270, ha = "left",
                va = "center", fontsize = 6.0, color = INK)
        ax.set_ylabel("depth (m)")

        # ---- 1. water -------------------------------------------------
        ax = axs[r, 2]
        ax.plot(TF["θ_total"][zi], z[zi], color = WHITE, lw = 1.1)
        ax.plot(TF["θ_SMR"][zi], z[zi], color = WHITE, lw = 1.1, ls = (0, (4, 2)))
        ax.set_xlabel("water content (v/v)"); ax.set_xlim(0, 0.45)

        # ---- 2. clay --------------------------------------------------
        ax = axs[r, 3]
        ax.plot(TF["clay"][zi], z[zi], color = WHITE, lw = 1.1)
        ax.set_xlabel("clay fraction (v/v)"); ax.set_xlim(0, 1.0)

        # ---- 3. pore-water salinity -----------------------------------
        ax = axs[r, 4]
        ax.plot(TS["σw"][zi] .* MGL, z[zi], color = bc, lw = 1.3)
        ax.plot(TF["σw"][zi] .* MGL, z[zi], color = WHITE, lw = 1.1)
        for (v, nm) in ((1500.0, "1,500"), (3000.0, "3,000"))
            ax.axvline(v, color = RULE, lw = 0.9, ls = (0, (3, 3)), zorder = 1)
            r == 1 && ax.text(v, ZMAX - 1.5, " $nm", rotation = 90, fontsize = 5.8,
                              color = MUTE, va = "bottom", ha = "left")
        end
        ax.set_xlabel("pore-water salinity (mg/L)"); ax.set_xscale("log")
        ax.set_xlim(60, 2e4)

        # ---- 4. bulk conductivity -------------------------------------
        ax = axs[r, 5]
        m = Float64(TF["m_archie"])
        clayterm = max.(TF["σ_bulk"][zi] .- TF["σw"][zi] .* TF["θ_total"][zi] .^ m, 1e-6)
        ax.fill_betweenx(z[zi], 1e-3, clayterm, color = "#9aa6ad", alpha = 0.33, lw = 0)
        ax.plot(TS["σ_bulk"][zi], z[zi], color = bc, lw = 1.3)
        ax.plot(TF["σ_bulk"][zi], z[zi], color = WHITE, lw = 1.1)
        ax.set_xlabel("bulk conductivity (S/m)"); ax.set_xscale("log")
        ax.set_xlim(1e-3, 3)
        ax.set_xticks([1e-3, 1e-2, 1e-1, 1.0])
        ax.set_xticklabels(["0.001", "0.01", "0.1", "1"])

        for c in 1:5
            a = axs[r, c]
            a.set_ylim(ZMAX, 0)
            for w in win
                a.axhline(w, color = "#cfd8dd", lw = 0.8, ls = (0, (4, 3)), zorder = 5)
            end
            c > 1 && a.set_yticklabels([])
            for k in ("top","bottom","left","right"); a.spines[k].set_color(RULE); end
            a.text(0.0, 1.012, string("(", letters[li], ")"), transform = a.transAxes,
                   fontsize = 7.5, fontweight = "bold", color = INK, va = "bottom")
            li += 1
        end
        axs[r, 1].text(0.0, 1.075, label, transform = axs[r, 1].transAxes,
                       fontsize = 8.5, fontweight = "bold", color = INK, va = "bottom")
    end

    # ---- one shared key below --------------------------------------
    L = matplotlib.lines.Line2D; P = matplotlib.patches.Patch
    hs = Any[L([], [], color = WHITE, lw = 1.1),
             L([], [], color = WHITE, lw = 1.1, ls = (0, (4, 2))),
             L([], [], color = BRACK["3S"], lw = 1.3),
             L([], [], color = BRACK["2S"], lw = 1.3),
             P(facecolor = "#9aa6ad", alpha = 0.33, edgecolor = "none")]
    ls = ["total porosity; fresh pore water", "SNMR water content",
          "brackish partner, test one", "brackish partner, test two", "clay term"]
    present = [l for l in Lithology.LITH_ORDER if l in unique(allliths)]
    for l in present
        push!(hs, P(facecolor = Lithology.tintof(l; dark = true),
                    edgecolor = INK, linewidth = 0.4)); push!(ls, l)
    end
    leg = fig.legend(hs, ls, loc = "lower center", ncol = 5, frameon = false,
                     bbox_to_anchor = (0.5, -0.005), fontsize = 6.4,
                     handlelength = 1.6, columnspacing = 1.6)
    for t in leg.get_texts(); t.set_color(INK); end

    fig.subplots_adjust(wspace = 0.14, hspace = 0.30, left = 0.075, right = 0.985,
                        top = 0.935, bottom = 0.14)
    p = joinpath(OUT, "fig02_columns.png")
    fig.savefig(p, dpi = 300, facecolor = BG, bbox_inches = "tight", pad_inches = 0.03)
    fig.savefig(replace(p, ".png" => ".pdf"), facecolor = BG, bbox_inches = "tight", pad_inches = 0.03)
    println("  wrote ", p)
end
main()
