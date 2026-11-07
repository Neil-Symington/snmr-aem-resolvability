# fig02_columns.jl -- Figure 2. The two synthetic columns.
#
#   julia --project=. fig02_columns.jl
#
# Row per column (test one = 3P unconfined, test two = 2P confined). Panels share
# one depth axis 0-60 m: (a) total porosity and SNMR water content; (b) clay
# fraction; (c) pore-water salinity against the 1,000 mg/L line; (d) bulk
# conductivity split into its water and clay terms, stacked. Fresh and saline
# control overlaid in (c) and (d) only, which are the only panels where they
# differ -- the controls share theta and clay with their parents bit-for-bit.
using JLD2, Statistics, Printf, PyPlot
include("figstyle.jl"); using .FigStyle
const SYN = normpath(joinpath(@__DIR__, "..", "data"))
const MGL = 6400.0
ld(t) = load(joinpath(SYN, "truth_murray$t.jld2"))

# aquifer intervals, from each column's own unit table
function aquifers(T)
    out = Tuple{Float64,Float64}[]
    for u in T["units"]
        k = u isa Dict ? u : Dict(string(f) => getfield(u, f) for f in propertynames(u))
        occursin("AQUIFER", uppercase(string(k["3"]))) && push!(out, (Float64(k["1"]), Float64(k["2"])))
    end
    out
end

function main()
    cols = (("3P","3S","test one — unconfined"), ("2P","2S","test two — confined"))
    fig, axs = FigStyle.newfig(width = MM170, h = 5.2, nrows = 2, ncols = 4)
    letters = "abcdefgh"; li = 1
    for (r, (ft, st, label)) in enumerate(cols)
        TF, TS = ld(ft), ld(st)
        z = TF["z"]; m = FigStyle.archie_m(TF)
        aqs = aquifers(TF)
        zmax = 60.0; zi = findall(<=(zmax), z)

        # --- (a) porosity and SNMR water content ---
        ax = axs[r, 1]
        ax.plot(TF["θ_total"][zi], z[zi], color = C_TRUTH, lw = 1.1, label = "θ total")
        ax.plot(TF["θ_SMR"][zi],   z[zi], color = C_FRESH, lw = 1.1, label = "θ SNMR")
        ax.set_xlabel("water content (v/v)"); ax.set_ylabel("depth (m)")
        ax.set_xlim(0, 0.45)
        r == 1 && ax.legend(loc = "center right", bbox_to_anchor = (1.02, 0.62))

        # --- (b) clay fraction ---
        ax = axs[r, 2]
        ax.plot(TF["clay"][zi], z[zi], color = C_TRUTH, lw = 1.1, label = "clay fraction")
        ax.set_xlabel("clay fraction (v/v)"); ax.set_xlim(0, 1.0)

        # --- (c) pore-water salinity ---
        ax = axs[r, 3]
        # saline underneath, fresh on top: outside the aquifer the control is
        # identical to its parent, so whichever is drawn second is the one seen.
        ax.plot(TS["σw"][zi] .* MGL, z[zi], color = C_SALINE, lw = 1.6, label = "saline control")
        ax.plot(TF["σw"][zi] .* MGL, z[zi], color = C_FRESH,  lw = 1.1, label = "fresh")
        ax.axvline(1000, color = C_MGR, ls = "--", lw = 0.8, label = "1,000 mg/L")
        ax.set_xlabel("pore water (mg/L)"); ax.set_xscale("log"); ax.set_xlim(60, 6e3)
        r == 1 && ax.legend(loc = "lower right")

        # --- (d) bulk conductivity split into water and clay ---
        ax = axs[r, 4]
        clayterm = max.(TF["σ_bulk"][zi] .- TF["σw"][zi] .* TF["θ_total"][zi] .^ m, 1e-6)
        ax.fill_betweenx(z[zi], 1e-4, clayterm, color = "#555555", alpha = 0.30, lw = 0,
                         label = "clay term (shared)")
        for (T, cc, nm) in ((TS, C_SALINE, "saline"), (TF, C_FRESH, "fresh"))
            tot = T["σ_bulk"][zi]
            ax.fill_betweenx(z[zi], clayterm, tot, color = cc, alpha = 0.33, lw = 0,
                             label = "$nm, water term")
            ax.plot(tot, z[zi], color = cc, lw = 1.0, label = "$nm, total")
        end
        ax.set_xlabel("σ bulk (S/m)"); ax.set_xscale("log"); ax.set_xlim(1e-3, 3)
        r == 1 && ax.legend(loc = "lower left", fontsize = 6.0)

        for c in 1:4
            a = axs[r, c]
            for (z1, z2) in aqs; FigStyle.aquifer_band(a, z1, min(z2, zmax)); end
            a.set_ylim(zmax, 0)
            c > 1 && a.set_yticklabels([])
            FigStyle.panel(a, letters[li]); li += 1
        end
        axs[r, 1].text(-0.42, 0.5, label, transform = axs[r, 1].transAxes,
                       rotation = 90, va = "center", ha = "center", fontsize = 8)
    end
    fig.subplots_adjust(wspace = 0.12, hspace = 0.34, left = 0.10, right = 0.995,
                        top = 0.995, bottom = 0.095)
    FigStyle.finish(fig, "fig02_columns")
end
main()
