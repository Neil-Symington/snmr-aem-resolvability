# fig03_aem_posteriors.jl -- Figure 3. AEM conductivity posteriors (band form).
#
# FINAL FORM: p10-p90 band and p50 from the Stage 1 trans-D posterior, truth,
# aquifer shaded, and a narrow strip giving the Stage 1 resolution. The 2D
# density variants are stood down.
#
# Resolution is `dlog_for_1sigma` from aem_sens: the change in log10 sigma at
# that cell needed to move the data by one noise standard deviation. Small =
# well resolved. It has its own strip rather than a wash over the data, so it
# carries the number as well as the impression.
#
# Draws are the FULL 2002 post-burn-in T=1 models per column
# (aem_draws_full_*.jld2, written by synthetics/extract_full_draws.jl), not the
# 600-draw subsample. The p50 agrees between the two to 0.045 dex on test one
# and 0.072 dex on test two, so this is more samples, not a different answer.
using JLD2, Statistics, Printf, PyPlot
include("figstyle.jl"); using .FigStyle
const SYN = normpath(joinpath(@__DIR__, "..", "data"))

function aquifers(T)
    acc = Tuple{Float64,Float64}[]
    for u in T["units"]
        k = u isa Dict ? u : Dict(string(f)=>getfield(u,f) for f in propertynames(u))
        occursin("AQUIFER", uppercase(string(k["3"]))) && push!(acc,(Float64(k["1"]),Float64(k["2"])))
    end
    acc
end

function main()
    # Panel (c) appears as soon as the joint draws and the joint sensitivity
    # exist, so this script builds either way and needs no edit when they land.
    cols = Any[("3P","3M","test one — unconfined"), ("2P","2M","test two — confined")]
    (isfile(joinpath(SYN,"aem_draws_full_2PJ.jld2")) &&
     isfile(joinpath(SYN,"aem_sens_2PJ.jld2"))) &&
        push!(cols, ("2PJ","2PJ","test two + sounding"))
    np = length(cols)
    fig = figure(figsize = (MM170, 3.5))
    gs = fig.add_gridspec(1, 2np; width_ratios = repeat([1, 0.17], np), wspace = 0.07)
    zmax = 80.0
    for (i,(tag, src, label)) in enumerate(cols)
        ax  = fig.add_subplot(get(gs, (0, 2(i-1))))
        axr = fig.add_subplot(get(gs, (0, 2(i-1)+1)))
        T = load(joinpath(SYN, "truth_murray$(replace(tag, "J" => "")).jld2"))
        Σ = Matrix{Float64}(load(joinpath(SYN,"aem_draws_full_$src.jld2"), "sigma"))
        zb = collect(Float64, load(joinpath(SYN,"aem_draws_full_$src.jld2"), "zall_aem"))
        S = load(joinpath(SYN, "aem_sens_$tag.jld2"))
        d = S["dlog_for_1sigma"]; zs = S["zall"]; dref = S["dlog_ref"]

        for (z1,z2) in aquifers(T)
            FigStyle.aquifer_band(ax, z1, min(z2,zmax)); FigStyle.aquifer_band(axr, z1, min(z2,zmax))
        end

        # cells that OVERLAP the plotted range, so the curves reach the axis
        # bottom instead of stopping at the last cell above 80 m
        ctop = vcat(0.0, zb[1:end-1])
        k = findall(j -> ctop[j] < zmax, eachindex(zb))
        p10 = [quantile(Σ[j,:],0.10) for j in k]
        p50 = [quantile(Σ[j,:],0.50) for j in k]
        p90 = [quantile(Σ[j,:],0.90) for j in k]
        ax.fill_betweenx(zb[k], p10, p90, color = C_FRESH, alpha = 0.28, lw = 0,
                         label = "AEM p10–p90")
        ax.step(p50, zb[k], where = "pre", color = C_FRESH, lw = 1.1, label = "AEM p50")
        zt = T["z"]; kt = findall(<=(zmax + 2), zt)
        ax.plot(T["σ_bulk"][kt], zt[kt], color = C_TRUTH, lw = 1.0, label = "truth")

        ax.set_xscale("log"); ax.set_xlim(1e-3, 3); ax.set_ylim(zmax, 0)
        ax.set_xlabel("σ bulk (S/m)")
        i == 1 ? ax.set_ylabel("depth (m)") : ax.set_yticklabels([])
        FigStyle.panel(ax, "abc"[i]); ax.set_title(label, fontsize = 8, pad = 3)
        i == 1 && ax.legend(loc = "lower right", fontsize = 6.4)

        kk = findall(j -> zs[j] <= zmax, eachindex(zs))
        axr.fill_betweenx(zs[kk], 1e-3, d[kk], color = "#b03030", alpha = 0.45, lw = 0)
        axr.plot(d[kk], zs[kk], color = "#b03030", lw = 0.8)
        axr.axvline(dref, color = C_MGR, ls = "--", lw = 0.7)
        axr.set_xscale("log"); axr.set_xlim(1e-2, 3); axr.set_ylim(zmax, 0)
        axr.set_yticklabels([]); axr.set_xticks([1e-2, 1])
        axr.set_xticklabels(["0.01","1"], fontsize = 6)
        axr.set_xlabel("Δlog σ\nfor 1σ", fontsize = 6.5, labelpad = 1)
        axr.tick_params(axis = "y", length = 0)
    end
    fig.subplots_adjust(left = 0.065, right = 0.995, top = 0.925, bottom = 0.175)
    FigStyle.finish(fig, "fig03_aem_posteriors")
end
main()
