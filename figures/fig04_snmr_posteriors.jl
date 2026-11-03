# fig04_snmr_posteriors.jl -- Figure 4. SNMR water-content posteriors (band form).
#
# FINAL FORM: p10-p90 band and p50 of theta_SMR from the Stage 2 trans-D
# posterior, the truth's theta_SNMR, the Stage 2 depth limit marked, aquifer
# shaded. The 2D density variants are stood down.
#
# THE DEPTH LIMIT IS PER COLUMN, NOT PER TIER -- 59.2 m on test one against
# 49.4 m on test two, both at tier "typical". The SNMR kernel is
# conductivity-sensitive, so each column's own conductivity sets it. The 70.1 m
# in stage1_typical.jld2 belongs to the BARE lineage on an 88.6 m loop and is
# not these columns; an earlier draft of this figure carried it in error.
#
# It is the SENSITIVITY DOI: per cell the smallest detectable water-content
# change is dmin = 1/sqrt(sum_q (w|K|/sigma-hat)^2), and the limit is where dmin
# crosses theta_detect = 0.05 v/v -- Christiansen & Auken in SNMR units, computed
# from the kernel and the noise alone. It is NEVER the trans-D posterior-width
# DOI, which tracks nucleus count rather than depth sensitivity and must not be
# quoted (doi_study/README.md).
using JLD2, Statistics, Printf, PyPlot
include("figstyle.jl"); using .FigStyle
const SYN = normpath(joinpath(@__DIR__, "..", "data"))
cellsof(z)=(mid=(z[1:end-1].+z[2:end])./2; (top=vcat(0.0,mid), bot=vcat(mid,z[end]+(z[end]-z[end-1])/2)))

function aquifers(T)
    acc = Tuple{Float64,Float64}[]
    for u in T["units"]
        k = u isa Dict ? u : Dict(string(f)=>getfield(u,f) for f in propertynames(u))
        occursin("AQUIFER", uppercase(string(k["3"]))) && push!(acc,(Float64(k["1"]),Float64(k["2"])))
    end
    acc
end

"per-column sensitivity DOI: where the smallest detectable dtheta crosses 0.05 v/v"
function sensitivity_doi(S; θdetect = 0.05)
    K = S["kernel"]; dz = S["kernel_dz"]; z = S["zsmr"]; sh = S["typical_σ̂"]
    nz, nq = size(K)
    w = [0.5*((i>1 ? dz[i-1] : 0.0) + (i<nz ? dz[i] : 0.0)) for i in 1:nz]
    dmin = [1/sqrt(sum((w[i]*abs(K[i,q])/sh[q])^2 for q in 1:nq)) for i in 1:nz]
    i = findlast(<=(θdetect), dmin)
    i === nothing ? 0.0 : Float64(z[i])
end

function main()
    cols = (("3P","3PT","test one — unconfined"), ("2P","2PT","test two — confined"))
    fig, axs = FigStyle.newfig(width = MM170, h = 3.4, nrows = 1, ncols = 2)
    zmax = 80.0
    for (i,(tag, ens, label)) in enumerate(cols)
        ax = axs[i]
        T = load(joinpath(SYN, "truth_murray$tag.jld2"))
        E = load(joinpath(SYN, "cond_C2_typical$ens.jld2"))
        S = load(joinpath(SYN, "smr_data_murray$ens.jld2"))
        Θ = Matrix{Float64}(reduce(hcat, E["ens"])); zs = E["zsmr"]
        g = cellsof(zs); zlim = sensitivity_doi(S)

        for (z1,z2) in aquifers(T); FigStyle.aquifer_band(ax, z1, min(z2,zmax)); end
        k = findall(j -> g.top[j] < zmax, eachindex(zs))
        ax.fill_betweenx(zs[k], [quantile(Θ[j,:],0.10) for j in k],
                                [quantile(Θ[j,:],0.90) for j in k],
                         color = C_FRESH, alpha = 0.28, lw = 0, label = "θ SNMR p10–p90")
        ax.plot([quantile(Θ[j,:],0.50) for j in k], zs[k], color = C_FRESH, lw = 1.1,
                label = "θ SNMR p50")
        zt = T["z"]; kt = findall(<=(zmax + 2), zt)
        ax.plot(T["θ_SMR"][kt], zt[kt], color = C_TRUTH, lw = 1.0, label = "truth θ SNMR")
        ax.axhline(zlim, color = C_MGR, ls = "--", lw = 1.0,
                   label = @sprintf("sensitivity DOI, %.1f m", zlim))

        ax.set_xlim(0, 0.45); ax.set_ylim(zmax, 0)
        ax.set_xlabel("θ SNMR (v/v)")
        i == 1 ? ax.set_ylabel("depth (m)") : ax.set_yticklabels([])
        FigStyle.panel(ax, "ab"[i]); ax.set_title(label, fontsize = 8, pad = 3)
        # a legend per panel, because the DOI differs between them
        ax.legend(loc = "lower right", fontsize = 6.2)
    end
    fig.subplots_adjust(wspace = 0.08, left = 0.075, right = 0.995, top = 0.925, bottom = 0.14)
    FigStyle.finish(fig, "fig04_snmr_posteriors")
end
main()
