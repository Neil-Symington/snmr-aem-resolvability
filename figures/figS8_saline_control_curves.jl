# fig07_saline_controls.jl -- Figure 7. The saline controls.
#
# One panel per column pair: the Laplace L(C) at lambda-hat for the fresh column
# and its saline control, each truth marked, and the 1,000 mg/L line converted to
# conductance. Test one separates; test two does not. That contrast is the
# figure's whole purpose, so both panels share an x-range.
using JLD2, Statistics, Printf, PyPlot, Random
include("figstyle.jl"); using .FigStyle

# Verify every label against the column Stage 0 sealed. These scripts were
# left unwired when the others were done, and figS8 was one of three a
# reviewer drew from fabricated columns with no message at all.
foreach(FigStyle.assert_column, ("3P", "3S", "2P", "2S"))
const SYN = normpath(joinpath(@__DIR__, "..", "data"))
const OUT = normpath(joinpath(@__DIR__,"..","outputs"))
const MGL = 6400.0
cellsof(z)=(mid=(z[1:end-1].+z[2:end])./2; (top=vcat(0.0,mid), bot=vcat(mid,z[end]+(z[end]-z[end-1])/2)))
function winw(z,a,b); g=cellsof(z); w=[max(0.0,min(g.bot[i],b)-max(g.top[i],a)) for i in eachindex(z)]
    ix=findall(>(1e-12),w); (ix=ix,w=w[ix]); end

function curve(tag)
    G = load(joinpath(OUT,"synth_$tag","stage4_synth_$tag.jld2"))
    SW = G["sweeps"][string(G["lambda_hat"])]
    lL = -SW.phi./2 .- SW.ldH./2 .+ SW.ldP./2; lL[.!SW.good] .= NaN
    k = findall(isfinite, lL)
    pv = fill(NaN, length(lL)); isempty(k) || (pv[k] = exp.(lL[k] .- maximum(lL[k])))
    (G["Cg"], pv ./ maximum(filter(isfinite, pv)), G)
end
function truthC(tag, w)
    T = load(joinpath(SYN,"truth_murray$tag.jld2")); zt, sbt, dzt = T["z"], T["σ_bulk"], T["dz"]
    sum(sbt[findall(k -> zt[k] > w[1] && zt[k] <= w[2], eachindex(zt))] .* dzt)
end

function main()
    pairs = (("3P","3S","3PT","test one — unconfined"), ("2P","2S","2PT","test two — confined"))
    fig, axs = FigStyle.newfig(width = MM170, h = 3.3, nrows = 1, ncols = 2)
    for (i,(ft,st,ens,label)) in enumerate(pairs)
        ax = axs[i]
        for (tag, cc, nm) in ((ft, C_FRESH, "fresh ($ft)"), (st, C_SALINE, "saline control ($st)"))
            Cg, pv, G = curve(tag)
            w = (G["window"].top, G["window"].base)
            ax.plot(Cg, pv, color = cc, lw = 1.3, label = "L(C), $nm")
            ax.axvline(truthC(tag, w), color = cc, ls = "-", lw = 0.9, alpha = 0.55)
            ax.plot([], [], color = cc, lw = 0.9, alpha = 0.55, label = "truth, $nm")
        end
        # the 1,000 mg/L line in conductance, at the fresh column's median theta, m
        G = load(joinpath(OUT,"synth_$ft","stage4_synth_$ft.jld2"))
        w = (G["window"].top, G["window"].base)
        E = load(joinpath(SYN,"cond_C2_typical$ens.jld2")); ww = winw(E["zsmr"], w...)
        Θ = Matrix{Float64}(reduce(hcat, E["ens"])); Random.seed!(1)
        ms = 1.5 .+ 0.4 .* rand(4000); js = rand(1:size(Θ,2), 4000)
        K = median([sum((@view Θ[ww.ix,js[t]]) .^ ms[t] .* ww.w) for t in 1:4000])
        ax.axvline(1000/MGL*K, color = C_MGR, ls = "--", lw = 0.9, label = "1,000 mg/L")
        ax.set_xscale("log"); ax.set_xlim(2e-2, 20); ax.set_ylim(0, 1.12)
        ax.set_xlabel("unit conductance C (S)")
        i == 1 ? ax.set_ylabel("L(C) at λ̂, normalised") : ax.set_yticklabels([])
        FigStyle.panel(ax, "ab"[i]); ax.set_title(label, fontsize = 8, pad = 3)
        ax.legend(loc = i == 1 ? "upper right" : "center left", fontsize = 6.0)
    end
    fig.subplots_adjust(wspace = 0.07, left = 0.08, right = 0.995, top = 0.925, bottom = 0.145)
    FigStyle.finish(fig, "figS8_saline_control_curves")
end
main()
