# sensitivity_factors.jl -- what each unknown costs the bound on test one.
#
# Replaces the S4 bound-decomposition waterfall. The waterfall walked from the
# OLD gate to the new one, mixing a method change with a column change; this
# asks a single question of the FINAL run instead: holding everything else
# fixed, how much tighter would the bound be if one unknown were known?
#
#   data-only bound  B0 = p95 over draws of  C_hi / K  * 6400          mg/L
#                    K  = sum over window cells of theta_draw^m_draw dz
#   clay known       numerator becomes C_hi - C_clay(truth)
#   porosity known   theta from the truth profile instead of SNMR draws
#   m known          m fixed at the truth's Archie exponent, not U(1.5,1.9)
#
# Each factor is bound_known / bound_data_only, so a number below 1 is the
# tightening that knowledge would buy.
using JLD2, Statistics, Printf, Random
const SYN = joinpath(@__DIR__, "..", "data")
const OUT = normpath(joinpath(@__DIR__,"..","outputs"))
const MGL = 6400.0
const ND  = 200_000

cellsof(z)=(mid=(z[1:end-1].+z[2:end])./2; (top=vcat(0.0,mid), bot=vcat(mid,z[end]+(z[end]-z[end-1])/2)))
function winw(z,a,b); g=cellsof(z); w=[max(0.0,min(g.bot[i],b)-max(g.top[i],a)) for i in eachindex(z)]
    ix=findall(>(1e-12),w); (ix=ix,w=w[ix]); end

function main(tag="3P", ens="3PT")
    T  = load(joinpath(SYN,"truth_murray$tag.jld2"))
    E  = load(joinpath(SYN,"cond_C2_typical$ens.jld2"))
    P4 = load(joinpath(OUT,"synth_$tag","stage4post_synth_$tag.jld2"))
    w  = (P4["window"].top, P4["window"].base)
    Chi = P4["band_C"].hi
    mlo, mhi = P4["m_prior"]; mtrue = T["build"]["m_archie"]

    # truth, on its own 0.25 m grid
    zt, dzt = T["z"], T["dz"]
    kt = findall(i -> zt[i] > w[1] && zt[i] <= w[2], eachindex(zt))
    clayterm = T["σ_s"][kt] .* T["clay"][kt]
    # check the rock-physics identity before using the pieces separately
    resid = maximum(abs.(T["σ_bulk"][kt] .- (T["σw"][kt] .* T["θ_total"][kt].^mtrue .+ clayterm)))
    @printf("rock-physics identity residual over the window: %.3e S/m\n", resid)
    C_clay = sum(clayterm) * dzt
    C_true = sum(T["σ_bulk"][kt]) * dzt
    @printf("window %.2f-%.2f m   band upper edge C_hi = %.4f S   true C = %.4f S\n",
            w..., Chi, C_true)
    @printf("true clay conductance over the window = %.4f S  (%.0f%% of true C)\n",
            C_clay, 100*C_clay/C_true)

    # K from SNMR draws, and K from the truth profile
    Θ  = Matrix{Float64}(reduce(hcat, E["ens"])); ws = winw(E["zsmr"], w...)
    wt = winw(zt, w...)
    Random.seed!(20261103)
    ms = mlo .+ (mhi-mlo) .* rand(ND)
    js = rand(1:size(Θ,2), ND)
    K_draw  = [sum((@view Θ[ws.ix,js[i]]) .^ ms[i] .* ws.w) for i in 1:ND]
    K_truth = [sum(T["θ_total"][wt.ix] .^ ms[i] .* wt.w) for i in 1:ND]
    # m known: same theta source, m pinned
    K_draw_mfix = [sum((@view Θ[ws.ix,js[i]]) .^ mtrue .* ws.w) for i in 1:ND]

    B(num, K) = quantile(num ./ K .* MGL, 0.95)
    B0    = B(Chi,          K_draw)
    Bclay = B(Chi - C_clay, K_draw)
    Bthe  = B(Chi,          K_truth)
    Bm    = B(Chi,          K_draw_mfix)

    println("\n" * "="^66)
    @printf("  %-34s %10s %8s\n", "bound on pore-water TDS", "mg/L", "factor")
    @printf("  %-34s %10.0f %8s\n", "data only (the published bound)", B0, "1.00")
    @printf("  %-34s %10.0f %7.2fx\n", "clay term known",   Bclay, Bclay/B0)
    @printf("  %-34s %10.0f %7.2fx\n", "porosity known",    Bthe,  Bthe/B0)
    @printf("  %-34s %10.0f %7.2fx\n", "m known",           Bm,    Bm/B0)
    println("="^66)
    @printf("  true pore water %.0f mg/L\n", MGL * mean(T["σw"][kt]))
end
main()
