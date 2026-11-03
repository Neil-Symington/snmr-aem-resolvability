# dc_operator.jl -- a trans-D operator for the Schlumberger sounding.
#
# HiQGA ships no DC physics, so this is the FooPhysicsInversion pattern filled
# in with the Ghosh filter forward of dc_forward.jl: a subtype of Operator1D
# with a `get_misfit(m, opt, F)` method. Wrapped beside the SkyTEM operator in
# transD_GP's GenericJointInversion.JointOperator, whose get_misfit simply sums
# the two, the sounding becomes a second data term in the same trans-D chain.
#
# SPEED. The layer thicknesses are FIXED (the AEM grid), only the resistivities
# move, so every tanh(lambda_i * h_k) is precomputed once: 201 filter abscissae
# x 40 radii x 64 interfaces. A forward is then 515k recursion steps and costs
# about 4 ms against the AEM forward's 56 ms.
using HiQGA, HiQGA.transD_GP, LinearAlgebra
import HiQGA.transD_GP.AbstractOperator.get_misfit
include(joinpath(@__DIR__, "..", "src", "dc_forward.jl"))

mutable struct DCSchlum <: HiQGA.transD_GP.AbstractOperator.Operator1D
    d     :: Vector{Float64}       # observed apparent resistivity, Ohm.m
    σ     :: Vector{Float64}       # per-point sigma, Ohm.m
    W0    :: Vector{Float64}       # J0 filter weights
    TANH  :: Array{Float64,3}      # [nfilt, nradius, ninterface]
    invr  :: Vector{Float64}       # 1/r per radius
    pref  :: Vector{Float64}       # (L^2-l^2)/(2l) per spacing
    pred  :: Vector{Float64}       # workspace: predicted rho_a
    ndata :: Int
    useML :: Bool
end

"""Build the operator. `ztop` are the layer tops on the AEM grid (z[2:end]),
so there are length(ztop) layers and length(ztop)-1 interfaces."""
function DCSchlum(d, σ, L, l, ztop; useML = false)
    h  = diff(collect(Float64, ztop))
    nf = length(DCForward.Filter_base)
    r  = Float64[]
    for j in eachindex(L); push!(r, L[j]-l[j]); push!(r, L[j]+l[j]); end
    TANH = Array{Float64,3}(undef, nf, length(r), length(h))
    @inbounds for (jr, rj) in enumerate(r), i in 1:nf
        lam = DCForward.Filter_base[i]/rj
        for k in eachindex(h); TANH[i,jr,k] = tanh(lam*h[k]); end
    end
    DCSchlum(collect(Float64,d), collect(Float64,σ), copy(DCForward.Filter_J0),
             TANH, 1 ./ r, [(L[j]^2-l[j]^2)/(2l[j]) for j in eachindex(L)],
             zeros(length(d)), length(d), useML)
end

"""Apparent resistivity for log10-resistivity vector `m` on the AEM grid."""
function forward!(F::DCSchlum, m::AbstractArray)
    nf, nr, ni = size(F.TANH)
    ρ = 10 .^ m                                   # layer resistivities
    @inbounds for jr in 1:nr
        acc = 0.0
        for i in 1:nf
            T = ρ[ni+1]
            for k in ni:-1:1
                t = F.TANH[i,jr,k]
                T = (T + ρ[k]*t) / (1 + T*t/ρ[k])
            end
            acc += T * F.W0[i]
        end
        s = (jr + 1) >> 1
        jr % 2 == 1 ? (F.pred[s]  =  acc*F.invr[jr]) :
                      (F.pred[s] -=  acc*F.invr[jr])
    end
    @inbounds for s in eachindex(F.pred); F.pred[s] *= F.pref[s]; end
    F.pred
end

get_misfit(m::HiQGA.transD_GP.Model, opt::HiQGA.transD_GP.Options, F::DCSchlum) =
    get_misfit(m.fstar, opt, F)

function get_misfit(m::AbstractArray, opt::HiQGA.transD_GP.Options, F::DCSchlum)
    opt.debug && return 0.0
    forward!(F, m)
    r = (F.pred .- F.d) ./ F.σ
    F.useML ? 0.5*F.ndata*log(norm(r)^2) : 0.5*norm(r)^2
end

"chi2/N for the current model, for reporting."
chi2byN(F::DCSchlum, m) = (forward!(F, m); sum(((F.pred .- F.d)./F.σ).^2)/F.ndata)
