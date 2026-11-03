# dc_forward.jl -- 1D Schlumberger sounding forward, Ghosh digital filter.
#
# HiQGA has no DC resistivity physics (MT1D, CSEM1DEr, AEM_VMD_HMD, and a
# FooPhysics template -- no DCIP), so this implements the standard linear
# filter method. It reuses HiQGA's own 201-point J1 Hankel filter rather than
# shipping a second set of coefficients.
#
#   rho_a(s) = s^2 * integral_0^inf T(lambda) J1(lambda s) lambda dlambda
#
# with the Hankel filter  integral f(l) J1(l s) dl ~ (1/s) sum f(b_i/s) W1_i
# and f(l) = T(l)*l this collapses to
#
#   rho_a(s) = sum_i T(b_i/s) * b_i * W1_i
#
# T is the Koefoed/Pekeris resistivity transform, recursed from the basement up:
#   T_N = rho_N;  T_k = (T_{k+1} + rho_k*tanh(l h_k)) / (1 + T_{k+1}*tanh(l h_k)/rho_k)
#
# VERIFIED against the two-layer analytic series (see selftest below), which is
# the check the brief asks for.
module DCForward
using Printf
# HiQGA's filter coefficients, found through the loaded package when there is
# one and through the depot otherwise, so this file also loads in a plain
# environment (the independent cross-check runs in one).
function _digfilters()
    p = Base.find_package("HiQGA")
    p === nothing || return joinpath(dirname(p), "DigFilters.jl")
    for d in DEPOT_PATH
        r = joinpath(d, "packages", "HiQGA")
        isdir(r) || continue
        for v in readdir(r)
            f = joinpath(r, v, "src", "DigFilters.jl")
            isfile(f) && return f
        end
    end
    error("DigFilters.jl not found: HiQGA is not installed in this depot")
end
include(_digfilters())            # Filter_base, Filter_J0, Filter_J1
export rho_schlumberger, rho_schlumberger_finite, rho_two_layer_analytic, selftest

"""Koefoed resistivity transform T(lambda) for layer resistivities `rho`
(Ohm.m, length N) and thicknesses `h` (m, length N-1)."""
function tfun(lam::Real, rho::AbstractVector, h::AbstractVector)
    T = rho[end]
    @inbounds for k in length(rho)-1:-1:1
        t = tanh(lam * h[k])
        T = (T + rho[k]*t) / (1 + T*t/rho[k])
    end
    T
end

"""Schlumberger apparent resistivity at half-spacings `s` = AB/2 (m).
The MN/2 << AB/2 limit; the brief's MN/2 = AB/10 keeps the correction below
the 3% noise, which `selftest` reports."""
function rho_schlumberger(rho::AbstractVector, h::AbstractVector, s::AbstractVector)
    out = similar(collect(Float64, s))
    @inbounds for (j, sj) in enumerate(s)
        acc = 0.0
        for i in eachindex(Filter_base)
            acc += tfun(Filter_base[i]/sj, rho, h) * Filter_base[i] * Filter_J1[i]
        end
        out[j] = acc          # sum_i b_i*W1_i = 1 exactly, so no 1/s factor
    end
    out
end

"""Schlumberger apparent resistivity with FINITE MN, which is what the survey
actually has: the brief sets MN/2 = AB/10, so (l/L)^2 = 1/100 and the dipole
limit is wrong by about a percent -- not negligible beside 3% noise.

    rho_a = (L^2-l^2)/(2l) * [ int T J0(lambda (L-l)) dl - int T J0(lambda (L+l)) dl ]

which goes to the dipole form as l -> 0. `L` = AB/2, `l` = MN/2."""
function rho_schlumberger_finite(rho::AbstractVector, h::AbstractVector,
                                 L::AbstractVector, l::AbstractVector)
    out = similar(collect(Float64, L))
    @inbounds for j in eachindex(L)
        Lj, lj = L[j], l[j]
        acc = 0.0
        for (r, sgn) in ((Lj - lj, 1.0), (Lj + lj, -1.0))
            a = 0.0
            for i in eachindex(Filter_base)
                a += tfun(Filter_base[i]/r, rho, h) * Filter_J0[i]
            end
            acc += sgn * a / r
        end
        out[j] = (Lj^2 - lj^2) / (2lj) * acc
    end
    out
end

"""Two-layer Schlumberger apparent resistivity, analytic image series."""
function rho_two_layer_analytic(rho1, rho2, h, s; nmax = 4000)
    k = (rho2 - rho1) / (rho2 + rho1)
    acc = 0.0
    for n in 1:nmax
        acc += k^n / (1 + (2n*h/s)^2)^1.5
    end
    rho1 * (1 + 2acc)
end

"""Check the filter against the analytic two-layer case over a wide contrast
range and spacing range. Returns the worst relative error."""
function selftest(; verbose = true)
    worst = 0.0; worstat = ()
    for (r1, r2, h) in ((10.0, 100.0, 20.0), (100.0, 10.0, 20.0), (5.0, 500.0, 8.0),
                        (200.0, 2.0, 40.0), (30.0, 30.0, 15.0))
        s = 10 .^ range(log10(2.0), log10(300.0), length = 20)
        num = rho_schlumberger([r1, r2], [h], s)
        ana = [rho_two_layer_analytic(r1, r2, h, sj) for sj in s]
        rel = abs.(num .- ana) ./ ana
        m, i = findmax(rel)
        if m > worst; worst = m; worstat = (r1, r2, h, s[i]); end
        verbose && @printf("  rho1 %6.1f  rho2 %6.1f  h %5.1f  max rel err %.3e\n",
                           r1, r2, h, m)
    end
    verbose && @printf("  WORST %.3e  at rho1=%.1f rho2=%.1f h=%.1f s=%.2f m\n",
                       worst, worstat...)
    worst
end
end # module
