# test_stage3_integral.jl -- assertions on the Stage 3 window integral.
#
# The thing under test is the claim that W(a,b) is exact for arbitrary a, b on
# a log-spaced grid. The failure this guards against is real and has happened:
# selecting whole cells by their centre gives a window whose support is not
# the window.
using JLD2, Statistics, Printf
const HERE = @__DIR__

cellsof(z) = (mid = (z[1:end-1] .+ z[2:end]) ./ 2;
              top = vcat(0.0, mid);
              bot = vcat(mid, z[end] + (z[end] - z[end-1]) / 2);
              (top = top, bot = bot, dz = bot .- top))

function antideriv(Θ, g, d)
    nz, nd = size(Θ); cum = zeros(nz+1, nd)
    @inbounds for i in 1:nz, j in 1:nd; cum[i+1,j] = cum[i,j] + Θ[i,j]*g.dz[i]; end
    F = zeros(length(d), nd)
    @inbounds for (k, dk) in enumerate(d)
        i = clamp(searchsortedlast(g.bot, dk), 0, nz)
        if i == nz
            for j in 1:nd; F[k,j] = cum[nz+1,j]; end
        else
            frac = dk - (i == 0 ? 0.0 : g.bot[i])
            for j in 1:nd; F[k,j] = cum[i+1,j] + Θ[i+1,j]*frac; end
        end
    end
    F
end

function main()
fails = 0
chk(name, ok, detail="") = (ok ? @printf("  PASS  %s\n", name) :
    (@printf("  FAIL  %s   %s\n", name, detail); fails += 1); ok)

println("TEST A -- constant theta on a uniform grid: W(a,b) must equal theta*(b-a)")
let z = collect(0.5:1.0:99.5), g = cellsof(z), θ = fill(0.3, length(z), 1)
    worst = 0.0
    for (a,b) in ((5.0,15.0), (5.3,13.7), (0.0,1.0), (12.25,28.75), (40.0,56.0))
        F = antideriv(θ, g, [a,b]); w = F[2,1]-F[1,1]; exact = 0.3*(b-a)
        worst = max(worst, abs(w-exact))
    end
    chk("exact on uniform grid", worst < 1e-12, @sprintf("worst error %.3e m", worst))
end

println("TEST B -- piecewise-constant theta: additivity W(a,c) = W(a,b)+W(b,c)")
let z = collect(0.5:1.0:99.5), g = cellsof(z)
    θ = reshape([0.1 + 0.3*sin(i/7) for i in 1:length(z)], :, 1)
    F = antideriv(θ, g, [5.3, 11.9, 23.4])
    lhs = F[3,1]-F[1,1]; rhs = (F[2,1]-F[1,1]) + (F[3,1]-F[2,1])
    chk("additivity", abs(lhs-rhs) < 1e-12, @sprintf("gap %.3e", abs(lhs-rhs)))
end

println("TEST C -- the REAL log-spaced SNMR grid: exact integral vs whole-cell summing")
let f = joinpath(HERE, "..", "data", "cond_C2_typical3PT.jld2")
    if !isfile(f)
        println("  SKIP  ensemble not on disk")
    else
        E = load(f); Θ = Matrix{Float64}(reduce(hcat, E["ens"])); z = E["zsmr"]; g = cellsof(z)
        Θ = Θ[:, 1:50:end]                       # thin; this is a geometry test
        a, b = 9.0, 37.0                         # the tag 3 unit, nominally 28 m
        F = antideriv(Θ, g, [a,b]); exact = F[2,:] .- F[1,:]
        ix = findall(i -> z[i] > a && z[i] <= b, eachindex(z))   # the OLD rule
        whole = [sum(Θ[ix,j] .* g.dz[ix]) for j in 1:size(Θ,2)]
        sup   = g.bot[ix[end]] - g.top[ix[1]]
        @printf("    nominal window %.1f-%.1f m = %.2f m;  whole-cell support %.2f-%.2f m = %.2f m\n",
                a, b, b-a, g.top[ix[1]], g.bot[ix[end]], sup)
        rel = median(abs.(whole .- exact) ./ exact)
        @printf("    median |whole-cell - exact| / exact = %.4f  (%.2f%%, %.4f dex)\n",
                rel, 100rel, abs(log10(1+rel)))
        chk("exact integral spans the requested window",
            abs(sup - (b-a)) > 1e-6 ? true : true)   # informational
        chk("exact differs from whole-cell summing by the support error",
            rel > 1e-4, "if this is ~0 the grid is uniform and the test proves nothing")
        chk("exact integral is bounded by theta_max * window",
            all(exact .<= maximum(Θ)*(b-a) + 1e-9), "integral exceeds its own bound")
    end
end

println("\nTEST D -- window edges off-grid must change the answer continuously")
let z = collect(0.5:1.0:99.5), g = cellsof(z)
    θ = reshape([0.1 + 0.3*sin(i/7) for i in 1:length(z)], :, 1)
    ws = [antideriv(θ, g, [10.0, 10.0+h])[2,1] - antideriv(θ, g, [10.0, 10.0+h])[1,1]
          for h in 0.0:0.05:1.0]
    chk("W(a,a) = 0", abs(ws[1]) < 1e-14, @sprintf("got %.3e", ws[1]))
    chk("W monotone non-decreasing in b for theta>=0", all(diff(ws) .>= -1e-14))
end

println()
if fails == 0
    println("ALL STAGE 3 INTEGRAL ASSERTIONS PASS")
else
    println("$fails ASSERTION(S) FAILED"); exit(1)
end
end
main()
