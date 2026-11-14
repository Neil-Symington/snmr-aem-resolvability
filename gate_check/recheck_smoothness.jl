# Apply both smoothness tests to every stored profile. No solves.
using JLD2, Printf
include(joinpath(@__DIR__, "smoothness.jl")); using .Smoothness
for f in sort(filter(x -> startswith(x,"profile_") && endswith(x,".jld2"),
                     readdir(joinpath(@__DIR__,"results"))))
    occursin("smoke", f) && continue
    D = load(joinpath(@__DIR__,"results",f))
    x = D["lt"]; y = D["min_chi2"]; C = 10 .^ x .* D["thick"]
    r = both_flags(x, y)
    @printf("\n=== %s  (N=%d, threshold %.2f)\n", f, D["ndata"], D["ndata"]+3.84)
    @printf("%-5s %10s %11s %11s %7s %9s\n","i","C (S)","min chi2","quad pred","spike","curvature")
    for i in eachindex(y)
        @printf("%-5d %10.4f %11.3f %11s %7s %9s\n", i, C[i], y[i],
                isfinite(r.predicted[i]) ? @sprintf("%.3f", r.predicted[i]) : "-",
                r.spike[i] ? "FLAG" : "", r.curvature[i] ? "FLAG" : "")
    end
    @printf("  spike flags: %s\n", any(r.spike) ? string(findall(r.spike)) : "none")
    @printf("  curvature flags: %s\n", any(r.curvature) ? string(findall(r.curvature)) : "none")
end
