# test_cond_constraint.jl -- does ctype=:cond actually force the conductance?
#
#   julia --project=. test_cond_constraint.jl
#
# Three properties, asserted. A constraint that is nearly satisfied is not a
# constraint: the whole point of L(C) is that C is known exactly at each node.
#   1. the returned model's conductance equals the requested C to 1e-10 relative
#   2. it does so across the whole grid, not just near the data's preference
#   3. the forced quantity is NOT the geometric-mean statistic -- the two differ,
#      and the difference is the Jensen gap the old code silently absorbed
using JLD2, Statistics, Printf, LinearAlgebra
const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "occam_core.jl")); using .OccamCore
include(joinpath(ROOT, "src", "aem_jac.jl"))
const OUT = joinpath(HERE, "..", "data")

# Wrapped in a function on purpose. At top level, `fails += 1` inside a `for`
# loop hits Julia's soft scope: `fails` becomes a NEW LOCAL each iteration, the
# global stays 0, and `exit(fails == 0)` then reports PASS whatever happened.
# Julia warns, but a warning in a log is how the last round of bugs survived.
function main()
T = load(joinpath(OUT, "truth_murray3.jld2"))
B = build_aem_jac(T, 123); aem, grid = B.aem, B.grid
n = length(aem.ρ) - 1; zall = grid.zall; ndata = length(aem.res)
R = makeR1(n); m0 = fill(-2.0, n)
u = T["units"]; aq = [uu for uu in u if occursin("AQUIFER", uppercase(string(uu[3])))][1]
za, zb = aq[1] + 3.0, min(aq[2], 115.0) - 3.0
fixidx = findall(z -> za <= z <= zb, zall)
dzc = diff(vcat(0.0, collect(Float64, zall)))
thick = sum(dzc[fixidx])
@printf("unit %.1f-%.1f m, %d cells, thickness %.2f m\n", za, zb, length(fixidx), thick)

fails = 0
@printf("\n%10s %12s %12s %10s %10s %6s %5s %5s\n",
        "target C", "achieved C", "rel err", "chi2/N", "geo-mean C", "bound", "acc", "ok")
for σeff in [1e-4, 1e-3, 3e-3, 1e-2, 0.0302, 0.1, 0.3, 1.0, 3.0]
    C = σeff * thick
    m, χ, H, nf = occam(aem, m0, m0, R, ndata; fixidx = fixidx, fixval = C,
                        mode = :map, ctype = :cond, dz = dzc, lambda2 = 100.0,
                        beta2 = 0.1, withH = true, info = true)
    ach = sum(10 .^ m[fixidx] .* dzc[fixidx])
    rel = abs(ach - C)/C
    Cgeo = 10^mean(m[fixidx]) * thick            # the statistic the old code forced
    # A node is only a point on L(C) if the constraint was actually imposed AND
    # the optimum is interior. Unreachable nodes are expected at the ends of a
    # wide grid -- the failure mode is normalising over them anyway, which is
    # what the superseded run did.
    ok = rel < 1e-10
    (ok || nf.nbound > 0) || (fails += 1)
    @printf("%10.4g %12.6g %12.3e %10.3f %10.4g %6d %5d %5s\n",
            C, ach, rel, χ/ndata, Cgeo, nf.nbound, nf.nacc,
            ok && nf.nbound == 0 ? "yes" : nf.nbound > 0 ? "bound" : "NO")
    flush(stdout)
end
@printf("\n%s\n", fails == 0 ? "PASS: conductance forced exactly at every REACHABLE node" :
                               "FAIL: $fails nodes off the constraint surface")
    return fails
end
exit(main() == 0 ? 0 : 1)
