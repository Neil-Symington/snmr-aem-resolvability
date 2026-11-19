# test_edgefinder.jl -- the edge finder against a MOCK solver, no AEM code.
#
#   julia --project=. test/test_edgefinder.jl
#
# The mock has the one property that makes bisection necessary: the chi2 it
# returns at a node is inflated by how far the warm start was from that node.
# So a coarse grid over-reports chi2 near the edge and a refined one does not,
# which is the behaviour measured on the real AEM problem (47.878 at 0.1 dex
# against 42.943 at 0.05 dex, across the threshold).
#
# These cases would have caught both bugs found on 17 Nov 2026.
using Printf
include(joinpath(@__DIR__, "..", "gate_check", "edgefinder.jl")); using .EdgeFinder
const EF = EdgeFinder

const THR   = 45.84
const XKNEE = -1.6
const TRUE_EDGE = XKNEE + sqrt((THR - 25.0)/1000)     # f(x) = 25 + 1000 (x-XKNEE)^2
f(x) = 25.0 + 1000.0 * max(0.0, x - XKNEE)^2
const PEN = 60.0                                       # chi2 added per dex of start distance

nfail = 0
function check(name, cond, detail = "")
    global nfail
    cond || (nfail += 1)
    @printf("  %-4s %s%s\n", cond ? "ok" : "FAIL", name, isempty(detail) ? "" : "  [$detail]")
end

"mock solver: model is the abscissa it was solved at; distance inflates chi2"
mock(x, from) = (chi = f(x) + PEN * abs(x - from), m = x)

@printf("edge finder tests; true crossing at %.5f dex (threshold %.2f)\n\n", TRUE_EDGE, THR)

# ── 1. recovers the edge from a coarse grid, and converges ────────────────
println("1. coarse grid, 0.2 dex")
S = Dict{Float64,NamedTuple}()
for x in -2.4:0.2:-1.2; S[x] = (chi = f(x) + PEN*0.2, m = x); end   # as if each solved from a neighbour
r = EF.bisect_edge(S, mock, THR; maxref = 10)
@printf("   edge %.5f  converged %s  state %s  (%d solves)\n", r.edge, r.converged, r.state, r.nsolves)
check("CONVERGED reported", r.converged, r.state)
check("edge within 0.02 dex of truth", abs(r.edge - TRUE_EDGE) < 0.02,
      @sprintf("off by %.4f", abs(r.edge - TRUE_EDGE)))
check("every recorded bracket valid", all(h -> h.valid, r.hist))

# ── 2. a seeded node that LOOKS rejected becomes accepted on re-solve ─────
# This is bug 1: the bracket is invalidated by the re-solve. The finder must
# notice and extend, not silently report an edge below the true one.
println("\n2. rejected node becomes accepted on re-solve (bug 1)")
S2 = Dict{Float64,NamedTuple}()
for x in (-2.0, -1.8, -1.6, -1.5); S2[x] = (chi = f(x) + PEN*0.2, m = x); end
S2[-1.47] = (chi = 50.0, m = -2.4)        # truly f=41.9 (accepted) but seeded as rejected
@printf("   seeded chi2 at -1.47 = %.1f (looks rejected); true f(-1.47) = %.1f (accepted)\n",
        50.0, f(-1.47))
r2 = EF.bisect_edge(S2, mock, THR; maxref = 10)
@printf("   edge %.5f  converged %s  state %s\n", r2.edge, r2.converged, r2.state)
check("did not report an edge below the stale bracket", isnan(r2.edge) || r2.edge > -1.47,
      @sprintf("edge %.4f", r2.edge))
check("edge within 0.02 dex of truth, or bracket reported lost",
      (!isnan(r2.edge) && abs(r2.edge - TRUE_EDGE) < 0.02) || occursin("bracket lost", r2.state),
      r2.state)

# ── 3. bracket bookkeeping ────────────────────────────────────────────────
println("\n3. bracket states")
a, b, st = EF.bracket(Dict(-2.0 => (chi = 60.0, m = -2.0)), THR)
check("no accepted node detected", st == "no accepted node", st)
a, b, st = EF.bracket(Dict(-2.0 => (chi = 26.0, m = -2.0), -1.9 => (chi = 27.0, m = -1.9)), THR)
check("bracket lost detected when nothing is rejected", st == "bracket lost", st)
a, b, st = EF.bracket(Dict(-2.0 => (chi = 26.0, m = -2.0), -1.4 => (chi = 60.0, m = -1.4)), THR)
check("valid bracket found", st == "valid" && a == -2.0 && b == -1.4, "$st a=$a b=$b")

# ── 4. refinement must not be able to raise the stored chi2 ───────────────
println("\n4. the stored minimum never worsens")
S4 = Dict{Float64,NamedTuple}()
for x in -2.4:0.2:-1.2; S4[x] = (chi = f(x) + PEN*0.2, m = x); end
before = Dict(k => v.chi for (k,v) in S4)
EF.bisect_edge(S4, mock, THR; maxref = 6)
check("no stored node's chi2 increased",
      all(k -> S4[k].chi <= before[k] + 1e-9, collect(keys(before))))

# ── 5. BOTH CROSSINGS: the lower edge must not be the upper one ───────────
# The regression of 19 Nov 2026. `bracket` finds an upper crossing only, so a
# caller that passes the unreflected solved set gets the upper edge back and
# reports it as the lower one. On 3S the two came out 1.7957 S and 1.7963 S,
# which is the same edge twice. A V-shaped chi2 has both crossings, so the two
# must land a known distance apart.
println("\n5. a V-shaped chi2 has two crossings, and they are different")
const XC = -1.6
const HALF = sqrt((THR - 25.0)/1000.0)         # f2 crosses at XC +/- HALF
f2(x) = 25.0 + 1000.0 * (x - XC)^2
mock2(x, from) = (chi = f2(x) + PEN * abs(x - from), m = x)

S5 = Dict{Float64,NamedTuple}()
for x in -2.2:0.1:-1.0; S5[x] = (chi = f2(x) + PEN*0.1, m = x); end
ru = EF.bisect_edge(S5, mock2, THR; maxref = 10, side = :upper)
rl = EF.bisect_edge(S5, mock2, THR; maxref = 10, side = :lower)
@printf("   upper %.5f (truth %.5f)   lower %.5f (truth %.5f)\n",
        ru.edge, XC + HALF, rl.edge, XC - HALF)
check("upper edge near the upper crossing", abs(ru.edge - (XC + HALF)) < 0.02,
      "got $(round(ru.edge,digits=5))")
check("lower edge near the LOWER crossing", abs(rl.edge - (XC - HALF)) < 0.02,
      "got $(round(rl.edge,digits=5))")
check("the two edges are not the same edge", abs(ru.edge - rl.edge) > 0.5*2HALF,
      "separation $(round(abs(ru.edge - rl.edge),digits=5)), expected $(round(2HALF,digits=5))")
check("lower edge lies below the upper one", rl.edge < ru.edge)
check("reflected solves came back with the caller's sign",
      all(k -> k > -3.0 && k < 0.0, collect(keys(S5))),
      "keys span $(round(minimum(keys(S5)),digits=3)) to $(round(maximum(keys(S5)),digits=3))")


println()
if nfail == 0
    println("ALL EDGE FINDER TESTS PASSED")
else
    @printf("%d TEST(S) FAILED\n", nfail); exit(1)
end
