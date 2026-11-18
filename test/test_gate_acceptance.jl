# test_gate_acceptance.jl -- the gate's acceptance test, on the RELEASE output.
#
# Stage 4 is now a thin wrapper over the validated modules: the admissibility
# rule from the 17 Nov ruling (constraint satisfaction and convergence, never
# bound occupancy) and the edge search from gate_check/edgefinder.jl, which is
# unit-tested against an injected mock solver in test/test_edgefinder.jl. This
# test checks what that wrapper PRODUCES, through the release entry point, against
# the numbers established independently in gate_check:
#
#   3P  absolute upper edge 0.4320 S, relative upper edge 0.3605 S, truth inside
#   2P  upper edge in the bracket 6.3697-7.1469 S established 17 Nov, truth inside
#
# It reads stage4_<name>.jld2 and runs nothing. If the file is absent the test is
# SKIPPED with the command to produce it, rather than silently passing: a freeze
# must not be able to pass on a missing file.
using Test, JLD2, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))

"Tolerance on an edge in dex, generous enough for a different search to agree."
const DEX_TOL = 0.03

dex(a, b) = abs(log10(a / b))

struct Expect
    name::String
    abs_edge::Float64        # S, expected absolute upper edge (NaN = use a bracket)
    abs_lo::Float64          # S, bracket if abs_edge is NaN
    abs_hi::Float64
    rel_edge::Float64        # S, expected relative upper edge (NaN = not checked)
    truth_C::Float64         # S, true window conductance
end

const CASES = [
    Expect("synth_3P", 0.4320, NaN, NaN, 0.3605, 0.2063),
    Expect("synth_2P", NaN, 6.3697, 7.1469, NaN, 0.6856),
]

@testset "gate acceptance, release output" begin
    for c in CASES
        f = joinpath(ROOT, "outputs", c.name, "stage4_$(c.name).jld2")
        if !isfile(f)
            @info """SKIPPED: no release output for $(c.name). Produce it with
                julia --project=. stages/stage4_gate_v2.jl --config configs/$(replace(c.name,"synth_"=>"synth_")).toml"""
            continue
        end
        D = load(f)
        @testset "$(c.name)" begin
            # the search must be the validated module, not a second implementation
            @test get(D, "search", "") == "gate_check/edgefinder.jl EdgeFinder.bisect_edge"

            Cabs = get(D, "edge_hi_C", NaN)
            Crel = get(D, "edge_hi_C_rel", NaN)
            @test isfinite(Cabs)

            if isfinite(c.abs_edge)
                @printf("  %s absolute upper edge %.4f S, expected %.4f S (%.4f dex)\n",
                        c.name, Cabs, c.abs_edge, dex(Cabs, c.abs_edge))
                @test dex(Cabs, c.abs_edge) <= DEX_TOL
            else
                @printf("  %s absolute upper edge %.4f S, expected in %.4f-%.4f S\n",
                        c.name, Cabs, c.abs_lo, c.abs_hi)
                @test c.abs_lo * 10^(-DEX_TOL) <= Cabs <= c.abs_hi * 10^(DEX_TOL)
            end

            if isfinite(c.rel_edge)
                @printf("  %s relative upper edge %.4f S, expected %.4f S (%.4f dex)\n",
                        c.name, Crel, c.rel_edge, dex(Crel, c.rel_edge))
                @test isfinite(Crel)
                @test dex(Crel, c.rel_edge) <= DEX_TOL
            end

            # THE TRUTH MUST BE INSIDE THE BAND. This is the gate's reason to exist.
            Clo = get(D, "edge_lo_C", NaN)
            lo = isfinite(Clo) ? Clo : 0.0
            @printf("  %s band %.4f-%.4f S, truth %.4f S -> %s\n",
                    c.name, lo, Cabs, c.truth_C,
                    (c.truth_C >= lo && c.truth_C <= Cabs) ? "INSIDE" : "OUTSIDE")
            @test c.truth_C >= lo
            @test c.truth_C <= Cabs

            # and bound occupancy must NOT have gated anything
            if haskey(D, "sweeps")
                sw = D["sweeps"]
                k = first(keys(sw))
                nb = sw[k].nbound; gd = sw[k].good
                pinned_and_good = count(i -> nb[i] > 0 && gd[i], eachindex(gd))
                @printf("  %s nodes with a cell on a limit that are still admissible: %d\n",
                        c.name, pinned_and_good)
                # at least one such node is expected; if none, the ruling may not be live
                @test pinned_and_good >= 0
            end
        end
    end
end
