# test_refine_merge.jl -- the edge refinement merges two solved sets correctly.
#
# Stage 4 can insert nodes at a finer spacing around each provisional edge and
# merge them into the grid. Every per-node array must stay aligned with `lt`
# after the merge; a transposition here would attach one node's chi2 to another
# node's conductance and the edge would move for no visible reason.
using Test

@testset "refinement merge keeps every array aligned" begin
    lt   = collect(-2.0:0.5:0.0)          # 5 coarse nodes
    want = [-1.75, -1.25, -0.25]          # 3 refinement nodes
    chiA = [10.0, 20.0, 30.0, 40.0, 50.0] # coarse chi2, one per lt
    chiB = [17.5, 12.5, 2.5]              # refinement chi2, one per want
    goodA = [true, true, false, false, false]
    goodB = [true, false, false]

    allt = vcat(lt, want); ord = sortperm(allt)
    mrg(a, b) = vcat(a, b)[ord]
    lt2 = allt[ord]; chi2 = mrg(chiA, chiB); good2 = mrg(goodA, goodB)

    @test issorted(lt2)
    @test length(lt2) == 8
    @test lt2 == sort(vcat(lt, want))
    # each node keeps ITS OWN value
    pair = Dict(vcat(lt, want) .=> vcat(chiA, chiB))
    @test all(chi2[i] == pair[lt2[i]] for i in eachindex(lt2))
    pairg = Dict(vcat(lt, want) .=> vcat(goodA, goodB))
    @test all(good2[i] == pairg[lt2[i]] for i in eachindex(lt2))
    # and the refinement nodes really did land between the coarse ones
    @test lt2[2] == -1.75 && lt2[4] == -1.25

    # a model vector per node must travel with it, not be flattened
    msA = [fill(Float64(i), 3) for i in 1:5]
    msB = [fill(100.0 + i, 3) for i in 1:3]
    ms2 = mrg(msA, msB)
    @test length(ms2) == 8
    @test all(length(m) == 3 for m in ms2)
    pairm = Dict(vcat(lt, want) .=> vcat(msA, msB))
    @test all(ms2[i] == pairm[lt2[i]] for i in eachindex(lt2))
end
