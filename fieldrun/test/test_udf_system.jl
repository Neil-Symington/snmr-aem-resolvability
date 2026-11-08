# test_udf_system.jl -- the .des cross-check, and that it is not optional.
#
# check_against_config was defined and called NOWHERE for two rounds, while the
# module header and the README both said the gate table "is cross-checked". It
# now runs on load.
#
# WHAT THIS FILE DOES AND DOES NOT TEST. It asserts that the SHIPPED .des passes
# every identity. It does NOT construct a bad .des and assert that `include`
# fails -- an earlier version of this header said it did, which was a claim
# about a test that is not here, found by a reviewer reading the header against
# the code. Doing it properly needs a copy of the tree per mutation, which is
# what the reviewer's own harness does and what this suite is not set up for.
# The mutations are recorded in the round-7 and round-8 reports instead.
using Test
const FR = normpath(joinpath(@__DIR__, ".."))

@testset "the shipped .des passes every identity" begin
    n = UDFSystem.check_against_config()
    @test n.lm == 18
    @test n.hm == 23
    @test first(UDFSystem.LM_GATES) == 9 && last(UDFSystem.LM_GATES) == 26
    @test first(UDFSystem.HM_GATES) == 16 && last(UDFSystem.HM_GATES) == 38
    # the on-time limb: without it the forward is told the current had been on
    # forever, which over-predicted LM gate 26 by 52% on this floodplain
    @test UDFSystem.LM_ramp[1, 1] < 0
    @test UDFSystem.HM_ramp[1, 1] < 0
    @test maximum(UDFSystem.LM_ramp[:, 2]) ≈ 1.0
    # geometry from the config, which is DECLARED not verified: the .des does
    # not contain these values at all
    @test UDFSystem.rTx ≈ sqrt(342.0 / π)
    @test UDFSystem.rRx == 13.29
end
