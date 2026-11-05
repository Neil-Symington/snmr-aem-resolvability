# runtests.jl -- the checks that do not need a finished inversion.
#
#   julia --project=. test/runtests.jl
#
# Each of these guards a bug that actually happened in this project:
#   config              inheritance, and every validation rule firing
#   noise seed          rseed reaching makenoisydata!, which once silently did not
#   synthetics/         the truth generator's own suite, included from here
#   stage0 roundtrip    the sealed artefact rebuilding the operator exactly
#   gate solver         the conductance constraint being met at the optimum
#   stage3 integral     the window integral against a direct sum
# The last three need stage outputs and are skipped, loudly, without them.
using Test, Printf
const ROOT = normpath(joinpath(@__DIR__, ".."))

@testset "snmr-aem-resolvability" begin
    @testset "config inheritance and validation" begin
        include(joinpath(ROOT, "src", "config.jl"))
        for f in readdir(joinpath(ROOT, "configs"); join = true)
            endswith(f, ".toml") || continue
            c = SiteConfig.resolve(f)
            i = SiteConfig.assertvalid(c; root = ROOT, require_files = false)
            @test i.m_prior[1] < i.m_prior[2]
            @test i.role in ("primary", "saline_control")
            # a control without its parent named has produced one retracted
            # finding in this project; the validator must refuse it
            i.role == "saline_control" && @test !isempty(i.parent)
        end
    end

    @testset "noise seed reaches makenoisydata!" begin
        include(joinpath(@__DIR__, "test_noise_seed.jl"))
    end

    # The truth generator's own suite, from the repository it used to be.
    # Running it from here is the point of consolidating: one install, one
    # command, and the generator can no longer drift from the stages that
    # consume its output.
    include(joinpath(ROOT, "synthetics", "test", "runtests.jl"))

    for (name, script, needs) in (
            ("Stage 0 roundtrip", "test_stage0_roundtrip.jl", "outputs"),
            ("gate solver",       "test_gate_solver.jl",      "outputs"),
            ("conductance constraint", "test_cond_constraint.jl", "outputs"),
            ("Stage 3 integral",  "test_stage3_integral.jl",  "outputs"))
        f = joinpath(@__DIR__, script)
        if isfile(f) && isdir(joinpath(ROOT, needs))
            @testset "$name" begin; include(f); end
        else
            @printf("  SKIPPED %-26s needs stage outputs (%s)\n", name, needs)
        end
    end
end
