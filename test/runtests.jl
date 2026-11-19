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
# Each names the files it needs and skips, loudly, when they are absent.
using Test, Printf
const ROOT = normpath(joinpath(@__DIR__, ".."))

# The field workflow's control-file reader, loaded at TOP LEVEL: `using` is a
# top-level-only expression, and putting it inside the @testset made the whole
# block contribute zero tests without erroring -- which is exactly the kind of
# silent pass this suite exists to prevent, and it nearly shipped.
# Absolute, and captured here: synthetics/test/runtests.jl defines its OWN ROOT
# (as synthetics/), which clobbers this one, so a later joinpath(ROOT, ...)
# resolved to synthetics/fieldrun/... and the include failed.
const FIELDRUN_DIR = joinpath(ROOT, "fieldrun")
const FIELDCONFIG_TEST = let f = joinpath(FIELDRUN_DIR, "test", "test_fieldconfig.jl")
    isfile(f) ? f : ""
end
const UDFSYSTEM_TEST = let f = joinpath(FIELDRUN_DIR, "test", "test_udf_system.jl")
    isfile(f) ? f : ""
end
if FIELDCONFIG_TEST != ""
    include(joinpath(FIELDRUN_DIR, "common", "fieldconfig.jl"))
    using .FieldConfig
    # udf_system.jl parses the .des at load time and now cross-checks it there,
    # so loading it IS the test of the cross-check; the testset then asserts the
    # parsed values. Top level, because `using` is top-level only -- a lesson
    # from wiring the fieldconfig test inside a @testset and getting zero tests
    # with no error.
    FieldConfig.load_config(joinpath(FIELDRUN_DIR, "config", "udf_example.toml"))
    include(joinpath(FIELDRUN_DIR, "stages", "udf_system.jl"))
    using .UDFSystem
end

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
        # The seed test is behavioural now, so it needs the operator builders.
        # These includes are what the integration tests use too.
        include(joinpath(ROOT, "src", "aem_setup.jl"))
        include(joinpath(ROOT, "src", "aem_jac.jl"))
        include(joinpath(@__DIR__, "test_noise_seed.jl"))
    end

    # Stage 0 on a committed config. Nothing ran Stage 0 -- not CI, not this
    # suite -- which is how a seal that threw KeyError on every committed config
    # shipped and killed every real run for a round.
    @testset "Stage 0 seal on a committed config" begin
        include(joinpath(@__DIR__, "test_stage0_seed.jl"))
    end

    # The truth generator's own suite, from the repository it used to be.
    # Running it from here is the point of consolidating: one install, one
    # command, and the generator can no longer drift from the stages that
    # consume its output.
    include(joinpath(ROOT, "synthetics", "test", "runtests.jl"))

    # The field workflow's control-file reader. Its README claimed its failure
    # paths were "each negative-tested" while nothing was committed; this is
    # that claim made true. No stage outputs or field data needed.
    FIELDCONFIG_TEST == "" ||
        @testset "fieldrun control file" begin include(FIELDCONFIG_TEST) end
    UDFSYSTEM_TEST == "" ||
        @testset "fieldrun SkyTEM system" begin include(UDFSYSTEM_TEST) end

    # Each integration test names the FILES it needs, not a directory. The gate
    # used to be `isdir(ROOT/outputs)`, which is not what any of them read: as
    # soon as a pipeline run created outputs/, two of them fired and died on a
    # missing data/truth_murray3.jld2. Two guard themselves internally and two
    # do not, and the directory check hid the difference.
    for (name, script, needs) in (
            ("Stage 0 roundtrip", "test_stage0_roundtrip.jl", String[]),
            ("gate solver",       "test_gate_solver.jl",      ["data/truth_murray3.jld2"]),
            ("conductance constraint", "test_cond_constraint.jl", ["data/truth_murray3.jld2"]),
            ("Stage 3 integral",  "test_stage3_integral.jl",  String[]))
        f = joinpath(@__DIR__, script)
        missing_inputs = [n for n in needs if !isfile(joinpath(ROOT, n))]
        if !isfile(f)
            @printf("  SKIPPED %-26s script not present\n", name)
        elseif !isempty(missing_inputs)
            @printf("  SKIPPED %-26s needs %s\n", name, join(missing_inputs, ", "))
        else
            @testset "$name" begin; include(f); end
        end
    end
end
