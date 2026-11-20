# test_run_folders.jl -- a stage run twice must leave TWO results, not one.
#
# The ruling of 21 Nov 2026, after the third artefact lost to an overwrite:
# Stage 2 run A, the control nuisance files, and the 2P production gate, which a
# three-node smoke replaced because it inherited the same config and therefore
# the same output path. Each run now writes its own timestamped folder and the
# flat path becomes a pointer.
#
# Stage 0 is the stage used here because it is seconds, not hours. The property
# under test belongs to RunOut and is shared by every stage.
using Test, Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "runout.jl")); using .RunOut

@testset "a run never overwrites a run" begin
    base = mktempdir()

    d1 = RunOut.run_dir(base, "stageX", "col")
    write(joinpath(d1, "result.jld2"), "first")
    RunOut.publish_run(base, d1)

    d2 = RunOut.run_dir(base, "stageX", "col")
    write(joinpath(d2, "result.jld2"), "second")
    RunOut.publish_run(base, d2)

    runs = sort(readdir(joinpath(base, "runs")))
    @test length(runs) == 2                       # two runs, two folders
    @test d1 != d2
    @test read(joinpath(d1, "result.jld2"), String) == "first"   # the first survives
    @test read(joinpath(d2, "result.jld2"), String) == "second"

    # the flat path is a POINTER, and it points at the newest run
    link = joinpath(base, "result.jld2")
    @test islink(link)
    @test read(link, String) == "second"
    @test strip(read(joinpath(base, "LATEST"), String)) == basename(d2)

    # and a result file itself refuses to be written over
    @test_throws ErrorException RunOut.save_result(p -> write(p, "third"),
                                                   joinpath(d2, "result.jld2"))
    @test read(joinpath(d2, "result.jld2"), String) == "second"

    # same second, different folder
    many = [RunOut.run_dir(base, "stageY", "col") for _ in 1:3]
    @test length(unique(many)) == 3
end
