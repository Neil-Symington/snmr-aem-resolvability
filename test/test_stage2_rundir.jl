# test_stage2_rundir.jl -- the Stage 2 run folder must carry --label.
#
# A spliced-kernel run reaches Stage 2's DEFAULT kernel path through --label:
# the kernel file is data/aem<TTAG><LABEL>/aem_posterior.jld2. Without the
# label in the folder name such a run lands in <tag>_pooled_<stamp>,
# indistinguishable from a genuine pooled run, which is how four duplicate
# runs already went unnoticed once (22 Nov). With no --label the name must be
# exactly what it was.
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))

"Evaluate stage2_snmr.jl with `args` and return its RUNDIR, without running it."
function rundir_for(args::Vector{String})
    f = joinpath(ROOT, "stages", "stage2_snmr.jl")
    src = replace(read(f, String), r"\n\s*main\(\)\s*\n?\z" => "\n")  # define, do not run
    mod = Module(Symbol("RunDirProbe_", hash(string(args))))
    # a fresh Module has no `include` of its own, and the stage uses it
    Core.eval(mod, :(include(q::AbstractString) = Base.include($mod, q)))
    Core.eval(mod, :(const ARGS = $args))
    include_string(mod, src, abspath(f))     # abspath so @__DIR__ resolves
    Core.eval(mod, :(RUNDIR))
end

@testset "Stage 2 run folder carries the label" begin
    plain = rundir_for(["--tag", "2Pv2T"])
    @test occursin("2Pv2T_pooled_", basename(plain))
    # the stamp is all that follows; no stray separator where a label would go
    @test match(r"^2Pv2T_pooled_\d{8}-\d{6}$", basename(plain)) !== nothing

    lab = rundir_for(["--tag", "2Pv2T", "--label", "_splicewindow"])
    @test occursin("2Pv2T_pooled_splicewindow_", basename(lab))
    @test match(r"^2Pv2T_pooled_splicewindow_\d{8}-\d{6}$", basename(lab)) !== nothing

    # and the two are distinguishable, which is the whole point
    @test basename(plain) != basename(lab)

    # --truekernel still names the kind, with and without a label
    tk = rundir_for(["--tag", "2Pv2T", "--truekernel"])
    @test occursin("2Pv2T_truekernel_", basename(tk))
    tkl = rundir_for(["--tag", "2Pv2T", "--truekernel", "--label", "_x"])
    @test occursin("2Pv2T_truekernel_x_", basename(tkl))
end
