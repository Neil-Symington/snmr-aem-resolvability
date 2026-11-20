# test_stage_loads.jl -- every stage script must LOAD.
#
# WHY THIS EXISTS. stages/stage2_snmr.jl could not run at all from 17 to 18 Nov.
# A docstring added for `pooled_sigma` was inserted between `inputsig`'s
# docstring and `inputsig` itself, leaving one docstring documenting another.
# Julia PARSES that happily -- Meta.parseall returns a clean tree -- and so does
# Meta.lower on each top-level expression. It fails only when the file is
# actually evaluated:
#
#     ERROR: LoadError: cannot document the following expression:
#
# so neither a syntax check nor a lowering check catches it. Nothing else did
# either: tools/m4_run.sh skips Stage 2 deliberately, and
# test/test_stage2_pool.jl exercises `pooled_sigma` without loading the stage
# around it. The stage was broken for a day and the only thing that found it was
# trying to run it. This test is the cheapest thing that would have.
#
# A plain `include` is no good, because every stage ends in `main()` and
# including one would run the whole stage. So the trailing `main()` call is
# stripped and the rest is evaluated: every definition, every `using`, every
# `include`, with nothing executed. That is exactly the step the docstring bug
# fails, and it is verified below against the real regression's shape rather
# than assumed.
using Test

const STAGEDIR = joinpath(@__DIR__, "..", "stages")

"""
Evaluate `src` as if it were the file `asfile`, with a trailing bare `main()`
stripped. Returns the error message, or nothing if it loaded.

`asfile` is passed to `include_string` so `@__DIR__` inside the source resolves
the way it would on disk -- the stages build their paths from it, so getting
this wrong turns a clean load into a spurious missing-file error.
"""
function load_error(src::AbstractString, asfile::AbstractString;
                    args::Vector{String} = String[])
    src = replace(src, r"\n\s*main\(\)\s*\n?\z" => "\n")
    mod = Module(Symbol("LoadProbe_", hash(asfile * string(hash(src)))))
    # a fresh Module has no `include` of its own
    Core.eval(mod, :(include(q::AbstractString) = Base.include($mod, q)))
    # Most stages validate their arguments at TOP LEVEL and call error() with a
    # usage message when there is no --config. That is correct behaviour and
    # must not read as a load failure, so the probe is given a real config and
    # the script's `ARGS` resolves to this module-local one.
    isempty(args) || Core.eval(mod, :(const ARGS = $args))
    try
        include_string(mod, src, abspath(asfile))
        return nothing
    catch e
        e isa LoadError && (e = e.error)
        return sprint(showerror, e)
    end
end

load_error(path::AbstractString; kw...) = load_error(read(path, String), path; kw...)

# A real config, so the stages get past their own argument checks and on to the
# part of the file this test is about.
const CFG = joinpath(@__DIR__, "..", "configs", "synth_3P.toml")

# stage4_floor.jl checks at load that a PRIOR stage's output exists, so it can
# only be load-tested after a full run. Excluded by name, not by silence.
const NEEDS_PRIOR_OUTPUT = ["stage4_floor.jl"]

@testset "the check catches the shape that broke Stage 2" begin
    bad = "\"\"\"\ndocuments whatever is next\n\"\"\"\n" *
          "\"\"\"\nbut next is another docstring\n\"\"\"\n" *
          "function f(x)\n    x + 1\nend\n"
    err = load_error(bad, joinpath(STAGEDIR, "probe_bad.jl"))
    @test err !== nothing
    @test occursin("cannot document", err)

    good = "\"\"\"\ndocuments the function\n\"\"\"\n" *
           "function f(x)\n    x + 1\nend\n"
    @test load_error(good, joinpath(STAGEDIR, "probe_good.jl")) === nothing
end

@testset "every stage script loads" begin
    scripts = [joinpath(STAGEDIR, f) for f in sort(readdir(STAGEDIR)) if endswith(f, ".jl")]
    @test !isempty(scripts)
    for p in scripts
        if basename(p) in NEEDS_PRIOR_OUTPUT
            @info "skipped: checks for a prior stage's output at load" script = basename(p)
            continue
        end
        # Stage 2 refuses to load without --tag, by design: its run directory is
        # named from it and six untagged runs once landed in data/stage2/x_*.
        # The load test supplies one, because what it tests is that the file
        # EVALUATES, not that the argument check works -- that has its own test.
        a = ["--config", CFG]
        basename(p) == "stage2_snmr.jl" && append!(a, ["--tag", "loadtest"])
        err = load_error(p; args = a)
        err === nothing || @info "load failure" script = basename(p) error = first(split(err, "\n"))
        @test err === nothing
    end
end
