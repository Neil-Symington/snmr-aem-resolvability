# test_stage0_seed.jl -- Stage 0's seal, run on a COMMITTED config.
#
# WHY A COMMITTED CONFIG AND NOT A GENERATED ONE. The round-7 seal re-parsed the
# config file from disk to avoid trusting a local variable, and in doing so
# bypassed `extends` inheritance. All seven committed configs inherit `seed`
# from ../config.toml and carry none of their own, so Stage 0 ran the forward
# model and the noise draw and then died with
#
#     KeyError: key "seed" not found
#
# on EVERY committed config. It killed every real run for a round. Neither CI
# nor the suite ran Stage 0, so nothing caught it -- and the prototype that
# "verified" the seal wrote `seed = <n>` into a config it generated itself,
# which is precisely the case that cannot fail this way. A reviewer found both
# the bug and that gap in their own prototype.
#
# So this test does the one thing that would have caught it: it runs Stage 0 as
# a subprocess, on a config from configs/, with inheritance in play.
using Test, JLD2

const S0_ROOT = normpath(joinpath(@__DIR__, ".."))
const S0_CFG  = joinpath(S0_ROOT, "configs", "smoke_2P.toml")

# The inputs Stage 0 needs are stage products and are not committed, so this is
# a skip and not a failure when they are absent -- said out loud, because a
# check that silently compares nothing is how "FIREWALL ROUND-TRIP OK" came to
# print after skipping every column.
let need = [joinpath(S0_ROOT, "data", "truth_murray2P.jld2"),
            joinpath(S0_ROOT, "data", "cond_C2_typical2PT.jld2")]
    missing_inputs = filter(!isfile, need)
    if !isempty(missing_inputs)
        @info "SKIPPED test_stage0_seed: Stage 0 inputs absent" missing_inputs
    else
        art = joinpath(S0_ROOT, "outputs", "synth_2P", "stage0_synth_2P.jld2")
        rm(art, force = true)
        p = run(pipeline(Cmd(`$(Base.julia_cmd()) --project=$S0_ROOT
                              $(joinpath(S0_ROOT, "stages", "stage0_inputs.jl"))
                              --config $S0_CFG`), stdout = devnull, stderr = devnull),
                wait = false)
        wait(p)

        @testset "Stage 0 runs on a committed config" begin
            # the whole point: inheritance must not break the seal
            @test success(p)
            @test isfile(art)
        end

        if isfile(art)
            A = load(art)
            @testset "the seal records what the config asked for" begin
                # config.toml's [data].seed, inherited by every configs/*.toml
                @test Int(A["seed"]) == 123
                @test String(A["column"]) == "2P"
                # and the data are the draw at that seed, which is what the
                # seal inside Stage 0 asserts on every real run
                @test haskey(A, "d_low") || haskey(A, "dlow") ||
                      any(k -> occursin("low", String(k)), keys(A))
            end
        end
    end
end
