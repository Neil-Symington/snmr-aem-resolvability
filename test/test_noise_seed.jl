# test_noise_seed.jl -- the seed you pass must be the seed that draws.
#
# HiQGA's makenoisydata! calls Random.seed!(rseed) itself, so a caller that
# seeds beforehand and omits the keyword gets HiQGA's default, 123, whatever it
# thinks it set. Three tools passed rseed_noise = SEED = 20261012 into a
# function that dropped it, and every published AEM number came from the 123
# draw.
#
# FOUR VERSIONS OF THIS FILE HAVE BEEN DEFEATED, thirty-nine mutations in total
# across rounds 4 to 7. The history is the design rationale:
#
#   v1-v3  matched the source text of the call. Beaten by spelling: a keyword in
#          a comment, in a string, behind an alias, in an unscanned directory.
#   v4     went behavioural but at THREE fixed seeds, and kept a textual half to
#          supply "discovery". A reviewer then showed that `rseed_noise &
#          0x7fffffff` is the identity at all three seeds and wrong above 2^31;
#          that the discovery globbed tools/ only, so callers in figures/,
#          run_site.jl and synthetics/scripts/forward.jl were invisible; that
#          `= 123 + 1` satisfied a `=\s*123\b` prefix match; and that the
#          argument-parsing half wrongly REFUSED three valid refactors -- an old
#          call inside a #= =# block, a docstring quoting the call, and a call
#          broken across lines. It bought neither proof nor discovery.
#
# v5, this one, follows the reviewer's recommendation exactly: behavioural over a
# seed set no mask or modulus can fake, a behavioural test of the default, and
# ONE structural rule that parses nothing -- which files may call the function.
using Random

# Seeds chosen so that `& 0x7fffffff`, `% 2^31-1`, `% 1000` and a special case on
# any small set all change at least one of them. The last three straddle 2^31
# and 2^32, where a 31-bit mask starts lying.
const SEED_SET = (0, 1, 123, 2^31 - 1, 2^31, 2^32 + 5, 2^40 + 17,
                  999_983, 20_261_012, typemax(Int) - 3)

let
    zb = collect(0.0:5.0:150.0)
    σ  = fill(0.05, length(zb)); σ[zb .>= 40.0] .= 0.5
    halt(aem) = vcat(vec(aem.dlow), vec(aem.dhigh))

    # A reference draw at seed s, made by calling makenoisydata! directly. If the
    # builder's data differ from this, the builder is not drawing at s --
    # whatever it says it is doing.
    function reference(op, s)
        r = deepcopy(op.aem)
        transD_GP.SkyTEM1DInversion.makenoisydata!(r,
            log10.(op.ρ_aem_truth[2:end]); rseed = s, showplot = false,
            noisefrac = 0.03, σ_halt_low = LM_noise, σ_halt_high = HM_noise,
            units = 1e-12)
        halt(r)
    end

    @testset "build_skytem_operator honours rseed_noise" begin
        for s in SEED_SET
            op = build_skytem_operator(σ, zb; rseed_noise = s)
            @test halt(op.aem) == reference(op, s)
        end
        # distinct seeds must give distinct data, or the keyword is inert
        d = [halt(build_skytem_operator(σ, zb; rseed_noise = s).aem) for s in SEED_SET]
        @test length(unique(d)) == length(SEED_SET)
    end

    @testset "build_aem_jac honours its seed" begin
        T = Dict{String,Any}("z" => zb, "σ_bulk" => σ)
        for s in SEED_SET
            J = build_aem_jac(T, s)
            r = deepcopy(J.aem)
            transD_GP.SkyTEM1DInversion.makenoisydata!(r,
                log10.((1.0 ./ J.σt)[2:end]); rseed = s, showplot = false,
                noisefrac = 0.03, σ_halt_low = LM_noise, σ_halt_high = HM_noise,
                units = 1e-12)
            @test halt(J.aem) == halt(r)
        end
    end

    @testset "the default is the published realisation" begin
        # Behavioural, not a regex on the signature. The default was 47 -- a
        # realisation no run ever took -- and a reviewer changed it to 5, and
        # later to `123 + 1`, without the textual check noticing either.
        @test halt(build_skytem_operator(σ, zb).aem) ==
              halt(build_skytem_operator(σ, zb; rseed_noise = 123).aem)
    end
end

# --- the one structural rule ------------------------------------------------
# Not argument parsing. Just: which files may CALL makenoisydata!. Everything
# else about a call is covered behaviourally above, and a new caller anywhere in
# the repository is the one thing behaviour cannot see -- a file nothing
# exercises has no behaviour to test.
@testset "makenoisydata! is called from exactly two files" begin
    ALLOWED = Set(["src/aem_setup.jl", "src/aem_jac.jl"])
    # Comments AND docstrings. A comment mentioning the function is fine and a
    # reviewer checked it; a DOCSTRING mentioning it is a string, not a comment,
    # so it would otherwise trip the rule -- a false positive on documentation,
    # which is the direction that breaks working code.
    strip_noncode(s) = begin
        s = replace(s, r"#=(?s:.*?)=#" => "")               # block comments
        s = replace(s, Regex("\"\"\"(?s:.*?)\"\"\"") => "")   # docstrings
        join([replace(l, r"#.*$" => "") for l in split(s, '\n')], '\n')
    end
    callers = String[]
    for (root, _, files) in walkdir(ROOT)
        occursin(r"(^|/)(\.git|outputs|data|output)(/|$)", relpath(root, ROOT)) && continue
        for f in files
            endswith(f, ".jl") || continue
            rel = replace(relpath(joinpath(root, f), ROOT), '\\' => '/')
            startswith(rel, "test/") && continue          # this file, and runtests
            # ANY mention, not just `makenoisydata!(`. A reviewer aliased it --
            # `const mnd! = ...makenoisydata!` then `mnd!(aem, x; ...)` -- and the
            # aliasing file contains the name without a paren, so a call-shaped
            # match never saw it. The name may only appear in these two files.
            occursin("makenoisydata", strip_noncode(read(joinpath(root, f), String))) &&
                push!(callers, rel)
        end
    end
    @test Set(callers) == ALLOWED
    isempty(symdiff(Set(callers), ALLOWED)) ||
        @info "makenoisydata! callers changed" found = sort(callers) allowed = sort(collect(ALLOWED))
end
