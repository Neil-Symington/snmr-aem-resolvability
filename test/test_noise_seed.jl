# test_noise_seed.jl -- the noise seed must REACH the draw. Behavioural.
#
# HiQGA's makenoisydata! calls Random.seed!(rseed) itself, immediately before
# drawing, so a caller that seeds beforehand and omits the keyword gets HiQGA's
# default realisation, 123, whatever seed it thinks it set. Three tools passed
# rseed_noise = SEED = 20261012 into a function that dropped it, so every
# published AEM number came from the 123 draw. The realisations differ by a
# median of about 3% and up to 9-14% per gate.
#
# WHY THIS IS NO LONGER A TEXT CHECK. Three successive versions of this file
# matched the source of the call site, and a reviewer defeated every one of
# them -- sixteen distinct mutations in total: `rseed = seed+1` (matched
# `seed\b`), `rseed_noise += 1` (not `=`), a decoy keyword inside a nested
# argument, an aliased makenoisydata!, a new call site in a subdirectory
# readdir never visited, a new call site in a directory the test did not scan.
# Each patch closed the exact spelling used and left the class open, which is
# the signature of testing the wrong thing.
#
# So this builds the operator and compares the DATA. Any mutation that stops
# the caller's seed reaching the draw -- dropped, hardcoded, reassigned,
# aliased, or routed through a new call site in the same path -- changes the
# answer below. Spelling is irrelevant.
using Random

let
    # A minimal two-layer section. No truth file, no config, no scenario: the
    # question is only whether the seed propagates.
    zb = collect(0.0:5.0:150.0)
    σ  = fill(0.05, length(zb)); σ[zb .>= 40.0] .= 0.5

    data(op) = vcat(vec(op.aem.dlow), vec(op.aem.dhigh))

    a1 = data(build_skytem_operator(σ, zb; rseed_noise = 123))
    a2 = data(build_skytem_operator(σ, zb; rseed_noise = 123))
    b  = data(build_skytem_operator(σ, zb; rseed_noise = 20261012))

    # the same seed must give the same realisation ...
    @test a1 == a2
    # ... and a different seed a different one. If this fails, the keyword is
    # not reaching makenoisydata! and every "different seed" is the same draw.
    @test a1 != b
    # not merely different: the realisations are worth percent, not epsilon
    @test maximum(abs.((a1 .- b) ./ max.(abs.(b), eps()))) > 0.01

    # AND IT MUST BE THE REALISATION ASKED FOR, not merely a reproducible one.
    # "seeds matter" is too weak: a review mutated the call to `rseed = seed+1`,
    # which keeps different seeds giving different data while silently drawing
    # the wrong realisation. So the operator built at 123 is compared against
    # makenoisydata! invoked at 123 directly. Version-independent -- no
    # hardcoded fingerprint to rot when HiQGA changes.
    # AT MORE THAN ONE SEED. Checking identity only at 123 let a reviewer pass
    # `rseed = rseed_noise % 1000`: it is the identity at 123, and the variation
    # check ran at a different seed, so a function that agrees at 123 and
    # differs elsewhere satisfied both halves. 20261012 is included because it
    # is the seed the tools used to pass, and %1000 maps it to 12.
    for sd in (123, 20261012, 999_983)
        op = build_skytem_operator(σ, zb; rseed_noise = sd)
        ref = deepcopy(op.aem)
        transD_GP.SkyTEM1DInversion.makenoisydata!(ref,
            log10.(op.ρ_aem_truth[2:end]); rseed = sd, showplot = false,
            noisefrac = 0.03, σ_halt_low = LM_noise, σ_halt_high = HM_noise,
            units = 1e-12)
        @test vec(ref.dlow) == vec(op.aem.dlow)
        @test vec(ref.dhigh) == vec(op.aem.dhigh)
    end

    # the Jacobian operator is a SECOND path to the same draw, and a review
    # mutated it without this file noticing
    T = Dict{String,Any}("z" => zb, "σ_bulk" => σ)
    J = build_aem_jac(T, 123)
    j1 = J.aem
    j2 = build_aem_jac(T, 20261012).aem
    @test vec(j1.dlow) != vec(j2.dlow)

    # The SAME identity check as above, and it is needed: the first version of
    # this file checked only that different seeds differ, and a review's
    # `rseed = Int(seed)+1` mutation passed it -- different seeds still differ
    # while the realisation drawn is not the one asked for. Checking variation
    # is not checking identity.
    for sd in (123, 20261012, 999_983)
        Js = build_aem_jac(T, sd)
        jref = deepcopy(Js.aem)
        transD_GP.SkyTEM1DInversion.makenoisydata!(jref,
            log10.((1.0 ./ Js.σt)[2:end]); rseed = sd, showplot = false,
            noisefrac = 0.03, σ_halt_low = LM_noise, σ_halt_high = HM_noise,
            units = 1e-12)
        @test vec(jref.dlow) == vec(Js.aem.dlow)
        @test vec(jref.dhigh) == vec(Js.aem.dhigh)
    end
end

# --- policy, which is genuinely textual -----------------------------------
# NOT a complete defence, and the header should say so: a reviewer defeated the
# text half with `@eval const $(Symbol("NOISE_RSE","ED")) = 7`, which no textual
# check can see, and with a `#` inside a string literal that the naive comment
# stripper cuts at. What this half buys is DISCOVERY -- finding call sites and
# callers that no behavioural test knows to exercise -- not proof.
# Which realisation the TOOLS pin is a decision, not a behaviour: 123 is the
# draw behind every published number, so the tools must pin it and must not
# pass their MCMC seed. Recursive, because a reviewer added a call site in a
# subdirectory that a flat readdir never saw.
let
    strip_comments(s) = join([replace(l, r"#.*$" => "") for l in split(s, '\n')], '\n')
    jlfiles(dir) = let out = String[]
        for (root, _, files) in walkdir(dir), f in files
            endswith(f, ".jl") && push!(out, joinpath(root, f))
        end
        out
    end

    for d in ("src", "stages", "tools", "synthetics")
        for f in jlfiles(joinpath(ROOT, d))
            txt = strip_comments(read(f, String))
            # the inert seed must not come back anywhere
            @test !occursin(r"Random\.seed!\(\s*rseed_noise\s*\)", txt)
        end
    end

    # DISCOVER the callers rather than naming three. A reviewer added a new tool
    # and changed a fourth that the hardcoded list did not mention, and both
    # went unnoticed; and the per-file loop asserted NOTHING when a file had
    # zero rseed_noise= matches, so simply deleting the keyword passed.
    callers = [f for f in jlfiles(joinpath(ROOT, "tools"))
               if occursin("build_skytem_operator", strip_comments(read(f, String)))]
    @test !isempty(callers)      # if this empties, the check has stopped checking

    for f in callers
        s = strip_comments(read(f, String))
        name = basename(f)
        # every caller must pass the keyword, and must pass NOISE_RSEED
        ms = collect(eachmatch(r"rseed_noise\s*=\s*([^,;\)\n]+)", s))
        @test !isempty(ms)       # closes "drop the keyword and inherit the default"
        for m in ms
            @test strip(m.captures[1]) == "NOISE_RSEED"
        end
        # \b after 123, or `const NOISE_RSEED = 1230` satisfies it
        @test occursin(r"const NOISE_RSEED\s*=\s*123\b", s)
        @test length(collect(eachmatch(r"NOISE_RSEED\s*=", s))) == 1
    end

    # The builder's own default is the realisation behind every published
    # number, so a caller that omits the keyword gets the right one rather than
    # a draw nobody used. It was 47 -- a realisation no run ever took -- which
    # a reviewer changed to 5 without the suite noticing.
    setup = strip_comments(read(joinpath(ROOT, "src", "aem_setup.jl"), String))
    @test occursin(r"rseed_noise::Int\s*=\s*123\b", setup)

    # Every makenoisydata! call site anywhere must pass rseed. This is the
    # DISCOVERY the behavioural test cannot do: a new file calling it directly
    # has no behavioural assertion at all.
    nsites = 0
    for d in ("src", "stages", "tools", "synthetics")
        for f in jlfiles(joinpath(ROOT, d))
            txt = strip_comments(read(f, String))
            for m in eachmatch(r"makenoisydata!\((?s:.*?)\)\s*\n", txt)
                nsites += 1
                @test occursin(r"\brseed\s*=", m.match)
            end
        end
    end
    @test nsites >= 2
end
