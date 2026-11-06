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
    op = build_skytem_operator(σ, zb; rseed_noise = 123)
    ref = deepcopy(op.aem)
    transD_GP.SkyTEM1DInversion.makenoisydata!(ref,
        log10.(op.ρ_aem_truth[2:end]); rseed = 123, showplot = false,
        noisefrac = 0.03, σ_halt_low = LM_noise, σ_halt_high = HM_noise,
        units = 1e-12)
    @test vec(ref.dlow) == vec(op.aem.dlow)
    @test vec(ref.dhigh) == vec(op.aem.dhigh)

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
    jref = deepcopy(j1)
    transD_GP.SkyTEM1DInversion.makenoisydata!(jref,
        log10.((1.0 ./ J.σt)[2:end]); rseed = 123, showplot = false,
        noisefrac = 0.03, σ_halt_low = LM_noise, σ_halt_high = HM_noise,
        units = 1e-12)
    @test vec(jref.dlow) == vec(j1.dlow)
    @test vec(jref.dhigh) == vec(j1.dhigh)
end

# --- policy, which is genuinely textual -----------------------------------
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

    for t in ("aem_HYP.jl", "aem_draws.jl", "aem_prior_test.jl")
        s = strip_comments(read(joinpath(ROOT, "tools", t), String))
        @test occursin(r"const NOISE_RSEED\s*=\s*123", s)
        @test length(collect(eachmatch(r"NOISE_RSEED\s*=", s))) == 1
        for m in eachmatch(r"rseed_noise\s*=\s*([^,;\)\n]+)", s)
            @test strip(m.captures[1]) == "NOISE_RSEED"
        end
    end
end
