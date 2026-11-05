# test_noise_seed.jl -- guards the inert-seed bug at EVERY call site.
#
# HiQGA's makenoisydata! calls Random.seed!(rseed) ITSELF, immediately before
# drawing. A caller that seeds beforehand and omits the keyword gets HiQGA's
# default realisation, 123, whatever seed it thinks it set. That is what
# happened here: three tools passed rseed_noise = SEED = 20261012 into
# build_skytem_operator, which dropped it. Every published AEM number came
# from the 123 draw, and the two realisations differ by a median of about 3%
# and up to 9-14% per gate, so the pin to 123 is deliberate.
#
# The first version of this test read one hardcoded file and matched the call
# text including its comments. A review defeated it three ways: by leaving the
# keyword in a comment inside the call, by reassigning rseed_noise just before
# the call, and by dropping the keyword at the SECOND call site in
# src/aem_jac.jl, which the test never looked at. All three are closed below.
#
# The check is textual on purpose -- the defect is a dropped keyword argument,
# which needs neither HiQGA nor a truth file to detect -- but it now discovers
# its own call sites rather than trusting a filename. Comment stripping is
# line-wise and naive; a `#` inside a string literal would confuse it, and
# there are none in these files.
let
    srcdirs = ("src", "stages", "tools")
    strip_comments(s) = join([replace(l, r"#.*$" => "") for l in split(s, '\n')], '\n')

    # --- find every call site, rather than naming one ---------------------
    sites = Tuple{String,String}[]
    for d in srcdirs, f in readdir(joinpath(ROOT, d); join = true)
        endswith(f, ".jl") || continue
        txt = strip_comments(read(f, String))
        occursin("makenoisydata!", txt) || continue
        for m in eachmatch(r"makenoisydata!\((?s:.*?)\)\s*\n", txt)
            push!(sites, (relpath(f, ROOT), m.match))
        end
    end

    # If this drops to zero the test has stopped testing anything.
    @test length(sites) >= 2

    for (file, call) in sites
        # The keyword must reach the noise draw, in code and not in a comment --
        # AND it must carry the caller's argument. A review defeated an earlier
        # version by hard-coding `rseed = 999`, which satisfies "rseed is passed"
        # while ignoring the argument entirely: the original bug in a new form.
        # So the value is checked, not just the keyword.
        @test occursin(r"\brseed\s*=\s*(rseed_noise\b|Int\(seed\)|seed\b)", call)
    end

    # --- and the seed that does nothing must not come back ----------------
    for d in srcdirs, f in readdir(joinpath(ROOT, d); join = true)
        endswith(f, ".jl") || continue
        txt = strip_comments(read(f, String))
        @test !occursin(r"Random\.seed!\(\s*rseed_noise\s*\)", txt)
    end

    # --- nor may the argument be overwritten before it is used ------------
    setup = strip_comments(read(joinpath(ROOT, "src", "aem_setup.jl"), String))
    body = match(r"function build_skytem_operator(?s:.*?)\nend", setup)
    @test body !== nothing
    if body !== nothing
        # an assignment to rseed_noise anywhere in the body is a reassignment:
        # the keyword itself is in the signature, which this match excludes
        # only after the first newline, so allow the signature's default
        inner = replace(body.match, r"^function build_skytem_operator(?s:.*?)\)\s*\n" => "")
        @test !occursin(r"\brseed_noise\s*=", inner)
    end

    # --- the tools pin the published realisation, not their MCMC seed -----
    for t in ("aem_HYP.jl", "aem_draws.jl", "aem_prior_test.jl")
        s = strip_comments(read(joinpath(ROOT, "tools", t), String))
        @test occursin(r"const NOISE_RSEED\s*=\s*123", s)
        # Not merely "not SEED": every rseed_noise= in these tools must be
        # NOISE_RSEED. A review passed a third value and the test stayed green.
        for m in eachmatch(r"rseed_noise\s*=\s*([A-Za-z_0-9]+)", s)
            @test m.captures[1] == "NOISE_RSEED"
        end
        # and NOISE_RSEED itself must not be reassigned after it is defined
        @test length(collect(eachmatch(r"NOISE_RSEED\s*=", s))) == 1
    end
end
