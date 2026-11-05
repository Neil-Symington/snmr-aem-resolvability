# test_noise_seed.jl -- guards the inert-seed bug.
#
# HiQGA's makenoisydata! calls Random.seed!(rseed) ITSELF, immediately before
# drawing (SkyTEM1DInversion.jl). A caller that seeds beforehand and omits the
# keyword therefore gets HiQGA's default realisation, 123, no matter what seed
# it thinks it set. That is what happened here: every tool passed
# rseed_noise = SEED = 20261012 into a function that dropped it, so every
# published AEM number came from the 123 draw. The two realisations differ by
# up to 9% per gate, so the bug was not cosmetic -- and it survived because a
# second copy of aem_setup.jl in the generator repo was fixed and this one was
# not.
#
# The check is textual on purpose. The defect is a dropped keyword argument in
# one of two copies of a file, which no numerical assertion can localise and
# which needs neither HiQGA nor a truth file to detect. It proves the call site
# is still correct, not that the physics is.
let src = read(joinpath(ROOT, "src", "aem_setup.jl"), String)
    call = match(r"makenoisydata!\((?s:.*?)\)\n", src)
    @test call !== nothing
    body = call === nothing ? "" : call.match

    # the keyword must reach the noise draw
    @test occursin(r"rseed\s*=\s*rseed_noise", body)

    # and the seed that does nothing must not have come back
    @test !occursin(r"Random\.seed!\(\s*rseed_noise\s*\)", src)

    # the tools pin the published realisation rather than passing their MCMC seed
    for t in ("aem_HYP.jl", "aem_draws.jl", "aem_prior_test.jl")
        s = read(joinpath(ROOT, "tools", t), String)
        @test occursin(r"const NOISE_RSEED\s*=\s*123", s)
        @test !occursin(r"rseed_noise\s*=\s*SEED", s)
    end
end
