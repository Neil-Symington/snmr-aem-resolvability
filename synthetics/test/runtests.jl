# runtests.jl -- the checks that need no install beyond this project.
#
#   julia --project=. test/runtests.jl
#
# These were the generator repository's own tests before the two repositories
# were consolidated. They are run by the suite at the repository root; each
# guards a bug that actually happened:
#
#   config keys     an unknown key being ignored, so a typo changed nothing
#   determinism     the same config giving the same NUMBERS (not the same bytes)
using Test, Printf
const ROOT = normpath(joinpath(@__DIR__, ".."))

@testset "snmr-aem-synthetics" begin

    # The noise-seed guard that used to live here is gone, and deliberately.
    # It existed because aem_setup.jl existed TWICE, in this repository and in
    # the workflow one, and the fix landed in only one of them. Consolidation
    # left a single copy at ../src/aem_setup.jl, which test/test_noise_seed.jl
    # at the repository root scans along with every other call site. Two guards
    # over one file is how the first one rotted.

    @testset "config keys" begin
        include(joinpath(ROOT, "src", "SNMRAEMSynthetics.jl"))
        S = SNMRAEMSynthetics

        # every shipped config must load; template.toml claims to show every
        # option, so it is the one most likely to drift from the allow-list
        for f in readdir(joinpath(ROOT, "config"); join = true)
            endswith(f, ".toml") || continue
            @test (S.load_config(f); true)
        end

        write_tmp(body) = (p = tempname() * ".toml"; write(p, body); p)
        base = read(joinpath(ROOT, "config", "minimal.toml"), String)

        # a typo must be refused, not ignored
        @test_throws Exception S.load_config(write_tmp(base *
            "\n[watertable]\ndepth = 8.0\nfring = 0.8\n"))
        @test_throws Exception S.load_config(write_tmp(base *
            "\n[variability]\nlog_sd_phii = 0.1\n"))
        # a key valid elsewhere, in the wrong section
        @test_throws Exception S.load_config(write_tmp(base * "\n[grid]\nphi = 0.3\n"))

        # but a config carrying the forward-model sections must still load:
        # scripts/forward.jl reads [aem] and [smr] off the stored config, and
        # refusing them broke its own inputs
        @test (S.load_config(write_tmp(base *
            "\n[aem]\nnoise_scale = 1.0\nnoise_frac = 0.03\nrseed = 123\n" *
            "[smr]\nloop_side = 88.6\nBe_nT = 58000.0\nnoise_nV = 20.0\n")); true)

        # [unit] instead of [[unit]] parses as one table, and iterating it used
        # to throw MethodError from inside the checker
        @test_throws Exception S.load_config(write_tmp(
            "name = \"t\"\n[unit]\ntop = 0.0\nphi = 0.3\n"))
    end

    @testset "determinism, by value" begin
        # By VALUE, not by hash. The .jld2 embeds the Julia version string and a
        # serialised Dict, so its SHA-256 differs between machines and between
        # Julia patch versions while every number is identical -- confirmed on
        # three machines. Comparing hashes would make this test fail for the
        # wrong reason.
        include(joinpath(ROOT, "src", "SNMRAEMSynthetics.jl"))
        S = SNMRAEMSynthetics
        cfg = joinpath(ROOT, "config", "minimal.toml")
        a = S.build_column(S.load_config(cfg)...)
        b = S.build_column(S.load_config(cfg)...)
        for f in (:σ_bulk, :θ_total, :clay, :σw, :S, :φ)
            x, y = getfield(a, f), getfield(b, f)
            @test maximum(abs.(x .- y)) == 0.0
        end
    end

end
