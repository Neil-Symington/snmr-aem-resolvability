# runtests.jl -- the checks that need no install beyond this project.
#
#   julia --project=. test/runtests.jl
#
# This repository had no tests at all until a review pointed out that its
# src/aem_setup.jl -- the file whose FIXED copy exposed the broken one in
# snmr-aem-resolvability -- was guarded by nothing. Each test below guards a
# bug that actually happened:
#
#   noise seed      rseed reaching makenoisydata!, which once silently did not
#   config keys     an unknown key being ignored, so a typo changed nothing
#   determinism     the same config giving the same NUMBERS (not the same bytes)
using Test, Printf
const ROOT = normpath(joinpath(@__DIR__, ".."))

@testset "snmr-aem-synthetics" begin

    @testset "noise seed reaches makenoisydata!" begin
        # Same guard as the workflow repository's, and for the same reason: a
        # caller that omits rseed gets HiQGA's own default, 123, whatever seed
        # it thinks it set. The two copies of aem_setup.jl diverged once
        # BECAUSE only one of them was tested.
        strip_comments(s) = join([replace(l, r"#.*$" => "") for l in split(s, '\n')], '\n')
        sites = 0
        for f in readdir(joinpath(ROOT, "src"); join = true)
            endswith(f, ".jl") || continue
            txt = strip_comments(read(f, String))
            @test !occursin(r"Random\.seed!\(\s*rseed_noise\s*\)", txt)
            for m in eachmatch(r"makenoisydata!\((?s:.*?)\)\s*\n", txt)
                sites += 1
                @test occursin(r"\brseed\s*=", m.match)
            end
        end
        @test sites >= 1
    end

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
