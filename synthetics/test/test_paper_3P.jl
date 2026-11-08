# test_paper_3P.jl -- the one paper column that ships must still rebuild.
#
# config/paper_3P.toml reproduces the paper's stored truth_murray3P.jld2 to
# 1.1e-16 across every field. That is a property worth locking: it took four
# recovered parameters to get there (the fringe, the Simandoux form, the
# per-unit surface conductivity, and the logit jitter), and any of them drifting
# would silently turn an exact reproduction into an approximate one.
#
# The stored truth is NOT in the repository -- it is one of the unpublished
# products -- so this skips, loudly, when it is absent. Point PAPER_3P_TRUTH at
# it to run the check.
using Test, JLD2

let ref = get(ENV, "PAPER_3P_TRUTH",
              joinpath(ROOT, "..", "synthetics", "outputs", "HYP_SYNTH",
                       "truth_murray3P.jld2"))
    if !isfile(ref)
        @info "SKIPPED test_paper_3P: the stored truth is not available" ref
    else
        units, cfg = SNMRAEMSynthetics.load_config(
            joinpath(ROOT, "config", "paper_3P.toml"))
        col = SNMRAEMSynthetics.build_column(units, cfg)
        t = load(ref)
        @testset "paper_3P rebuilds the stored column" begin
            for f in (:S, :φ, :clay, :σw, :θ_total, :θ_SMR, :σ_bulk, :z)
                @test maximum(abs.(getfield(col, f) .- t[String(f)])) < 1e-10
            end
            # and the four recovered parameters, so a drift says which one
            @test cfg["variability"]["jitter"] == "logit"
            @test cfg["watertable"]["depth"] == 6.0
            @test cfg["watertable"]["fringe"] == 0.8
            @test cfg["petrophysics"]["p_clay"] == 3.0
            @test [u.σ_clay for u in units] == [0.15, 0.15, 1.80, 0.30, 0.10]
        end
    end
end
