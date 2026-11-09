# test_paper_3P.jl -- the one paper column that ships must still rebuild.
#
# config/paper_3P.toml reproduces the paper's stored truth_murray3P.jld2 to
# about 1e-9 RELATIVE. It took four recovered parameters to get there (the
# fringe, the Simandoux form, the per-unit surface conductivity, and the logit
# jitter), and any of them drifting would turn that into a visible difference.
#
# THE TOLERANCE IS 1e-6 RELATIVE, AND NOT TIGHTER, FOR A MEASURED REASON. The
# first version of this file asserted 1e-10 absolute, after the author measured
# 1.1e-16 on one machine and called the reproduction exact. A reviewer ran it on
# theirs and got 6 failures out of 89: the true agreement is 1.7e-9 relative,
# and 1.1e-16 was a property of one LAPACK and thread configuration.
#
# The cause is in gp_draw: it Cholesky-factorises a squared-exponential kernel
# with a 1e-8 nugget, and cond(K) is about 1e9. Changing only
# OPENBLAS_NUM_THREADS from 4 to 1 moves phi's difference from 1.44e-10 to
# 1.75e-10. So the rebuilt field is machine- and thread-dependent at the 1e-10
# level, and any tolerance at that level tests the hardware rather than the
# config. 1e-6 is four orders above the observed spread and three below the
# 0.35 standard deviation the jitter actually applies, so it still fails if a
# recovered parameter is wrong.
#
# This is the same lesson as the truth-file SHA in round 5: compare values with
# a tolerance that survives a different machine, or the check belongs to one.
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
            # relative, not absolute: sigma_w spans 0.018 to 1.6 S/m and
            # sigma_bulk more, so one absolute bound cannot serve both.
            for f in (:S, :φ, :clay, :σw, :θ_total, :θ_SMR, :σ_bulk, :z)
                a, b = getfield(col, f), t[String(f)]
                scale = max(maximum(abs.(b)), eps())
                @test maximum(abs.(a .- b)) / scale < 1e-6
            end
            # S, z and the unit boundaries are not random-field derived, so they
            # ARE bit-identical on every machine -- asserted exactly, because a
            # drift there means a parameter changed rather than a rounding path.
            @test col.S == t["S"]
            @test col.z == t["z"]
            # and the four recovered parameters, so a drift says which one
            @test cfg["variability"]["jitter"] == "logit"
            @test cfg["watertable"]["depth"] == 6.0
            @test cfg["watertable"]["fringe"] == 0.8
            @test cfg["petrophysics"]["p_clay"] == 3.0
            @test [u.σ_clay for u in units] == [0.15, 0.15, 1.80, 0.30, 0.10]
        end
    end
end
