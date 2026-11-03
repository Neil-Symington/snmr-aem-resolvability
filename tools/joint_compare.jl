# joint_compare.jl -- the AEM-only / AEM+sounding comparison table.
#
#   julia --project=. joint_compare.jl 2P 2PJ 2S 2SJ
#
# Pairs are given as (aem_only_tag, joint_tag) ... and every row is read from
# the stage outputs, never recomputed, so the AEM-only column is literally the
# file that produced 5.2-5.5.
#
# The Stage 1 conductance row is taken over ONE FIXED window per pair -- the
# AEM-only column's representative window -- so the comparison is not
# confounded by Stage 3 choosing a different window on the joint draws. The
# joint column's own window is printed beside it.
using JLD2, Statistics, Printf
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE,".."))
const SYN  = joinpath(HERE, "..", "data")
const OUT  = joinpath(ROOT, "outputs")
const MGL  = 6400.0

function cond_over(drawsfile, a, b)
    Σ  = Matrix{Float64}(load(drawsfile, "sigma"))
    zb = collect(Float64, load(drawsfile, "zall_aem"))
    ct = vcat(0.0, zb[1:end-1])
    ov = [max(0.0, min(zb[i],b) - max(ct[i],a)) for i in eachindex(zb)]
    ix = findall(>(1e-9), ov)
    [sum(Σ[ix,j] .* ov[ix]) for j in 1:size(Σ,2)]
end

truth_cond(tag, a, b) = let T = load(joinpath(SYN,"truth_murray$tag.jld2"))
    z, dz = T["z"], T["dz"]
    sum(T["σ_bulk"][findall(i -> z[i] > a && z[i] <= b, eachindex(z))]) * dz
end

function row(tag; window = nothing)
    base = joinpath(OUT, "synth_$tag")
    G  = load(joinpath(base, "stage4_synth_$tag.jld2"))
    P4 = load(joinpath(base, "stage4post_synth_$tag.jld2"))
    S5 = load(joinpath(base, "stage5_synth_$tag.jld2"))
    cfgcol = G["column"]
    w = window === nothing ? (G["window"].top, G["window"].base) : window
    dr = joinpath(SYN, tag == cfgcol ? "aem_draws_$(cfgcol == "2P" ? "2M" : cfgcol).jld2" :
                       "aem_draws_$(tag).jld2")
    (post = P4, s5 = S5, window_own = (G["window"].top, G["window"].base),
     column = cfgcol, drawsfile = dr)
end

verdicts(bound) = join([(bound < X ? "PASS" : "UNDET") for X in (1000.0,1500.0,3000.0)], " / ")

function compare(pairs)
    for (a, j) in pairs
        ra = row(a); rj = row(j)
        w = ra.window_own                      # the fixed comparison window
        Ct = truth_cond(ra.column, w...)
        println("\n", "="^78)
        @printf("COLUMN %s   comparison window %.2f-%.2f m   true conductance %.4f S\n",
                ra.column, w..., Ct)
        println("="^78)
        @printf("  %-34s %22s %22s\n", "", "AEM only", "AEM + sounding")
        for (lab, fa, fj) in (
            ("Stage 3 window (m)",
             @sprintf("%.2f-%.2f", ra.window_own...), @sprintf("%.2f-%.2f", rj.window_own...)),)
            @printf("  %-34s %22s %22s\n", lab, fa, fj)
        end
        for (lab, src) in (("Stage 1 conductance p10-p90 (S)", :c), ("  p50 / truth", :r))
            va = cond_over(ra.drawsfile, w...); vj = cond_over(rj.drawsfile, w...)
            if src == :c
                @printf("  %-34s %22s %22s\n", lab,
                        @sprintf("%.4f-%.4f", quantile(va,.1), quantile(va,.9)),
                        @sprintf("%.4f-%.4f", quantile(vj,.1), quantile(vj,.9)))
            else
                @printf("  %-34s %22s %22s\n", lab,
                        @sprintf("%.2fx", median(va)/Ct), @sprintf("%.2fx", median(vj)/Ct))
            end
        end
        ba, bj = ra.post["band_C"], rj.post["band_C"]
        oa, oj = ra.post["band_open"], rj.post["band_open"]
        @printf("  %-34s %22s %22s\n", "Stage 4 band (S)",
                @sprintf("%.4f-%.4f", ba.lo, ba.hi), @sprintf("%.4f-%.4f", bj.lo, bj.hi))
        @printf("  %-34s %22s %22s\n", "  edges",
                (oa.lo ? "open" : "closed")*" / "*(oa.hi ? "open" : "closed"),
                (oj.lo ? "open" : "closed")*" / "*(oj.hi ? "open" : "closed"))
        Ba, Bj = ra.post["tds_data_only_p95"], rj.post["tds_data_only_p95"]
        @printf("  %-34s %22s %22s\n", "Data-only bound (mg/L)",
                @sprintf("%.0f", Ba), @sprintf("%.0f", Bj))
        @printf("  %-34s %22s %22s\n", "Stage 5 TDS p50 (p5-p95)",
                @sprintf("%.0f (%.0f-%.0f)", ra.s5["TDS_p50"], ra.s5["TDS_p5"], ra.s5["TDS_p95"]),
                @sprintf("%.0f (%.0f-%.0f)", rj.s5["TDS_p50"], rj.s5["TDS_p5"], rj.s5["TDS_p95"]))
        @printf("  %-34s %22s %22s\n", "  P(TDS < 1000)",
                @sprintf("%.4f", ra.s5["P_criterion"]), @sprintf("%.4f", rj.s5["P_criterion"]))
        @printf("  %-34s %22s %22s\n", "  acceptance",
                @sprintf("%.1f%%", 100ra.s5["acceptance"]), @sprintf("%.1f%%", 100rj.s5["acceptance"]))
        @printf("  %-34s %22s %22s\n", "Verdict at 1,000/1,500/3,000",
                verdicts(Ba), verdicts(Bj))
        if haskey(load(joinpath(SYN, "aem_draws_$(j).jld2")), "chi2byN_dc")
            d = load(joinpath(SYN, "aem_draws_$(j).jld2"))
            @printf("  %-34s %22s %22s\n", "Posterior-predictive chi2/N", "-",
                    @sprintf("AEM %.2f  DC %.2f", median(d["chi2byN_aem"]), median(d["chi2byN_dc"])))
        end
    end
end

ps = [(ARGS[i], ARGS[i+1]) for i in 1:2:length(ARGS)-1]
compare(isempty(ps) ? [("2P","2PJ"), ("2S","2SJ")] : ps)
