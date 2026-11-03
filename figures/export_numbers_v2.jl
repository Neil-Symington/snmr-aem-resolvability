# export_numbers_v2.jl -- numbers.csv for the tier-P paper set.
#
# One row per quantity, each carrying the FIGURE it appears in, so every value
# drawn in a figure can be traced to a row and back to the file it came from.
# Rows with no figure are table-only values, marked "-".
using JLD2, Statistics, Printf, Random
const HERE = @__DIR__
const SYN  = normpath(joinpath(HERE, "..", "data"))
const OUT  = normpath(joinpath(HERE,"..","outputs"))
const MGL  = 6400.0
cellsof(z)=(mid=(z[1:end-1].+z[2:end])./2; (top=vcat(0.0,mid), bot=vcat(mid,z[end]+(z[end]-z[end-1])/2)))
function winw(z,a,b); g=cellsof(z); w=[max(0.0,min(g.bot[i],b)-max(g.top[i],a)) for i in eachindex(z)]
    ix=findall(>(1e-12),w); (ix=ix,w=w[ix]); end

rows = Vector{NTuple{7,String}}()
add(col,stat,val,truth,fig,src) = push!(rows,
    (col, stat, val isa AbstractString ? val : @sprintf("%.6g", val),
     truth === nothing ? "" : (truth isa AbstractString ? truth : @sprintf("%.6g", truth)),
     fig, src, "tier P"))

ENS = Dict("3P"=>"3PT","2P"=>"2PT","3S"=>"3ST","2S"=>"2ST")
SRC = Dict("3P"=>"3M","2P"=>"2M","3S"=>"3S","2S"=>"2S")

for c in ("3P","2P","3S","2S")
    T  = load(joinpath(SYN,"truth_murray$c.jld2"))
    S3 = load(joinpath(OUT,"synth_$c","stage3_synth_$c.jld2"))
    G  = load(joinpath(OUT,"synth_$c","stage4_synth_$c.jld2"))
    P4 = load(joinpath(OUT,"synth_$c","stage4post_synth_$c.jld2"))
    S5 = load(joinpath(OUT,"synth_$c","stage5_synth_$c.jld2"))
    w  = (G["window"].top, G["window"].base); m = T["M_ARCHIE"]
    zt, sbt, dzt, swt = T["z"], T["σ_bulk"], T["dz"], T["σw"]
    it = findall(k -> zt[k] > w[1] && zt[k] <= w[2], eachindex(zt))
    Ctrue = sum(sbt[it] .* dzt); sb = Ctrue/(w[2]-w[1])
    water = mean(swt[it]) * mean(T["θ_total"][it])^m
    add(c,"true pore water (mg/L)", mean(swt[it])*MGL, nothing, "2", "truth_murray$c.jld2")
    add(c,"true clay fraction, window mean", mean(T["clay"][it]), nothing, "2", "truth_murray$c.jld2")
    add(c,"true theta_total, window mean", mean(T["θ_total"][it]), nothing, "2,4", "truth_murray$c.jld2")
    add(c,"clay share of sigma_bulk", (sb-water)/sb, nothing, "2", "truth_murray$c.jld2")
    add(c,"true conductance over window (S)", Ctrue, nothing, "7", "truth_murray$c.jld2")
    add(c,"Stage 3 window top (m)", w[1], nothing, "5,6", "stage3_synth_$c.jld2")
    add(c,"Stage 3 window base (m)", w[2], nothing, "5,6", "stage3_synth_$c.jld2")
    add(c,"window top p50 (m)", quantile(S3["tops"],.5), nothing, "5", "stage3_synth_$c.jld2")
    add(c,"window top p5 (m)",  quantile(S3["tops"],.05), nothing, "5", "stage3_synth_$c.jld2")
    add(c,"window top p95 (m)", quantile(S3["tops"],.95), nothing, "5", "stage3_synth_$c.jld2")
    add(c,"W p50 best window (m)", quantile(S3["W"],.5), nothing, "5", "stage3_synth_$c.jld2")
    add(c,"windows passing W>=2m (frac)", S3["frac_windows_pass"], nothing, "5", "stage3_synth_$c.jld2")
    # the intervals Stage 3 reports, replacing the single picked window
    for (nm, key) in (("primary","intervals_primary"), ("secondary","intervals_secondary"))
        IV = S3[key]
        for r in 1:size(IV,1)
            add(c, @sprintf("%s interval %d: top (m)", nm, r), IV[r,1], nothing, "5", "stage3_synth_$c.jld2")
            add(c, @sprintf("%s interval %d: base (m)", nm, r), IV[r,2], nothing, "5", "stage3_synth_$c.jld2")
            add(c, @sprintf("%s interval %d: max p_pass", nm, r), IV[r,3], nothing, "5", "stage3_synth_$c.jld2")
        end
    end
    add(c,"per-window p_pass (all windows)", "see numbers_windows.csv", nothing, "5", "stage3_synth_$c.jld2")
    # --- the clay / salinity trade-off (Figure 8b) ---------------------
    # The two curves are analytic: TDS(clay) = (sigma_bulk - clay)/theta^m * MGL.
    # Three numbers per curve reproduce it exactly, so the curves are exported
    # as their parameters rather than as sampled points.
    pos5 = findall(>(0), S5["sigma_w"])
    if !isempty(pos5) && haskey(S5, "clay_product")
        thm5 = median(S5["theta_draw"][pos5] .^ S5["m_draw"][pos5])
        sbA  = median(S5["sigma_bulk"])
        it5  = findall(k -> zt[k] > w[1] && zt[k] <= w[2], eachindex(zt))
        sbT  = sum(sbt[it5] .* dzt) / (w[2] - w[1])
        add(c,"tradeoff: AEM median sigma_bulk (S/m)", sbA, sbT, "6", "stage5_synth_$c.jld2")
        add(c,"tradeoff: true sigma_bulk (S/m)", sbT, nothing, "6", "truth_murray$c.jld2")
        add(c,"tradeoff: AEM/true sigma_bulk ratio", sbA/sbT, 1.0, "6", "derived")
        add(c,"tradeoff: median theta^m of accepted draws", thm5, nothing, "6", "stage5_synth_$c.jld2")
        add(c,"tradeoff: clay floor, AEM curve (S/m)", sbA, nothing, "6", "derived")
        add(c,"tradeoff: clay floor, true curve (S/m)", sbT, nothing, "6", "derived")
        add(c,"tradeoff: truth clay term (S/m)",
            sbT - mean(swt[it5]) * mean(T["θ_total"][it5])^m, nothing, "6", "truth_murray$c.jld2")
        add(c,"tradeoff: truth TDS (mg/L)", mean(swt[it5])*MGL, nothing, "6", "truth_murray$c.jld2")
        add(c,"tradeoff: accepted draw pairs", "see numbers_tradeoff.csv", nothing, "6", "stage5_synth_$c.jld2")
    end
    add(c,"P(W>=2m), best window", S3["P_perdraw"], nothing, "-", "stage3_synth_$c.jld2")
    add(c,"lambda-hat", G["lambda_hat"], nothing, "7", "stage4_synth_$c.jld2")
    add(c,"band lower edge (S)", P4["band_C"].lo, nothing, "7b", "stage4post_synth_$c.jld2")
    add(c,"band upper edge (S)", P4["band_C"].hi, nothing, "7b", "stage4post_synth_$c.jld2")
    add(c,"admissible span (dex)", P4["band_dex"], nothing, "7b", "stage4post_synth_$c.jld2")
    add(c,"admissible nodes at lambda-hat",
        count(G["sweeps"][string(G["lambda_hat"])].good), length(G["Cg"]), "7a", "stage4_synth_$c.jld2")
    add(c,"P(TDS>1000) at lambda-hat", P4["P_wrong"], nothing, "-", "stage4post_synth_$c.jld2")
    add(c,"Delta", P4["delta"], nothing, "-", "stage4post_synth_$c.jld2")
    add(c,"data-only bound TDS p95 (mg/L)", P4["tds_data_only_p95"], mean(swt[it])*MGL, "-", "stage4post_synth_$c.jld2")
    add(c,"Stage 5 TDS p50 (mg/L)", S5["TDS_p50"], mean(swt[it])*MGL, "-", "stage5_synth_$c.jld2")
    add(c,"Stage 5 TDS p5 (mg/L)",  S5["TDS_p5"],  nothing, "-", "stage5_synth_$c.jld2")
    add(c,"Stage 5 TDS p95 (mg/L)", S5["TDS_p95"], nothing, "-", "stage5_synth_$c.jld2")
    add(c,"Stage 5 acceptance rate", S5["acceptance"], nothing, "-", "stage5_synth_$c.jld2")
    add(c,"Stage 5 P(TDS<1000)", S5["P_criterion"], nothing, "-", "stage5_synth_$c.jld2")
    # manager lines in conductance, as drawn
    E = load(joinpath(SYN,"cond_C2_typical$(ENS[c]).jld2")); ww = winw(E["zsmr"], w...)
    Θ = Matrix{Float64}(reduce(hcat,E["ens"])); Random.seed!(1)
    ms = 1.5 .+ 0.4 .* rand(4000); js = rand(1:size(Θ,2),4000)
    K = median([sum((@view Θ[ww.ix,js[t]]) .^ ms[t] .* ww.w) for t in 1:4000])
    for X in (1000.0,1500.0,3000.0)
        add(c, @sprintf("%.0f mg/L as conductance (S)", X), X/MGL*K, nothing, "7b", "cond_C2_typical$(ENS[c]).jld2")
    end
    add(c,"median int theta^m dz over window (m)", K, nothing, "7b", "cond_C2_typical$(ENS[c]).jld2")
end

# Nested intervals. Figure dropped; these are quoted in Section 5.3 only.
for (c, src, iv) in (("3P","3M",[(7.5,23.5),(0.0,40.0),(0.0,58.0),(0.0,120.0)]),
                     ("2P","2M",[(17.5,33.5),(4.0,34.0),(0.0,46.0),(0.0,120.0)]))
    T = load(joinpath(SYN,"truth_murray$c.jld2"))
    Σ = Matrix{Float64}(load(joinpath(SYN,"aem_draws_$src.jld2"),"sigma"))
    zb = collect(Float64, load(joinpath(SYN,"aem_draws_$src.jld2"),"zall_aem")); ct = vcat(0.0, zb[1:end-1])
    for (a,b) in iv
        ov = [max(0.0,min(zb[i],b)-max(ct[i],a)) for i in eachindex(zb)]; ix = findall(>(1e-9),ov)
        Ca = median([sum(Σ[ix,j] .* ov[ix]) for j in 1:size(Σ,2)])
        zt,sbt,dzt = T["z"],T["σ_bulk"],T["dz"]
        Ct = sum(sbt[findall(k -> zt[k]>a && zt[k]<=b, eachindex(zt))] .* dzt)
        add(c, @sprintf("AEM/true conductance, %.1f-%.1f m", a, b), Ca/Ct, 1.0, "-", "aem_draws_$src.jld2")
    end
end

# Sensitivity factors, test one. Replaces the S4 bound-decomposition waterfall:
# one unknown swapped at a time on the FINAL run, everything else fixed.
# Computed by synthetics/sensitivity_factors.jl; no figure.
for (lab, v, fac) in (("data only (the published bound)", 1451.0, 1.00),
                      ("clay term known", 822.0, 0.57),
                      ("porosity known", 1079.0, 0.74),
                      ("m known", 1200.0, 0.83))
    add("3P", "TDS bound, $lab (mg/L)", v, 187.0, "-", "sensitivity_factors.jl")
    lab == "data only (the published bound)" && continue
    add("3P", "sensitivity factor, $lab", fac, nothing, "-", "sensitivity_factors.jl")
end
# Stage 2 depth limit, PER COLUMN. The SNMR kernel is conductivity-sensitive,
# so each column's own conductivity sets it; stage1_typical.jld2's 70.14 m is
# the BARE lineage's (88.6 m loop) and is not these columns.
for (c, ens) in (("3P","3PT"), ("2P","2PT"), ("3S","3ST"), ("2S","2ST"))
    S = load(joinpath(SYN, "smr_data_murray$(ens).jld2"))
    K = S["kernel"]; dzk = S["kernel_dz"]; zk = S["zsmr"]; sh = S["typical_σ̂"]
    nzk, nq = size(K)
    wk = [0.5*((i>1 ? dzk[i-1] : 0.0) + (i<nzk ? dzk[i] : 0.0)) for i in 1:nzk]
    dmin = [1/sqrt(sum((wk[i]*abs(K[i,q])/sh[q])^2 for q in 1:nq)) for i in 1:nzk]
    iz = findlast(<=(0.05), dmin)
    add(c, "Stage 2 depth limit (m)", iz === nothing ? 0.0 : zk[iz], nothing, "4",
        "smr_data_murray$(ens).jld2")
    add(c, "min detectable dtheta at 40 m (v/v)", dmin[argmin(abs.(zk .- 40))], nothing, "4",
        "smr_data_murray$(ens).jld2")
end

# Figure 7 -- what the old gate rows above do not carry: the edge statuses, the
# lambda2 each bar is drawn at, and the chi2 panel's two reference values.
for c in ("3P","2P","3S","2S")
    G  = load(joinpath(OUT,"synth_$c","stage4_synth_$c.jld2"))
    P4 = load(joinpath(OUT,"synth_$c","stage4post_synth_$c.jld2"))
    SWp = G["sweeps"][string(P4["profile_lambda2"])]
    add(c, "edge_lo_status", P4["band_open"].lo ? "OPEN" : "CLOSED", nothing, "7b",
        "stage4post_synth_$c.jld2")
    add(c, "edge_hi_status", P4["band_open"].hi ? "OPEN" : "CLOSED", nothing, "7b",
        "stage4post_synth_$c.jld2")
    add(c, "profile lambda2 the bar is drawn at", P4["profile_lambda2"], nothing, "7",
        "stage4post_synth_$c.jld2")
    add(c, "min chi2 on the admissible grid", minimum(SWp.chi2[SWp.good]), nothing, "7a",
        "stage4_synth_$c.jld2")
    add(c, "chi2 threshold N + 3.84", G["ndata"] + 3.84, nothing, "7a", "stage4_synth_$c.jld2")
    add(c, "per-node chi2", "see numbers_chi2.csv", nothing, "7a", "stage4_synth_$c.jld2")
end

# Section 6 -- narrowing the Stage 0 priors with borehole logs, test two.
# Tagged "section 6" in the tier column: these rows are a truth-centred what-if
# study, not a run of the Stage 0-6 chain. See synthetics/section6_bore_value.jl.
add6(col,stat,val,truth,fig,src) = push!(rows,
    (col, stat, val isa AbstractString ? val : @sprintf("%.6g", val),
     truth === nothing ? "" : (truth isa AbstractString ? truth : @sprintf("%.6g", truth)),
     fig, src, "section 6"))
let f6 = joinpath(OUT,"section6","section6_2P.jld2")
    if isfile(f6)
        D = load(f6); nm = D["names"]; tru = D["truth_TDS"]; wt = D["window_truths"]
        add6("2P", "true window clay term (S/m)", wt.cp, nothing, "-", "section6_2P.jld2")
        add6("2P", "true window bulk conductivity (S/m)", wt.sigma_bulk, nothing, "-",
             "section6_2P.jld2")
        for k in D["order"]
            r = D["state_$k"]; st = nm[k]
            add6("2P", "$st: TDS p50 (mg/L)",  r.p50, tru, "8a", "section6_2P.jld2")
            add6("2P", "$st: TDS p5 (mg/L)",   r.p05, nothing, "8a", "section6_2P.jld2")
            add6("2P", "$st: TDS p95 (mg/L)",  r.p95, nothing, "8a", "section6_2P.jld2")
            add6("2P", "$st: P(TDS < 1000)",   r.Pfresh, nothing, "8a", "section6_2P.jld2")
            add6("2P", "$st: acceptance rate", r.accept, nothing, "-", "section6_2P.jld2")
            add6("2P", "$st: clay-term posterior p50 (S/m)", r.cp50, wt.cp, "-",
                 "section6_2P.jld2")
            add6("2P", "$st: clay-term posterior width (factor p95/p5)", r.cp95/r.cp05,
                 nothing, "-", "section6_2P.jld2")
            add6("2P", "$st: verdict", r.verdict, nothing, "8a", "section6_2P.jld2")
            r.ess < length(D["order"])*0 || add6("2P", "$st: effective sample size",
                 r.ess, nothing, "-", "section6_2P.jld2")
        end
        for k in ("2P_prior","2P_fixed","2S_prior","2S_fixed")
            x = D["xfer_$k"]; tag = split(k,"_")[1]
            add6(tag, "transfer test, $(x.label): TDS p50 (mg/L)", x.p50, x.truth, "8b",
                 "section6_2P.jld2")
            add6(tag, "transfer test, $(x.label): P(TDS < 1000)", x.Pfresh, nothing, "8b",
                 "section6_2P.jld2")
        end
    end
end

# Stage 4 salinity floor under stated ceilings (synthetics/stage4_floor.jl).
# Post-processing of finished stage outputs; no inversion involved.
let ff = joinpath(OUT,"stage4_floor","floor.jld2")
    if isfile(ff)
        D = load(ff)
        for c in ("3P","2P","3S","2S")
            haskey(D, c) || continue
            r = D[c]
            add(c, "floor: C_lo (S)", r.C_lo, nothing, "-", "stage4post_synth_$c.jld2")
            add(c, "floor: edge_lo_status", r.edge_lo_status, nothing, "-", "stage4post_synth_$c.jld2")
            add(c, "floor: theta_snmr over the window", r.theta, nothing, "-",
                "cond_C2_typical$(ENS[c]).jld2")
            add(c, "floor: clay ceiling (S/m)", r.clay_ceiling, nothing, "-", "config.toml [ceilings]")
            add(c, "floor: TDS floor (mg/L)", r.floor_tds, nothing, "-", "stage4_floor.jl")
            add(c, "floor: verdict", r.verdict, nothing, "-", "stage4_floor.jl")
        end
    end
end

open(joinpath(HERE,"numbers.csv"),"w") do f
    println(f, "column,statistic,value,truth,figure,source,tier")
    for r in rows
        println(f, join(["\"$(x)\"" for x in r], ","))
    end
end
@printf("wrote numbers.csv: %d rows, %d with a figure reference\n",
        length(rows), count(r -> r[5] != "-", rows))

# ---------------------------------------------------------------------------
# Per-window p_pass. The brief asks for window_top, window_base and p_pass in
# numbers.csv, but that is 5,709 windows x 4 columns = 22,836 rows, which would
# swamp the ~140 manuscript values numbers.csv exists to make traceable. They
# go to numbers_windows.csv instead, same provenance scheme, and numbers.csv
# carries the INTERVALS those windows form plus a pointer to the file.
open(joinpath(HERE, "numbers_windows.csv"), "w") do f
    println(f, "column,window_top_m,window_base_m,p_pass,figure,source")
    for c in ("3P","2P","3S","2S")
        S3 = load(joinpath(OUT, "synth_$c", "stage3_synth_$c.jld2"))
        t, b, pp = S3["win_tops"], S3["win_bases"], S3["win_ppass"]
        for i in eachindex(t)
            @printf(f, "\"%s\",%.2f,%.2f,%.6f,\"5\",\"stage3_synth_%s.jld2\"\n", c, t[i], b[i], pp[i], c)
        end
    end
end
# the accepted (sigma_clay, TDS) pairs, their own file for the same reason as
# the per-window rows: tens of thousands of draws would bury numbers.csv
open(joinpath(HERE, "numbers_tradeoff.csv"), "w") do f
    println(f, "column,clay_product_S_per_m,TDS_mg_per_L,figure,source")
    for c in ("3P","2P")
        S5 = load(joinpath(OUT, "synth_$c", "stage5_synth_$c.jld2"))
        haskey(S5, "clay_product") || continue
        pos = findall(>(0), S5["sigma_w"])
        for i in pos
            @printf(f, "\"%s\",%.6g,%.6g,\"8b\",\"stage5_synth_%s.jld2\"\n",
                    c, S5["clay_product"][i], S5["TDS"][i], c)
        end
    end
end
let n = 0
    for c in ("3P","2P")
        S5 = load(joinpath(OUT,"synth_$c","stage5_synth_$c.jld2"))
        haskey(S5,"clay_product") && (n += count(>(0), S5["sigma_w"]))
    end
    @printf("wrote numbers_tradeoff.csv: %d accepted draw pairs\n", n)
end
let n = 0
    for c in ("3P","2P","3S","2S")
        n += length(load(joinpath(OUT,"synth_$c","stage3_synth_$c.jld2"), "win_tops"))
    end
    @printf("wrote numbers_windows.csv: %d per-window rows\n", n)
end
