# diagnose_production.jl -- Part A on the 2P PRODUCTION run, in place.
#
# Run this on the machine that holds the production chains. It changes nothing;
# it only reads. Nothing is written outside convergence_test/results/.
#
#   julia --project=. convergence_test/diagnose_production.jl [stage1_dir] [condfile]
#
# Both arguments are auto-discovered if omitted.
#
# Stage 1 production ran ONE ensemble with nchainsatone = 1, so there is a
# single T=1 sequence and a between-ensemble R-hat does not exist. This splits
# that sequence into four consecutive segments, which detects DRIFT but cannot
# detect a mode the chain never visited. That limitation is printed with the
# number, every time, because the number looks like an R-hat and is not one.
using JLD2, Statistics, Printf

const HERE = normpath(@__DIR__)
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(ROOT, "data")
const GADI_SRC = joinpath(ROOT, "src")
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "prerun_config.jl"))
include(joinpath(HERE, "prerun_lib.jl"))
using HiQGA.transD_GP
using .PrerunLib
const CFG = PrerunConfig

const WIN = (17.5, 33.5)     # test two: AQUIFER 1, trimmed
const PKG = (4.0, 34.0)      # test two: clay cap + aquifer

function discover_stage1()
    length(ARGS) >= 1 && return ARGS[1]
    for c in ["aem2P", "aem2PT", "aem2"]
        d = joinpath(OUT, c)
        isdir(d) && isfile(joinpath(d, "aem_setup.jld2")) && return d
    end
    cands = filter(d -> startswith(basename(d), "aem2") &&
                        isfile(joinpath(d, "aem_setup.jld2")),
                   [joinpath(OUT, x) for x in readdir(OUT)])
    isempty(cands) ? "" : first(cands)
end
function discover_cond()
    length(ARGS) >= 2 && return ARGS[2]
    cands = filter(f -> startswith(basename(f), "cond_C2_") && endswith(f, ".jld2") &&
                        occursin("2P", basename(f)),
                   [joinpath(OUT, x) for x in readdir(OUT)])
    isempty(cands) && (cands = filter(f -> startswith(basename(f), "cond_C2_") && endswith(f, ".jld2"),
                                      [joinpath(OUT, x) for x in readdir(OUT)]))
    isempty(cands) ? "" : first(cands)
end

cell_edges(ztop, zcen) = (n = length(ztop);
    [i <= n ? ztop[i] : ztop[n] + 2*(zcen[n]-ztop[n]) for i in 1:n+1])

function conductance(s, e, a, b)
    t = 0.0
    for i in eachindex(s)
        ov = min(e[i+1], b) - max(e[i], a); ov > 0 && (t += s[i]*ov)
    end
    t
end

function chainstats(costs)
    !isfile(costs) && return nothing
    rows = [split(l) for l in eachline(costs) if !isempty(strip(l))]
    (isempty(rows) || length(rows[1]) < 10) && return nothing
    bnd = [r for r in rows if (parse(Int, r[2]) - 1) % 100 == 0]
    mn(rs, c) = (v = Float64[]; for r in rs; x = tryparse(Float64, r[c]);
                 x !== nothing && isfinite(x) && push!(v, x); end;
                 isempty(v) ? NaN : mean(v))
    isT1(r) = isapprox(parse(Float64, r[10]), 1.0; atol = 1e-12)
    hotT = maximum(parse(Float64, r[10]) for r in rows)
    isHot(r) = isapprox(parse(Float64, r[10]), hotT; atol = 1e-12)
    t1 = sort([r for r in rows if isT1(r)], by = r -> parse(Int, r[2]))
    nodes = [parse(Int, r[8]) for r in t1]
    i0 = max(1, floor(Int, 0.5*length(nodes))); np = nodes[i0:end]
    # how often does the T=1 SLOT change hands? a lower bound on swaps, because
    # two swaps inside one save_freq interval look like none.
    occ = [parse(Int, r[1]) for r in t1]
    swap = length(occ) < 2 ? NaN : count(!=(0), diff(occ)) / (length(occ)-1)
    (allT = [mn(bnd,c) for c in 3:7],
     T1   = [mn(filter(isT1,  bnd), c) for c in 3:7],
     hot  = [mn(filter(isHot, bnd), c) for c in 3:7],
     hotT = hotT, nodes = np,
     nchange = length(np) < 2 ? NaN : count(!=(0), diff(np))/(length(np)-1),
     swap = swap, nT1 = length(t1), nrows = length(rows))
end

S1 = discover_stage1(); CF = discover_cond()
println("="^78)
println("PART A on the 2P PRODUCTION run")
println("stage 1 dir : ", isempty(S1) ? "NOT FOUND" : S1)
println("cond file   : ", isempty(CF) ? "NOT FOUND" : CF)
@printf("gate        : R-hat < %.2f and ESS >= %d\n", CFG.GATE_RHAT, CFG.GATE_MIN_ESS)
println("="^78)

results = String[]

# ── Stage 1 ─────────────────────────────────────────────────────────────────
if !isempty(S1)
    S = load(joinpath(S1, "aem_setup.jld2"))
    znall = S["znall_aem"]; zall = S["zall_aem"]
    edges = cell_edges(S["zboundaries_aem"], zall)
    nmaxk = haskey(S, "nmax_k") ? S["nmax_k"] : 40
    @printf("\nStage 1: nsamples=%s nchains=%s Tmax=%s nmax=%d\n",
            string(get(S,"nsamples","?")), string(get(S,"nchains","?")),
            string(get(S,"Tmax","?")), nmaxk)
    opt = make_aem_opt(znall, "hypsynth_aem_"; nmax = nmaxk)
    opt.xall[:] .= znall
    old = pwd(); cd(S1)
    F = try
        transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
    finally; cd(old); end
    sigma = 10 .^ (-hcat([f[:] for f in F]...))
    n = size(sigma, 2)
    @printf("assembled %d post-burn-in T=1 draws (ONE ensemble)\n", n)
    q2 = [conductance(view(sigma,:,j), edges, WIN...) for j in 1:n]
    q3 = [conductance(view(sigma,:,j), edges, PKG...) for j in 1:n]
    @printf("\nQ2 window conductance  (S): p05 %.4f  p50 %.4f  p95 %.4f\n",
            quantile(q2,0.05), median(q2), quantile(q2,0.95))
    @printf("Q3 package conductance (S): p05 %.4f  p50 %.4f  p95 %.4f\n",
            quantile(q3,0.05), median(q3), quantile(q3,0.95))
    seg(v, k=4) = (m = length(v) ÷ k; [collect(Float64, v[(i-1)*m+1:i*m]) for i in 1:k])
    println("\nSEGMENT R-hat -- NOT a between-ensemble R-hat.")
    println("  One chain split into 4 consecutive segments. Detects drift only;")
    println("  it CANNOT detect a mode the chain never visited, so it is a")
    println("  lower bound on the real disagreement.")
    @printf("%-28s %8s %10s %10s\n", "quantity", "segR-hat", "ESS bulk", "ESS tail")
    for (nm, v) in [("Q2 window conductance", q2), ("Q3 package conductance", q3)]
        r = rhat_and_ess(seg(v))
        @printf("%-28s %8.3f %10.0f %10.0f\n", nm, r.rhat, r.ess_bulk, r.ess_tail)
        push!(results, @sprintf("1,2P,production,%s,%.4f,%.1f,%.1f", nm, r.rhat, r.ess_bulk, r.ess_tail))
    end
    ls = log10.(sigma)
    println("\nper-depth log10 conductivity (segment R-hat):")
    @printf("%-10s %10s %10s\n", "depth", "segR-hat", "ESS bulk")
    for zq in CFG.GATE_DEPTHS
        zq > maximum(zall) && continue
        i = argmin(abs.(zall .- zq)); r = rhat_and_ess(seg(vec(ls[i,:])))
        @printf("%-10s %10.3f %10.0f\n", @sprintf("%.0f m", zall[i]), r.rhat, r.ess_bulk)
    end
    st = chainstats(joinpath(S1, "misfits_hypsynth_aem_s.bin"))
    if st !== nothing
        @printf("\ndimension mixing at T=1: median k %.1f, range %d-%d, k changes %.2f%% of steps\n",
                median(st.nodes), minimum(st.nodes), maximum(st.nodes), 100*st.nchange)
        @printf("T=1 slot changes hands %.2f%% of logged steps (lower bound on swaps)\n", 100*st.swap)
        @printf("\n%-16s %8s %8s %8s %8s %8s\n","acceptance %","birth","death","posn","prop","dc")
        @printf("%-16s %8.2f %8.2f %8.2f %8.2f %8.2f\n","all T", st.allT...)
        @printf("%-16s %8.2f %8.2f %8.2f %8.2f %8.2f\n","T=1", st.T1...)
        @printf("%-16s %8.2f %8.2f %8.2f %8.2f %8.2f\n",
                @sprintf("T=%.2f (hot)", st.hotT), st.hot...)
    end
else
    println("\nStage 1 production directory not found -- pass it as argument 1.")
end

# ── Stage 2 ─────────────────────────────────────────────────────────────────
if !isempty(CF) && isfile(CF)
    D = load(CF); E = D["ens"]; zsmr = D["zsmr"]
    @printf("\n\nStage 2: %d ensembles, %s draws each, nsamples=%s\n",
            length(E), size(E[1],2), string(get(D,"nsamples","?")))
    idx = findall(z -> WIN[1] <= z <= WIN[2], zsmr)
    Q1 = [vec(mean(view(M, idx, :), dims=1)) for M in E]
    @printf("\n-- Q1 mean water content over %.1f-%.1f m (v/v) --\n", WIN...)
    @printf("%-10s %10s %10s %10s\n","ensemble","p05","p50","p95")
    for (i,v) in enumerate(Q1)
        @printf("%-10d %10.4f %10.4f %10.4f\n", i, quantile(v,0.05), median(v), quantile(v,0.95))
    end
    r = rhat_and_ess(Q1)
    meds = [median(v) for v in Q1]
    ok = isfinite(r.rhat) && r.rhat < CFG.GATE_RHAT && r.ess_bulk >= CFG.GATE_MIN_ESS
    @printf("\nQ1: R-hat %.3f  ESS bulk %.0f  ESS tail %.0f  median spread %.2fx  %s\n",
            r.rhat, r.ess_bulk, r.ess_tail, maximum(meds)/minimum(meds), ok ? "PASS" : "FAIL")
    push!(results, @sprintf("2,2P,production,Q1 window water content,%.4f,%.1f,%.1f",
                            r.rhat, r.ess_bulk, r.ess_tail))
    println("\nper-depth water content (true between-ensemble R-hat):")
    @printf("%-10s %10s %10s %8s\n","depth","R-hat","ESS bulk","pass")
    for zq in CFG.GATE_DEPTHS
        zq > maximum(zsmr) && continue
        i = argmin(abs.(zsmr .- zq)); rr = rhat_and_ess([vec(M[i,:]) for M in E])
        pp = isfinite(rr.rhat) && rr.rhat < CFG.GATE_RHAT && rr.ess_bulk >= CFG.GATE_MIN_ESS
        @printf("%-10s %10.3f %10.0f %8s\n", @sprintf("%.0f m", zsmr[i]),
                rr.rhat, rr.ess_bulk, pp ? "PASS" : "FAIL")
    end
    allr = Float64[]; alle = Float64[]
    for i in eachindex(zsmr)
        rr = rhat_and_ess([vec(M[i,:]) for M in E])
        isfinite(rr.rhat) && push!(allr, rr.rhat); isfinite(rr.ess_bulk) && push!(alle, rr.ess_bulk)
    end
    !isempty(allr) && @printf("all %d depths: median R-hat %.3f, max %.3f; median ESS %.0f\n",
                              length(allr), median(allr), maximum(allr), median(alle))
    # node counts per ensemble, from the Stage 2 chain files
    rundir = joinpath(OUT, "smr_runs")
    if isdir(rundir)
        pats = filter(f -> startswith(f,"misfits_C2_") && endswith(f,"_s.bin") &&
                           occursin("2P", f), readdir(rundir))
        if !isempty(pats)
            println("\nStage 2 per-ensemble dimension mixing and acceptance:")
            @printf("%-34s %8s %8s %10s %8s %8s\n","file","med k","k chg%","swap% T1","birth@1","death@1")
            for f in sort(pats)
                st = chainstats(joinpath(rundir, f)); st === nothing && continue
                isempty(st.nodes) && continue
                @printf("%-34s %8.1f %8.2f %10.2f %8.2f %8.2f\n", f, median(st.nodes),
                        100*st.nchange, 100*st.swap, st.T1[1], st.T1[2])
            end
        end
    end
else
    println("\n\nStage 2 cond file not found -- pass it as argument 2.")
end

mkpath(joinpath(HERE, "results"))
open(joinpath(HERE, "results", "production_2P.csv"), "w") do io
    println(io, "stage,tag,setting,quantity,rhat_or_segrhat,ess_bulk,ess_tail")
    for r in results; println(io, r); end
end
println("\nwrote ", joinpath(HERE, "results", "production_2P.csv"))
