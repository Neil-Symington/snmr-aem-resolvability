# diagnose_stage1.jl -- convergence diagnostics for Stage 1 (AEM), across
# independent ensembles. Part A/C of the sampler convergence brief.
#
#   julia --project=. convergence_test/diagnose_stage1.jl <TAG> <SETTING> [NENS]
#
# Reads data/aem<TAG>_<SETTING>e<N>/ for N = 1..NENS, re-assembles ALL T=1
# post-burn-in draws from the chain files (aem_posterior.jld2 keeps only 30, far
# too few for an ESS), and reports, per ensemble and across ensembles:
#
#   Q2  conductance of the representative window      (S)
#   Q3  conductance of the package                    (S)
#   per-depth log10 sigma at the gate depths
#   node count at T=1: median, and the FRACTION OF CONSECUTIVE T=1 SAMPLES AT
#       WHICH THE NODE COUNT CHANGES -- the exact measure of dimension mixing
#   windowed acceptance rates by move type
#
# Conductance, not resistivity; factors, not decades.
using JLD2, Statistics, Printf, Random

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

const TAG     = length(ARGS) >= 1 ? ARGS[1] : "2P"
const SETTING = length(ARGS) >= 2 ? ARGS[2] : "b1"
const NENS    = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 4

# Windows from the brief. Test two's window is AQUIFER 1 (18-34 m) trimmed; its
# package is clay cap + aquifer (4-34 m), which is what AEM actually pins.
# Test one is UNCONFINED: there is no clay cap, so its package is the aquifer
# itself (6-40 m). Stated here because it is a judgement, not a given.
const WINDOWS = Dict(
    "2P" => (win = (17.5, 33.5), pkg = (4.0, 34.0)),
    "1P" => (win = ( 7.5, 23.5), pkg = (6.0, 40.0)),
)
const W = get(WINDOWS, TAG, (win = (17.5, 33.5), pkg = (4.0, 34.0)))

# zboundaries_aem holds 65 cell TOPS, not 66 edges, so the deepest cell's base
# has to be reconstructed from its centre. The windows here are all shallower
# than 40 m, so that last edge never affects a reported number -- it is done
# properly so the function is not quietly wrong if a window ever reaches it.
"Cell edges (n+1) from cell tops and cell centres."
function cell_edges(ztop::AbstractVector, zcen::AbstractVector)
    n = length(ztop)
    e = Vector{Float64}(undef, n + 1)
    e[1:n] .= ztop
    e[n+1] = ztop[n] + 2 * (zcen[n] - ztop[n])
    e
end

"Conductance in S: integrate sigma over [a,b] using cell edges, clipped."
function conductance(sigma_cells::AbstractVector, edges::AbstractVector, a, b)
    tot = 0.0
    for i in eachindex(sigma_cells)
        lo, hi = edges[i], edges[i+1]
        ov = min(hi, b) - max(lo, a)
        ov > 0 && (tot += sigma_cells[i] * ov)
    end
    tot
end

"Windowed acceptance rates and the T=1 node trace from a costs file."
function chainstats(costs; nchains = 4)
    !isfile(costs) && return nothing
    rows = [split(l) for l in eachline(costs) if !isempty(strip(l))]
    isempty(rows) && return nothing
    length(rows[1]) < 10 && return nothing
    # window boundaries only: stat_window = 100, so these are the rows at which
    # HiQGA recomputed the rate before zeroing its counters.
    bnd = [r for r in rows if (parse(Int, r[2]) - 1) % 100 == 0]
    acc = Dict{String, Vector{Float64}}()
    for (k, c) in zip(["birth","death","position","property","dc"], 3:7)
        v = Float64[]
        for r in bnd
            x = tryparse(Float64, r[c])
            x !== nothing && isfinite(x) && push!(v, x)
        end
        acc[k] = v
    end
    # the same, restricted to windows that ENDED at T=1
    acc1 = Dict{String, Vector{Float64}}()
    b1 = [r for r in bnd if isapprox(parse(Float64, r[10]), 1.0; atol = 1e-12)]
    for (k, c) in zip(["birth","death","position","property","dc"], 3:7)
        v = Float64[]
        for r in b1
            x = tryparse(Float64, r[c])
            x !== nothing && isfinite(x) && push!(v, x)
        end
        acc1[k] = v
    end
    # T=1 node trace, in iteration order
    t1 = [r for r in rows if isapprox(parse(Float64, r[10]), 1.0; atol = 1e-12)]
    sort!(t1, by = r -> parse(Int, r[2]))
    nodes = [parse(Int, r[8]) for r in t1]
    i0 = max(1, floor(Int, 0.5 * length(nodes)))
    nodes_post = nodes[i0:end]
    nchange = length(nodes_post) < 2 ? NaN :
              count(!=(0), diff(nodes_post)) / (length(nodes_post) - 1)
    (acc = acc, acc_T1 = acc1, nodes = nodes_post,
     nchange = nchange, nrows = length(rows), nT1 = length(t1))
end

println("="^78)
@printf("Stage 1 (AEM) convergence diagnostics -- tag=%s setting=%s\n", TAG, SETTING)
@printf("window %.1f-%.1f m   package %.1f-%.1f m\n", W.win..., W.pkg...)
@printf("gate: R-hat < %.2f and ESS >= %d (src/prerun_config.jl)\n",
        CFG.GATE_RHAT, CFG.GATE_MIN_ESS)
println("="^78)

Q2 = Vector{Vector{Float64}}()   # window conductance, per ensemble
Q3 = Vector{Vector{Float64}}()   # package conductance
LS = Vector{Matrix{Float64}}()   # log10 sigma, cells x draws
zall = Float64[]; zbnd = Float64[]
stats = []

for e in 1:NENS
    dir = joinpath(OUT, "aem$(TAG)_$(SETTING)e$(e)")
    if !isdir(dir)
        @printf("ensemble %d: MISSING (%s)\n", e, dir); continue
    end
    S = load(joinpath(dir, "aem_setup.jld2"))
    znall = S["znall_aem"]
    global zall = S["zall_aem"]
    global zbnd = cell_edges(S["zboundaries_aem"], S["zall_aem"])
    nmaxk = haskey(S, "nmax_k") ? S["nmax_k"] : 40
    opt = make_aem_opt(znall, "hypsynth_aem_"; nmax = nmaxk)
    opt.xall[:] .= znall
    old = pwd(); cd(dir)
    F = try
        transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
    finally
        cd(old)
    end
    log10rho = hcat([f[:] for f in F]...)
    sigma = 10 .^ (-log10rho)                       # S/m per AEM cell
    push!(LS, log10.(sigma))
    q2 = [conductance(view(sigma, :, j), zbnd, W.win...) for j in axes(sigma, 2)]
    q3 = [conductance(view(sigma, :, j), zbnd, W.pkg...) for j in axes(sigma, 2)]
    push!(Q2, q2); push!(Q3, q3)
    st = chainstats(joinpath(dir, "misfits_hypsynth_aem_s.bin"))
    push!(stats, st)
    @printf("ensemble %d: %d T=1 draws, nmax=%d\n", e, length(q2), nmaxk)
end

isempty(Q2) && (println("\nNO ENSEMBLES FOUND -- nothing to diagnose."); exit(1))

# ── per-ensemble decision quantities ────────────────────────────────────────
function spread(vs)   # factor: largest / smallest
    m = [median(v) for v in vs]
    (maximum(m) / minimum(m), m)
end
println("\n-- Q2 window conductance (S) --")
@printf("%-10s %10s %10s %10s\n", "ensemble", "p05", "p50", "p95")
for (i, v) in enumerate(Q2)
    @printf("%-10d %10.4f %10.4f %10.4f\n", i, quantile(v,0.05), median(v), quantile(v,0.95))
end
println("\n-- Q3 package conductance (S) --")
@printf("%-10s %10s %10s %10s\n", "ensemble", "p05", "p50", "p95")
for (i, v) in enumerate(Q3)
    @printf("%-10d %10.4f %10.4f %10.4f\n", i, quantile(v,0.05), median(v), quantile(v,0.95))
end

# ── R-hat / ESS across ensembles ────────────────────────────────────────────
println("\n-- convergence across ensembles --")
@printf("%-28s %8s %10s %10s %10s %8s\n",
        "quantity", "R-hat", "ESS bulk", "ESS tail", "med spread", "pass")
function line(name, chains)
    r = rhat_and_ess(chains)
    f, _ = spread(chains)
    ok = isfinite(r.rhat) && r.rhat < CFG.GATE_RHAT &&
         isfinite(r.ess_bulk) && r.ess_bulk >= CFG.GATE_MIN_ESS
    @printf("%-28s %8.3f %10.0f %10.0f %10.2fx %8s\n",
            name, r.rhat, r.ess_bulk, r.ess_tail, f, ok ? "PASS" : "FAIL")
    (name = name, rhat = r.rhat, ess_bulk = r.ess_bulk, ess_tail = r.ess_tail,
     spread = f, pass = ok)
end
rows = Any[]
push!(rows, line("Q2 window conductance", Q2))
push!(rows, line("Q3 package conductance", Q3))

# per-depth, reported separately as the brief asks
println("\n-- per-depth log10 conductivity (reported separately; not decisive) --")
@printf("%-28s %8s %10s %10s %8s\n", "depth", "R-hat", "ESS bulk", "ESS tail", "pass")
for zq in CFG.GATE_DEPTHS
    zq > maximum(zall) && continue
    i = argmin(abs.(zall .- zq))
    ch = [vec(M[i, :]) for M in LS]
    r = rhat_and_ess(ch)
    ok = isfinite(r.rhat) && r.rhat < CFG.GATE_RHAT &&
         isfinite(r.ess_bulk) && r.ess_bulk >= CFG.GATE_MIN_ESS
    @printf("%-28s %8.3f %10.0f %10.0f %8s\n",
            @sprintf("%.0f m", zall[i]), r.rhat, r.ess_bulk, r.ess_tail, ok ? "PASS" : "FAIL")
end

# ── dimension mixing and acceptance ─────────────────────────────────────────
println("\n-- dimension mixing at T=1 (the mechanism) --")
@printf("%-10s %12s %12s %14s\n", "ensemble", "median k", "k range", "k changes/step")
kmed = Float64[]
for (i, st) in enumerate(stats)
    st === nothing && continue
    isempty(st.nodes) && continue
    push!(kmed, median(st.nodes))
    @printf("%-10d %12.1f %12s %13.3f%%\n", i, median(st.nodes),
            "$(minimum(st.nodes))-$(maximum(st.nodes))", 100 * st.nchange)
end
if length(kmed) >= 2
    rk = rhat_and_ess([Float64.(st.nodes) for st in stats if st !== nothing && !isempty(st.nodes)])
    @printf("node count across ensembles: R-hat %.3f, medians %s, spread %.2fx\n",
            rk.rhat, string(round.(kmed, digits=1)), maximum(kmed)/minimum(kmed))
end

println("\n-- windowed acceptance rates, %% (mean over stat windows) --")
@printf("%-10s %9s %9s %9s %9s %9s   %s\n",
        "ensemble", "birth", "death", "position", "property", "dc", "scope")
for (i, st) in enumerate(stats)
    st === nothing && continue
    g(d,k) = isempty(d[k]) ? NaN : mean(d[k])
    @printf("%-10d %9.2f %9.2f %9.2f %9.2f %9.2f   %s\n", i,
            g(st.acc,"birth"), g(st.acc,"death"), g(st.acc,"position"),
            g(st.acc,"property"), g(st.acc,"dc"), "all T")
    @printf("%-10s %9.2f %9.2f %9.2f %9.2f %9.2f   %s\n", "",
            g(st.acc_T1,"birth"), g(st.acc_T1,"death"), g(st.acc_T1,"position"),
            g(st.acc_T1,"property"), g(st.acc_T1,"dc"), "T=1 windows")
end

# ── machine-readable ────────────────────────────────────────────────────────
mkpath(joinpath(HERE, "results"))
csv = joinpath(HERE, "results", "stage1_$(TAG)_$(SETTING).csv")
open(csv, "w") do io
    println(io, "stage,tag,setting,quantity,rhat,ess_bulk,ess_tail,med_spread_factor,pass")
    for r in rows
        @printf(io, "1,%s,%s,%s,%.4f,%.1f,%.1f,%.4f,%s\n",
                TAG, SETTING, r.name, r.rhat, r.ess_bulk, r.ess_tail, r.spread, r.pass)
    end
end
println("\nwrote ", csv)
