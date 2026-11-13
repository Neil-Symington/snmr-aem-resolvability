# diagnose_stage2.jl -- convergence diagnostics for Stage 2 (SNMR), across the
# four independent ensembles the cond file already stores separately.
#
#   julia --project=. convergence_test/diagnose_stage2.jl <condfile> <TAG> <SETTING>
#
# Q1 = mean water content over the representative window (v/v).
# Per-depth diagnostics are printed separately: they are expected to fail below
# the depth limit and decide nothing.
using JLD2, Statistics, Printf

const HERE = normpath(@__DIR__)
const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "prerun_config.jl"))
include("/home/symo/Documents/GA/AEM_SMR_Resolvability_Quickrun/ek_prerun/prerun_lib.jl")
using .PrerunLib
const CFG = PrerunConfig

const CONDF   = length(ARGS) >= 1 ? ARGS[1] : "data/cond_C2_typical2PT.jld2"
const TAG     = length(ARGS) >= 2 ? ARGS[2] : "2P"
const SETTING = length(ARGS) >= 3 ? ARGS[3] : "b1"

const WINDOWS = Dict("2P" => (17.5, 33.5), "1P" => (7.5, 23.5))
const WIN = get(WINDOWS, TAG, (17.5, 33.5))

isfile(CONDF) || (println("missing cond file: ", CONDF); exit(1))
D = load(CONDF)
E = D["ens"]; zsmr = D["zsmr"]

println("="^78)
@printf("Stage 2 (SNMR) convergence diagnostics -- tag=%s setting=%s\n", TAG, SETTING)
@printf("file %s\n", CONDF)
@printf("window %.1f-%.1f m   %d ensembles   gate: R-hat < %.2f, ESS >= %d\n",
        WIN..., length(E), CFG.GATE_RHAT, CFG.GATE_MIN_ESS)
println("="^78)

idx = findall(z -> WIN[1] <= z <= WIN[2], zsmr)
isempty(idx) && (println("window covers no cells"); exit(1))
@printf("window spans %d of %d cells (%.1f-%.1f m actual)\n",
        length(idx), length(zsmr), zsmr[first(idx)], zsmr[last(idx)])

# Q1 per ensemble: mean theta over the window, one value per draw
Q1 = [vec(mean(view(M, idx, :), dims = 1)) for M in E]
@printf("\n-- Q1 mean water content over window (v/v) --\n")
@printf("%-10s %8s %10s %10s %10s\n", "ensemble", "ndraws", "p05", "p50", "p95")
for (i, v) in enumerate(Q1)
    @printf("%-10d %8d %10.4f %10.4f %10.4f\n",
            i, length(v), quantile(v,0.05), median(v), quantile(v,0.95))
end

meds = [median(v) for v in Q1]
p05s = [quantile(v,0.05) for v in Q1]
p95s = [quantile(v,0.95) for v in Q1]
fac(x) = maximum(x) / minimum(x)

r = rhat_and_ess(Q1)
ok = isfinite(r.rhat) && r.rhat < CFG.GATE_RHAT &&
     isfinite(r.ess_bulk) && r.ess_bulk >= CFG.GATE_MIN_ESS
println("\n-- convergence across ensembles --")
@printf("%-28s %8s %10s %10s %10s %10s %8s\n",
        "quantity", "R-hat", "ESS bulk", "ESS tail", "med sprd", "p05 sprd", "pass")
@printf("%-28s %8.3f %10.0f %10.0f %9.2fx %9.2fx %8s\n",
        "Q1 window water content", r.rhat, r.ess_bulk, r.ess_tail,
        fac(meds), fac(p05s), ok ? "PASS" : "FAIL")
@printf("  ensemble medians: %s\n", string(round.(meds, digits=4)))
@printf("  p95 spread: %.2fx\n", fac(p95s))

println("\n-- per-depth water content (separate; expected to fail below the depth limit) --")
@printf("%-10s %8s %10s %10s %8s\n", "depth", "R-hat", "ESS bulk", "ESS tail", "pass")
for zq in CFG.GATE_DEPTHS
    zq > maximum(zsmr) && continue
    i = argmin(abs.(zsmr .- zq))
    ch = [vec(M[i, :]) for M in E]
    rr = rhat_and_ess(ch)
    pp = isfinite(rr.rhat) && rr.rhat < CFG.GATE_RHAT &&
         isfinite(rr.ess_bulk) && rr.ess_bulk >= CFG.GATE_MIN_ESS
    @printf("%-10s %8.3f %10.0f %10.0f %8s\n",
            @sprintf("%.0f m", zsmr[i]), rr.rhat, rr.ess_bulk, rr.ess_tail, pp ? "PASS" : "FAIL")
end
# all depths summary
allr = Float64[]; alle = Float64[]
for i in eachindex(zsmr)
    rr = rhat_and_ess([vec(M[i, :]) for M in E])
    isfinite(rr.rhat) && push!(allr, rr.rhat)
    isfinite(rr.ess_bulk) && push!(alle, rr.ess_bulk)
end
!isempty(allr) && @printf("all %d depths: median R-hat %.3f, max %.3f; median ESS %.0f\n",
                          length(allr), median(allr), maximum(allr), median(alle))

mkpath(joinpath(HERE, "results"))
csv = joinpath(HERE, "results", "stage2_$(TAG)_$(SETTING).csv")
open(csv, "w") do io
    println(io, "stage,tag,setting,quantity,rhat,ess_bulk,ess_tail,med_spread_factor,p05_spread_factor,p95_spread_factor,pass")
    @printf(io, "2,%s,%s,Q1 window water content,%.4f,%.1f,%.1f,%.4f,%.4f,%.4f,%s\n",
            TAG, SETTING, r.rhat, r.ess_bulk, r.ess_tail,
            fac(meds), fac(p05s), fac(p95s), ok)
end
println("\nwrote ", csv)
