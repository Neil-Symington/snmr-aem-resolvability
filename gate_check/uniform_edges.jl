#!/usr/bin/env julia
# uniform_edges.jl -- phase 2 of the uniform-window gate: edges, repeat starts
# and the per-column report, from the shards phase 1 wrote.
#
#   julia --project=. gate_check/uniform_edges.jl --col 3P --out <run folder>
#
# EDGES. Wherever an admissible value sits next to a rejected one, the edge is
# bisected in log10 C to two significant figures under the ABSOLUTE rule; the
# relative-rule edge is computed from the same solved set. Each bisection step
# is a full ladder solve, so the bracket is re-established from everything
# solved, never from the midpoint's verdict alone.
#
# THE BAND RULE (ruled 23 Nov 2026). A solve stopped by the stall rule counts
# for the verdict only if its chi2 is more than 2.0 from the threshold. Inside
# that band three starts must agree to 0.5, or the value is reported
# UNDECIDED. A solve that stopped on the forward budget with chi2 still
# falling is an UPPER BOUND: it can prove admissible, never rejected.
#
# UNRESOLVED, NEVER OPEN. If the value beyond the last admissible one cannot be
# solved -- outside the bounds, or no ladder completes -- the edge is probed
# INWARD and reported `unresolved` if it still cannot be bracketed.
using JLD2, Printf, Statistics, LinearAlgebra, SparseArrays, Random
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE, ".."))
cd(ROOT)
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "occam_core.jl")); using .OccamCore
const SK = HiQGA.transD_GP.SkyTEM1DInversion

const ARG = Dict{String,String}()
let i=1; while i <= length(ARGS)
    startswith(ARGS[i],"--") && i < length(ARGS) ? (ARG[ARGS[i][3:end]]=ARGS[i+1]; i+=2) : (i+=1)
end end
arg(k,d)=get(ARG,k,d)
const COL = arg("col","3P"); const OUT = arg("out","")
const LADDER = (10.0, 1.0, 0.1, 0.01); const BAND = 2.0; const AGREE = 0.5

sh = filter(f -> startswith(f, "shard_$(COL)_"), readdir(OUT))
isempty(sh) && error("no shards for $COL in $OUT")
D = [load(joinpath(OUT, f)) for f in sh]
const L = D[1]["L"]; const THR = D[1]["threshold"]; const WIN = D[1]["window"]
const NDATA = D[1]["ndata"]; const LOB, HIB = D[1]["bounds"]

"final-rung chi2 and status per (C, start), from the shards"
function collect_final()
    best = Dict{Tuple{Float64,Int},NamedTuple}()
    for d in D, i in eachindex(d["C"])
        isnan(d["lambda2"][i]) && continue
        k = (d["C"][i], d["start"][i])
        r = (chi=d["chi2"][i], l2=d["lambda2"][i], status=d["status"][i],
             nfev=d["nfev"][i], steps=d["steps"][i], nbound=d["nbound"][i],
             secs=d["secs"][i])
        # the LAST rung of the ladder is the answer; ladder descends, so keep
        # the smallest lambda2 seen for that (C, start)
        (!haskey(best,k) || r.l2 < best[k].l2) && (best[k] = r)
    end
    best
end
const FIN = collect_final()
Cs = sort(unique(first.(keys(FIN))))

"verdict for one C from its start-1 result and, if needed, repeat starts"
function verdict(C, extra = Dict{Int,Float64}())
    r = get(FIN, (C,1), nothing)
    r === nothing && return ("unreachable", NaN, r)
    chi = r.chi
    if r.status == "budget" && chi > THR
        return ("upper-bound only", chi, r)        # cannot reject on a budget stop
    end
    if abs(chi - THR) <= BAND
        vals = vcat(chi, collect(values(extra)))
        length(vals) < 3 && return ("needs starts", chi, r)
        (maximum(vals) - minimum(vals)) > AGREE && return ("UNDECIDED", chi, r)
        chi = median(vals)
    end
    (chi <= THR ? "ADMISSIBLE" : "REJECTED", chi, r)
end

println("# uniform-window edges, column $COL, threshold $(round(THR,digits=2))")
@printf("# %-12s %-11s %-9s %-7s %s\n","C (S)","chi2","status","nfev","verdict")
adm = Float64[]; rej = Float64[]
for C in Cs
    v, chi, r = verdict(C)
    r === nothing || @printf("  %-12.6g %-11.4f %-9s %-7d %s\n", C, chi, r.status, r.nfev, v)
    r === nothing && @printf("  %-12.6g %-11s %-9s %-7s %s\n", C, "-", "-", "-", v)
    v == "ADMISSIBLE" && push!(adm, C)
    v == "REJECTED"   && push!(rej, C)
end
isempty(adm) && (println("\nNO ADMISSIBLE VALUE -- nothing to bracket"); exit(0))

lo_adm, hi_adm = minimum(adm), maximum(adm)
floorC = 10.0^LOB * L
@printf("\nadmissible as solved: %.6g to %.6g S\n", lo_adm, hi_adm)
@printf("lower end: floor x L = %.4f S -> %s\n", floorC,
        isapprox(lo_adm, floorC; rtol=1e-3) ? "SET BY THE STATED LIMIT" :
        (any(c -> c < lo_adm, Cs) ? "a lower value exists and was not admissible" :
         "no lower value was tried"))
above = filter(>(hi_adm), Cs)
if isempty(above)
    println("upper end: no value above the highest admissible one was tried -> UNRESOLVED")
else
    nxt = minimum(above); v, chi, _ = verdict(nxt)
    @printf("upper end: next value %.6g S is %s -> %s\n", nxt, v,
            v == "REJECTED" ? "BRACKETED, bisect" : "UNRESOLVED (probe inward)")
end
jldsave(joinpath(OUT, "edges_$(COL).jld2"); col=COL, C=Cs,
        chi2=[get(FIN,(c,1),(chi=NaN,)).chi for c in Cs],
        admissible=adm, rejected=rej, threshold=THR, window=WIN, L=L,
        lower_is_stated_limit = isapprox(lo_adm, floorC; rtol=1e-3))
println("\nwrote ", joinpath(OUT, "edges_$(COL).jld2"))
