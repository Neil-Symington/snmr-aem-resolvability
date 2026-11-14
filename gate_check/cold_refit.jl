# cold_refit.jl -- M1.1 item 3, staged cold refits with an audit.
#
#   julia --project=. gate_check/cold_refit.jl --profile gate_check/results/profile_3P_warm.jld2 \
#         --stage0 <path> --top 7.5 --base 23.5 --mbounds -4,1 --seed 20261114 --tag 3P
#
# Cold starts cost about nine times a warm one (~7 min against ~45 s), so they
# are spent only where they can change an answer:
#
#   1. EDGE nodes: for each rule and each edge, the last two accepted nodes and
#      the first two rejected ones.
#   2. Every node the smoothness test flags, and its two neighbours.
#   3. An AUDIT sample: three unflagged interior nodes, drawn with a logged seed.
#
# Escalation: if any AUDIT node improves by more than 1.0 in chi2 under a cold
# refit, the staged result cannot be trusted and every node on that column needs
# a cold refit. That is reported, not applied silently.
using JLD2, Statistics, Printf, LinearAlgebra, Random, SparseArrays
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "occam_core.jl")); using .OccamCore
const SK = HiQGA.transD_GP.SkyTEM1DInversion

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i],"--") && i < length(ARGS) ? (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (i += 1)
    end
end
arg(k,d) = get(ARG,k,d)

# copied from stages/stage4_gate_v2.jl:62-86; M2 factors this into one place
function operator_from_stage0(A)
    ρ0 = fill(100.0, length(A["z"])); ρ0[A["z"] .< 0] .= 1e12
    aem = SK.dBzdt(; timeslow=A["times_low"], ramplow=A["ramp_low"], zRxlow=A["z_rx"],
        zTxlow=A["z_tx"], timeshigh=A["times_high"], ramphigh=A["ramp_high"],
        zRxhigh=A["z_rx"], zTxhigh=A["z_tx"], rRx=A["r_rx"], rTx=A["r_tx"],
        z=copy(A["z"]), ρ=ρ0, lowpassfcs=A["lowpass_fc"], calcjacobian=true)
    aem.dlow=copy(A["d_low"]); aem.dhigh=copy(A["d_high"])
    aem.σlow=copy(A["sigma_low"]); aem.σhigh=copy(A["sigma_high"])
    aem.selectlow=copy(A["select_low"]); aem.selecthigh=copy(A["select_high"])
    aem.ndatalow=Int(A["ndata_low"]); aem.ndatahigh=Int(A["ndata_high"])
    sg = vcat(aem.σlow[aem.selectlow], aem.σhigh[aem.selecthigh])
    aem.W=sparse(Diagonal(1 ./ sg)); aem.res=zeros(length(sg))
    aem.J=zeros(length(sg), length(aem.ρ)-1); aem
end

haskey(ARG,"mbounds") && (b = parse.(Float64, split(ARG["mbounds"],",")); setbounds!(b[1], b[2]))
PR = load(arg("profile",""))
A  = load(arg("stage0",""))
aem = operator_from_stage0(A)
zall = A["zall"]; n = length(aem.ρ)-1; ndata = A["ndata"]
dzc = diff(vcat(0.0, collect(Float64, zall)))
R = makeR1(n); m0 = fill(-2.0, n); const BETA2 = 0.01

wtop = parse(Float64, arg("top","7.5")); wbot = parse(Float64, arg("base","23.5"))
celltop = vcat(0.0, collect(Float64,zall)[1:end-1]); cellbot = collect(Float64,zall)
ov = [max(0.0, min(cellbot[i],wbot) - max(celltop[i],wtop)) for i in 1:n]
fixidx = findall(>(1e-9), ov)
dzeff = copy(dzc); dzeff[fixidx] .= ov[fixidx]; thick = sum(dzeff[fixidx])

lt = PR["lt"]; warmchi = PR["min_chi2"]; chiref = PR["chiref"]
LAMS = PR["lams"]; NP = length(lt); TAG = arg("tag","col")
tabs = ndata + 3.84; trel = chiref + 3.84

function solve_at(C, l2, ms)
    m, chi, _, nf = occam(aem, ms, m0, R, ndata; fixidx=fixidx, fixval=C, mode=:map,
                          ctype=:cond, dz=dzeff, lambda2=l2, beta2=BETA2,
                          withH=false, info=true)
    (chi=chi, phi=nf.phi, nbound=nf.nbound, conres=nf.conres, m=m)
end

# ── the unconstrained reference model ──────────────────────────────────────
# profile_minchi2.jl stores chiref but not the MODEL, and cold_mref needs the
# model. Recomputing it multi-start costs about three minutes and also re-checks
# chiref against the value the profile recorded.
let t = time()
    Random.seed!(SEED)
    cands = Any[occam(aem, m0, m0, R, ndata; mode = :minimise)]
    for _ in 1:2
        push!(cands, occam(aem, m0 .+ 0.5 .* randn(n), m0, R, ndata; mode = :minimise))
    end
    cs = [c[2] for c in cands]
    global mref0 = cands[argmin(cs)][1]
    global chiref_here = minimum(cs)
    @printf("  chiref recomputed: %s -> %.4f (profile recorded %.4f, diff %+.4f) in %.1f s\n",
            string(round.(cs, digits=3)), chiref_here, chiref, chiref_here - chiref, time()-t)
    flush(stdout)
end

# ── node selection ─────────────────────────────────────────────────────────
const SMOOTH_TOL = 2.0
rough = falses(NP)
for i in 2:NP-1
    pred = warmchi[i-1] + (warmchi[i+1]-warmchi[i-1])*(lt[i]-lt[i-1])/(lt[i+1]-lt[i-1])
    warmchi[i] > pred + SMOOTH_TOL && (rough[i] = true)
end

"last two accepted and first two rejected nodes at each edge, for one threshold"
function edgenodes(thr)
    ok = warmchi .<= thr; out = Int[]
    any(ok) || return out
    i = findfirst(ok); j = findlast(ok)
    append!(out, filter(k -> 1 <= k <= NP, [i, i+1, i-1, i-2]))      # lower edge
    append!(out, filter(k -> 1 <= k <= NP, [j, j-1, j+1, j+2]))      # upper edge
    unique(out)
end
edge_abs = edgenodes(tabs); edge_rel = edgenodes(trel)
flagged  = unique(vcat([[i-1,i,i+1] for i in findall(rough)]...))
flagged  = filter(k -> 1 <= k <= NP, flagged)

SEED = parse(Int, arg("seed","20261114")); Random.seed!(SEED)
interior = setdiff(2:NP-1, union(edge_abs, edge_rel, flagged))
audit = length(interior) >= 3 ? sort(randperm(length(interior))[1:3] .|> i -> interior[i]) : interior

targets = sort(unique(vcat(edge_abs, edge_rel, flagged, audit)))
@printf("cold_refit  tag=%s  seed=%d  profile=%s\n", TAG, SEED, basename(arg("profile","")))
@printf("  N=%d  chiref=%.4f  absolute thr=%.2f  relative thr=%.2f\n", ndata, chiref, tabs, trel)
@printf("  bounds %s ; lambda2 %s ; %d nodes in the profile\n", string(getbounds()), string(LAMS), NP)
@printf("  edge nodes (absolute): %s\n", string(sort(edge_abs)))
@printf("  edge nodes (relative): %s\n", string(sort(edge_rel)))
@printf("  smoothness-flagged + neighbours: %s\n", isempty(flagged) ? "none" : string(sort(flagged)))
@printf("  AUDIT sample (unflagged interior, seed %d): %s\n", SEED, string(audit))
@printf("  -> %d nodes x %d lambda2 x 2 cold starts = %d solves\n",
        length(targets), length(LAMS), length(targets)*length(LAMS)*2)
flush(stdout)

# ── the refits ─────────────────────────────────────────────────────────────
coldchi = fill(NaN, NP); coldsrc = fill("", NP); coldnb = zeros(Int, NP)
t0 = time(); done = 0
for i in targets
    best = Inf; src = ""; nb = 0
    for l2 in LAMS, (sname, st) in (("cold_mref", mref0), ("cold_m0", m0))
        r = solve_at(10^lt[i]*thick, l2, st); done += 1
        if r.conres < 1e-8 && r.chi < best
            best = r.chi; src = sname*"@l2="*string(l2); nb = r.nbound
        end
    end
    coldchi[i] = best; coldsrc[i] = src; coldnb[i] = nb
    @printf("  node %2d  C=%9.4f S  warm %9.3f  cold %9.3f  delta %+8.3f  %-20s (%d solves, %.1f min)\n",
            i, 10^lt[i]*thick, warmchi[i], best, best-warmchi[i], src, done, (time()-t0)/60)
    flush(stdout)
end

# ── escalation ─────────────────────────────────────────────────────────────
auditdelta = [warmchi[i] - coldchi[i] for i in audit if isfinite(coldchi[i])]
worst = isempty(auditdelta) ? 0.0 : maximum(auditdelta)
@printf("\nAUDIT: improvements at %s = %s ; worst %.3f\n", string(audit),
        string(round.(auditdelta, digits=3)), worst)
if worst > 1.0
    @printf("**ESCALATE**: an audit node improved by %.3f > 1.0 under a cold refit.\n", worst)
    println("The staged result cannot be trusted on this column; cold refits are needed at every node.")
else
    println("No audit node improved by more than 1.0. The staged result stands and this is the evidence.")
end

merged = copy(warmchi)
for i in 1:NP; isfinite(coldchi[i]) && coldchi[i] < merged[i] && (merged[i] = coldchi[i]); end
band(c, thr) = (ok = c .<= thr; any(ok) ? (10^lt[findfirst(ok)]*thick, 10^lt[findlast(ok)]*thick) : (NaN,NaN))
@printf("\nband, warm only   : absolute %.4f-%.4f S ; relative %.4f-%.4f S\n",
        band(warmchi,tabs)..., band(warmchi,trel)...)
@printf("band, warm+cold   : absolute %.4f-%.4f S ; relative %.4f-%.4f S\n",
        band(merged,tabs)..., band(merged,trel)...)

jldsave(joinpath(HERE,"results","coldrefit_$(TAG).jld2");
    lt=lt, conductance=10 .^ lt .* thick, warm_chi2=warmchi, cold_chi2=coldchi,
    merged_chi2=merged, cold_source=coldsrc, cold_nbound=coldnb, rough=rough,
    targets=targets, audit=audit, audit_delta=auditdelta, seed=SEED,
    escalate=(worst>1.0), chiref=chiref, ndata=ndata, mbounds=collect(getbounds()))
@printf("\n%d solves, %.1f min\n", done, (time()-t0)/60)
