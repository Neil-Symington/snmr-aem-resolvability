# edge_bisect.jl -- locate each band edge by warm-started bisection, and SHOW
# that it converged.
#
#   julia --project=. gate_check/edge_bisect.jl --profile <coarse profile.jld2> \
#         --stage0 <path> --top 7.5 --base 23.5 --mbounds -4,1 \
#         --lam 0.01,0.1 --maxrefine 6 --seed 20261116 --tag 3P
#
# WHY BISECTION AND NOT A FINER GRID. The warm chain is a continuation method:
# each solve starts at the previous node's answer, so a shorter step gives a
# better start and a lower chi2. Measured on 3P, the same forced conductance gave
# 47.878 on a 0.1 dex grid and 42.943 on a 0.05 dex grid -- across the 45.84
# threshold. Cold starts did not find it (47.876). Step size is the knob, so the
# edge must be approached by halving the step, not by guessing a grid.
#
# Every loop that carries state lives in a function; see docs/M2_NOTES.md.
using JLD2, Statistics, Printf, LinearAlgebra, Random, SparseArrays
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE,".."))
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "occam_core.jl")); using .OccamCore
const SK = HiQGA.transD_GP.SkyTEM1DInversion

const ARG = Dict{String,String}()
let i=1; while i<=length(ARGS); startswith(ARGS[i],"--") && i<length(ARGS) ?
    (ARG[ARGS[i][3:end]]=ARGS[i+1]; i+=2) : (i+=1); end; end
arg(k,d)=get(ARG,k,d)

function operator_from_stage0(A)   # copied from stages/stage4_gate_v2.jl:62-86
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

let mb = get(ARG, "mbounds", "")
    isempty(mb) || (v = parse.(Float64, split(mb, ",")); setbounds!(v[1], v[2]))
end
const SEED = parse(Int, arg("seed","20261116")); Random.seed!(SEED)
A = load(arg("stage0","")); aem = operator_from_stage0(A)
zall = A["zall"]; n = length(aem.ρ)-1; ndata = A["ndata"]
dzc = diff(vcat(0.0, collect(Float64, zall)))
R = makeR1(n); m0 = fill(-2.0, n); const BETA2 = 0.01
wtop = parse(Float64, arg("top","7.5")); wbot = parse(Float64, arg("base","23.5"))
celltop = vcat(0.0, collect(Float64,zall)[1:end-1]); cellbot = collect(Float64,zall)
ov = [max(0.0, min(cellbot[i],wbot) - max(celltop[i],wtop)) for i in 1:n]
fixidx = findall(>(1e-9), ov)
dzeff = copy(dzc); dzeff[fixidx] .= ov[fixidx]; thick = sum(dzeff[fixidx])
const LAMS = parse.(Float64, split(arg("lam","0.01,0.1"), ","))
const MAXREF = parse(Int, arg("maxrefine","6"))
const TAG = arg("tag","col")

function solve_at(C, l2, ms)
    m, chi, _, nf = occam(aem, ms, m0, R, ndata; fixidx=fixidx, fixval=C, mode=:map,
                          ctype=:cond, dz=dzeff, lambda2=l2, beta2=BETA2,
                          withH=false, info=true)
    (chi=chi, nbound=nf.nbound, conres=nf.conres, m=m)
end

"Best over the lambda2 set at one C, warm-started from `from`."
function best_at(C, from)
    bc = Inf; bm = from; bnb = 0; bl = NaN; k = 0
    for l2 in LAMS
        r = solve_at(C, l2, from); k += 1
        if r.conres < 1e-8 && r.chi < bc
            bc = r.chi; bm = r.m; bnb = r.nbound; bl = l2
        end
    end
    (chi=bc, m=bm, nbound=bnb, lam=bl, nsolves=k)
end

# ── the solved set: log10 sigma -> (chi2, model, nbound) ──────────────────
# The bisection itself lives in gate_check/edgefinder.jl with the solver
# injected, so it is unit-tested against a mock in test/test_edgefinder.jl
# without any AEM code. Two bugs found on 17 Nov are fixed there: a bracket
# invalidated by its own re-solve, and a chi2 convergence test that compared
# different nodes to each other.
include(joinpath(HERE, "edgefinder.jl")); using .EdgeFinder
const EF = EdgeFinder

"Nearest already-solved node to x, in dex. Used by the refinement audit; the
edge finder has its own copy inside the module."
nearest(S, x) = argmin(k -> abs(k - x), collect(keys(S)))

"Solver for the edge finder: best over the lambda2 set at one C, warm-started."
function edge_solve(x, from)
    r = best_at(10^x * thick, from)
    (chi = r.chi, m = r.m, nbound = r.nbound, lam = r.lam)
end

"Refinement audit: insert a midpoint in three random interior intervals and
re-solve the upper node warm from it. Any improvement over 0.5 is a trigger."
function refinement_audit!(S, thr, nint)
    ks = sort(collect(keys(S)))
    ints = [(ks[i], ks[i+1]) for i in 1:length(ks)-1]
    interior = filter(p -> S[p[1]].chi <= thr && S[p[2]].chi <= thr, ints)
    isempty(interior) && return (picked = Tuple{Float64,Float64}[], deltas = Float64[], nsolves = 0)
    pick = interior[randperm(length(interior))[1:min(nint, length(interior))]]
    deltas = Float64[]; ns = 0
    for (a, b) in pick
        mid = (a + b)/2
        rm = best_at(10^mid * thick, S[nearest(S, mid)].m); ns += rm.nsolves
        rb = best_at(10^b * thick, rm.m); ns += rb.nsolves
        push!(deltas, S[b].chi - rb.chi)
        rb.chi < S[b].chi && (S[b] = (chi = rb.chi, m = rb.m, nbound = rb.nbound))
        S[mid] = (chi = rm.chi, m = rm.m, nbound = rm.nbound)
    end
    (picked = pick, deltas = deltas, nsolves = ns)
end

# ── drive ─────────────────────────────────────────────────────────────────
# --profile takes one or more comma-separated artefacts. The union is seeded,
# keeping the LOWEST chi2 at any node solved more than once, so a coarse grid
# can supply the rejected bracket while a finer one supplies the accepted side.
function seed_nodes(paths, m0)
    S = Dict{Float64,NamedTuple}(); cr = Inf
    for pth in paths
        PR = load(pth)
        lt = PR["lt"]; ch = PR["min_chi2"]; ms = get(PR, "models", nothing)
        cr = min(cr, PR["chiref"])
        for i in eachindex(lt)
            isfinite(ch[i]) || continue
            mm = ms === nothing ? m0 : (any(isnan, ms[i]) ? m0 : ms[i])
            if !haskey(S, lt[i]) || ch[i] < S[lt[i]].chi
                S[lt[i]] = (chi = ch[i], m = mm, nbound = 0)
            end
        end
    end
    (S, cr)
end
S, chiref_unc = seed_nodes(split(arg("profile",""), ","), m0)
chiref = min(chiref_unc, minimum(v.chi for v in values(S)))
RULES = (("absolute", ndata + 3.84), ("relative", chiref + 3.84))

@printf("edge_bisect  tag=%s  seed=%d  lambda2 %s  bounds %s\n", TAG, SEED, string(LAMS), string(getbounds()))
@printf("  seeded with %d solved nodes from %s\n", length(S), arg("profile",""))
@printf("  chiref: unconstrained %.4f, used %.4f\n\n", chiref_unc, chiref)
function drive_edges!(S, RULES, thick, ndata)
    results = Dict{String,Any}(); total = 0
    for (rname, thr) in RULES
        @printf("=== %s rule, threshold %.3f\n", rname, thr)
        ks = sort(collect(keys(S)))
        acc = [k for k in ks if S[k].chi <= thr]
        if isempty(acc)
            println("  no accepted node; nothing to bisect"); continue
        end
        for ename in ("upper", "lower")
            out = ename == "upper" ? filter(>(maximum(acc)), ks) : filter(<(minimum(acc)), ks)
            if isempty(out)
                @printf("  %s edge: no rejected node beyond the grid -> OPEN, not bisected\n", ename)
                results["$(rname)_$(ename)"] = (edge = NaN, converged = false, state = "open")
                continue
            end
            inner = ename == "upper" ? maximum(acc) : minimum(acc)
            outer = ename == "upper" ? minimum(out) : maximum(out)
            @printf("  %s edge bracketed by %.4f (chi2 %.3f) and %.4f (chi2 %.3f)\n",
                    ename, inner, S[inner].chi, outer, S[outer].chi)
            r = EF.bisect_edge(S, edge_solve, thr; maxref = MAXREF)
            total += r.nsolves
            @printf("    %-9s %-10s %-10s %-10s %-9s %s\n",
                    "step dex","edge dex","chi2 lo","chi2 hi","move","dchi")
            for h in r.hist
                @printf("    %-9.4f %-10.4f %-10.3f %-10.3f %-9s %s\n",
                        h.step, h.edge, h.chi_lo, h.chi_hi,
                        isnan(h.move) ? "-" : @sprintf("%.4f", h.move),
                        isnan(h.dchi) ? "-" : @sprintf("%.3f", h.dchi))
            end
            @printf("    -> C = %s  %s  state %s  (%d solves)\n",
                    isnan(r.edge) ? "(none)" : @sprintf("%.4f S", 10^r.edge*thick),
                    r.converged ? "CONVERGED" : "**NOT CONVERGED**", r.state, r.nsolves)
            if !isempty(r.raised)
                println("    REFINEMENT RAISED chi2 (should be impossible):")
                for m in r.raised; println("      ", m); end
            end
            results["$(rname)_$(ename)"] = (edge = isnan(r.edge) ? NaN : 10^r.edge*thick,
                                            converged = r.converged, state = r.state,
                                            raised = r.raised)
        end
        println()
    end
    (results, total)
end

results, total = drive_edges!(S, RULES, thick, ndata)

aud = refinement_audit!(S, ndata + 3.84, 3)
@printf("REFINEMENT AUDIT (seed %d): %d intervals, improvements %s\n", SEED,
        length(aud.picked), string(round.(aud.deltas, digits=3)))
worst = isempty(aud.deltas) ? 0.0 : maximum(aud.deltas)
println(worst > 0.5 ?
    "**TRIGGER**: an interior node improved by $(round(worst,digits=3)) > 0.5 on refinement; the grid is too coarse." :
    "No interior node improved by more than 0.5. The grid is fine enough away from the edges.")

ks = sort(collect(keys(S)))
jldsave(joinpath(HERE,"results","edges_$(TAG).jld2");
    lt = ks, chi2 = [S[k].chi for k in ks], nbound = [S[k].nbound for k in ks],
    models = [S[k].m for k in ks], thick = thick, ndata = ndata,
    chiref = chiref, chiref_unconstrained = chiref_unc, lams = LAMS, seed = SEED,
    audit_intervals = collect(aud.picked), audit_deltas = aud.deltas,
    results = results, mbounds = collect(getbounds()), zall = collect(Float64, zall))
@printf("\n%d solves total (%d in the audit). wrote results/edges_%s.jld2\n",
        total + aud.nsolves, aud.nsolves, TAG)
