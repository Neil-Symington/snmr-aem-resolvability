#!/usr/bin/env julia
# uniform_gate.jl -- the UNIFORM-WINDOW gate. One shard of trial conductances.
#
#   julia --project=. gate_check/uniform_gate.jl --col 3P --cs 0.05,0.1 \
#         --out gate_check/results/uniform_<stamp> --starts 1
#
# THE FORMULATION. For a trial conductance C the whole window is held at
# log10(C/L) and REMOVED from the unknowns; only the cells outside it are
# inverted. There is no constraint row, no linearisation of it and no
# projection, so the three things that failed in the constrained solver cannot
# occur: the achieved conductance is exactly C by construction. This is
# OccamCore's existing `ctype = :uniform`, so no solver code changes.
#
# Each value is solved down the ladder 10, 1, 0.1, 0.01, each rung warm-started
# from the last, with the stall rule (chi2 moving less than `chistall` over
# three accepted steps). A ladder whose chi2 is above 2x the threshold and
# whose last rung bought less than a factor of two stops early, marked
# REJECTED-FAR.
#
# EVERY MODEL IS SAVED, per value and per lambda2, so the report pages can show
# them. Record artefacts are never opened for writing.
using JLD2, Printf, LinearAlgebra, SparseArrays, Random, Statistics
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE, ".."))
cd(ROOT)
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "occam_core.jl")); using .OccamCore
const SK = HiQGA.transD_GP.SkyTEM1DInversion

const ARG = Dict{String,String}()
let i=1; while i <= length(ARGS)
    startswith(ARGS[i],"--") && i < length(ARGS) ? (ARG[ARGS[i][3:end]]=ARGS[i+1]; i+=2) : (i+=1)
end end
arg(k,d) = get(ARG,k,d)

const COL    = arg("col","3P")
const CS     = parse.(Float64, split(arg("cs",""), ","))
const OUT    = arg("out","")
const NSTART = parse(Int, arg("starts","1"))
const SHARD  = arg("shard","1")
const LADDER = (10.0, 1.0, 0.1, 0.01)
isempty(OUT) && error("--out is required")
mkpath(OUT)

const WIN = COL in ("3P","3S") ? (7.5, 23.5) : (17.5, 33.5)
A = load("outputs/synth_$COL/stage0_synth_$COL.jld2")
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
aem.J=zeros(length(sg), length(aem.ρ)-1)

setbounds!(-3.0, 0.30103)
zall=A["zall"]; n=length(aem.ρ)-1; ndata=A["ndata"]
dzc=diff(vcat(0.0,collect(Float64,zall))); R=makeR1(n); m0=fill(-2.0,n)
wt,wb = WIN
ct=vcat(0.0,collect(Float64,zall)[1:end-1]); cb=collect(Float64,zall)
ov=[max(0.0,min(cb[i],wb)-max(ct[i],wt)) for i in 1:n]
fixidx=findall(>(1e-9),ov); dzeff=copy(dzc); dzeff[fixidx].=ov[fixidx]
const L = sum(dzeff[fixidx]); const LO, HI = getbounds()
const THR = ndata + 3.84; const FAR = 2*THR

@printf("uniform_gate  col=%s shard=%s  window %.1f-%.1f m (L=%.4f m)\n", COL, SHARD, wt, wb, L)
@printf("  threshold %.2f  bounds %.1e-%.4f S/m  ladder %s  starts %d\n",
        THR, 10.0^LO, 10.0^HI, string(LADDER), NSTART); flush(stdout)

rows = NamedTuple[]; models = Dict{String,Vector{Float64}}()
for C in CS, st in 1:NSTART
    v = log10(C/L)
    if !(LO - 1e-12 <= v <= HI + 1e-12)
        @printf("  C=%.6g  uniform value %.4f is outside the bounds -- UNREACHABLE\n", C, v)
        push!(rows, (C=C, start=st, lambda2=NaN, chi2=NaN, steps=0, nfev=0,
                     converged=false, stalled=false, nbound=0, secs=0.0,
                     status="unreachable")); continue
    end
    Random.seed!(20261116 + 1000*st)
    m = st == 1 ? fill(-2.0, n) : (LO .+ (HI-LO) .* rand(n))
    prev = Inf; tot = 0.0; abandoned = false
    for l2 in LADDER
        t0 = time()
        # A SOLVE THAT THROWS MUST NOT TAKE THE SHARD WITH IT. 2Pv2 at 7.5 and
        # 8.9975 S died with "matrix contains Infs or NaNs" under the
        # constrained solver, and these node lists go to 28.45 S. The failure
        # is recorded as a status and the shard carries on; losing one value
        # is a data point, losing fourteen is an outage.
        local m2, chi, nf
        try
            m2, chi, _, nf = occam(aem, m, m0, R, ndata; fixidx=fixidx, fixval=v,
                                   mode=:map, ctype=:uniform, dz=dzeff, lambda2=l2,
                                   beta2=0.01, withH=false, info=true,
                                   nsteps=400, chistall=0.01)
            m = m2
        catch e
            el = time()-t0; tot += el
            msg = first(sprint(showerror, e), 60)
            @printf("  C=%-10.6g l2=%-6.2f SOLVE THREW: %s\n", C, l2, msg); flush(stdout)
            push!(rows, (C=C, start=st, lambda2=l2, chi2=NaN, steps=0, nfev=0,
                         converged=false, stalled=false, nbound=0, secs=el,
                         status="threw"))
            abandoned = true; break
        end
        el = time()-t0; tot += el
        nb = count(j -> m[j] >= HI-1e-9 || m[j] <= LO+1e-9, 1:n)
        stat = nf.converged ? "converged" : get(nf,:stalled,false) ? "stall" :
               nf.nfev >= 1500 ? "budget" : "other"
        push!(rows, (C=C, start=st, lambda2=l2, chi2=chi, steps=nf.nacc, nfev=nf.nfev,
                     converged=nf.converged, stalled=get(nf,:stalled,false),
                     nbound=nb, secs=el, status=stat))
        models["$(C)_$(st)_$(l2)"] = copy(m)
        @printf("  C=%-10.6g l2=%-6.2f chi2 %-11.4f steps %-5d nfev %-6d %-10s bnd %-3d %.0fs\n",
                C, l2, chi, nf.nacc, nf.nfev, stat, nb, el); flush(stdout)
        if chi > FAR && prev/chi < 2.0
            @printf("  C=%-10.6g REJECTED-FAR after l2=%.2f (chi2 %.1f, last rung x%.2f)\n",
                    C, l2, chi, prev/chi); abandoned = true; break
        end
        prev = chi
    end
    @printf("  C=%-10.6g start %d done in %.0f s%s\n\n", C, st, tot,
            abandoned ? "  [REJECTED-FAR]" : ""); flush(stdout)
end

f = joinpath(OUT, "shard_$(COL)_$(SHARD).jld2")
jldsave(f; col=COL, window=(top=wt,base=wb), L=L, ndata=ndata, threshold=THR,
    ladder=collect(LADDER), bounds=[LO,HI], starts=NSTART,
    C=[r.C for r in rows], start=[r.start for r in rows],
    lambda2=[r.lambda2 for r in rows], chi2=[r.chi2 for r in rows],
    steps=[r.steps for r in rows], nfev=[r.nfev for r in rows],
    converged=[r.converged for r in rows], stalled=[r.stalled for r in rows],
    nbound=[r.nbound for r in rows], secs=[r.secs for r in rows],
    status=[r.status for r in rows], models=models, zall=collect(Float64,zall))
println("wrote ", f)
