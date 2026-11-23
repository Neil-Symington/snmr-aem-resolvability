#!/usr/bin/env julia
# uniform_bisect.jl -- locate a uniform-window edge to two significant figures.
#
#   julia --project=. gate_check/uniform_bisect.jl --col 3S --out <folder> \
#         --upper 1.1246,1.6 --lower 0.284525,0.505964
#
# Each --end takes "admissible,rejected" -- the order is the bracket, not the
# magnitude, so an upper edge reads low,high and a lower edge reads high,low.
# Every trial is a full uniform-window ladder (10, 1, 0.1, 0.01) with the stall
# rule, exactly as the profile ran.
#
# THE BAND RULE (23 Nov). A trial within 2.0 of its threshold is re-solved from
# two further starts drawn inside the bounds; the three must agree to 0.5 or the
# trial is UNDECIDED and the bisection stops with the bracket it has.
#
# Both rules are located from the SAME cache of solves: absolute (N + 3.84) and
# relative (chi2ref + 3.84, chi2ref the profile's minimum for this column).
using JLD2, Printf, LinearAlgebra, SparseArrays, Random, Statistics
const HERE=@__DIR__; const ROOT=normpath(joinpath(HERE,"..")); cd(ROOT)
include(joinpath(ROOT,"src","aem_setup.jl"))
include(joinpath(ROOT,"src","occam_core.jl")); using .OccamCore
const SK = HiQGA.transD_GP.SkyTEM1DInversion

const ARG=Dict{String,String}()
let i=1; while i<=length(ARGS)
  startswith(ARGS[i],"--") && i<length(ARGS) ? (ARG[ARGS[i][3:end]]=ARGS[i+1]; i+=2) : (i+=1)
end end
arg(k,d)=get(ARG,k,d)
const COL=arg("col","3P"); const OUT=arg("out",""); isempty(OUT) && error("--out required")
const LADDER=(10.0,1.0,0.1,0.01); const BAND=2.0; const AGREE=0.5; const MAXIT=8

const WIN = COL in ("3P","3S") ? (7.5,23.5) : (17.5,33.5)
A=load("outputs/synth_$COL/stage0_synth_$COL.jld2")
ρ0=fill(100.0,length(A["z"])); ρ0[A["z"].<0].=1e12
aem=SK.dBzdt(; timeslow=A["times_low"],ramplow=A["ramp_low"],zRxlow=A["z_rx"],
  zTxlow=A["z_tx"],timeshigh=A["times_high"],ramphigh=A["ramp_high"],
  zRxhigh=A["z_rx"],zTxhigh=A["z_tx"],rRx=A["r_rx"],rTx=A["r_tx"],
  z=copy(A["z"]),ρ=ρ0,lowpassfcs=A["lowpass_fc"],calcjacobian=true)
aem.dlow=copy(A["d_low"]); aem.dhigh=copy(A["d_high"])
aem.σlow=copy(A["sigma_low"]); aem.σhigh=copy(A["sigma_high"])
aem.selectlow=copy(A["select_low"]); aem.selecthigh=copy(A["select_high"])
aem.ndatalow=Int(A["ndata_low"]); aem.ndatahigh=Int(A["ndata_high"])
sg=vcat(aem.σlow[aem.selectlow],aem.σhigh[aem.selecthigh])
aem.W=sparse(Diagonal(1 ./ sg)); aem.res=zeros(length(sg)); aem.J=zeros(length(sg),length(aem.ρ)-1)
setbounds!(-3.0,0.30103)
zall=A["zall"]; n=length(aem.ρ)-1; ndata=A["ndata"]
dzc=diff(vcat(0.0,collect(Float64,zall))); R=makeR1(n); m0=fill(-2.0,n)
wt,wb=WIN; ct=vcat(0.0,collect(Float64,zall)[1:end-1]); cb=collect(Float64,zall)
ov=[max(0.0,min(cb[i],wb)-max(ct[i],wt)) for i in 1:n]
fixidx=findall(>(1e-9),ov); dzeff=copy(dzc); dzeff[fixidx].=ov[fixidx]
const L=sum(dzeff[fixidx]); const LO,HI=getbounds()
const THR_ABS=ndata+3.84

# chi2ref from the profile this column already ran
sh=filter(f->startswith(f,"shard_$(COL)_"), readdir(OUT))
# IN A FUNCTION, NOT A TOP-LEVEL LOOP. `cref = min(cref, ...)` inside a bare
# `for` assigns to a NEW LOCAL -- Julia's soft scope -- so cref stayed Inf and
# the relative threshold would have been Inf, making every value admissible
# under that rule. It threw instead of running, which is the good outcome.
function profile_chiref(shards)
    c = Inf
    for f in shards
        d = load(joinpath(OUT, f))
        for i in eachindex(d["chi2"])
            isnan(d["lambda2"][i]) && continue
            d["lambda2"][i] < 0.011 && isfinite(d["chi2"][i]) && (c = min(c, d["chi2"][i]))
        end
    end
    c
end
cref = profile_chiref(sh)
isfinite(cref) || error("no finite profile chi2 for $COL; cannot set the relative threshold")
const THR_REL=cref+3.84
@printf("uniform_bisect %s: thresholds absolute %.2f, relative %.3f (chi2ref %.4f)\n",
        COL, THR_ABS, THR_REL, cref); flush(stdout)

MODELS=Dict{String,Vector{Float64}}(); CACHE=Dict{Tuple{Float64,Int},Float64}()
function ladder(C, st)
  haskey(CACHE,(C,st)) && return CACHE[(C,st)]
  v=log10(C/L)
  (LO-1e-12<=v<=HI+1e-12) || (CACHE[(C,st)]=NaN; return NaN)
  Random.seed!(20261116+1000*st)
  m = st==1 ? fill(-2.0,n) : (LO .+ (HI-LO).*rand(n))
  chi=NaN
  for l2 in LADDER
    try
      m,chi,_,nf = occam(aem,m,m0,R,ndata; fixidx=fixidx,fixval=v,mode=:map,
        ctype=:uniform,dz=dzeff,lambda2=l2,beta2=0.01,withH=false,info=true,
        nsteps=400,chistall=0.01)
      MODELS["$(C)_$(st)_$(l2)"]=copy(m)
    catch e
      @printf("    C=%.6g start %d l2=%.2f THREW: %s\n",C,st,l2,first(sprint(showerror,e),50))
      chi=NaN; break
    end
  end
  CACHE[(C,st)]=chi
end

"chi2 at C with the band rule; returns (chi, decided, starts)"
function evaluate(C, thr)
  c1=ladder(C,1); isnan(c1) && return (NaN,false,[c1])
  abs(c1-thr) > BAND && return (c1,true,[c1])
  cs=[c1, ladder(C,2), ladder(C,3)]
  good=filter(isfinite,cs)
  length(good)<3 && return (c1,false,cs)
  (median(good), (maximum(good)-minimum(good))<=AGREE, cs)
end

sig2(x)=round(x, sigdigits=2)

"bisect between an admissible and a rejected conductance"
function bisect(Ca, Cr, thr, label)
  @printf("  %s edge, threshold %.3f: bracket %.6g (adm) / %.6g (rej)\n", label, thr, Ca, Cr)
  starts=Any[]; undecided=false
  for it in 1:MAXIT
    if sig2(Ca)==sig2(Cr)
      @printf("    bracket agrees to 2 significant figures at %.6g\n", sig2(Ca)); break
    end
    mid = 10^((log10(Ca)+log10(Cr))/2)
    chi, decided, cs = evaluate(mid, thr)
    push!(starts, (C=mid, chi2=chi, decided=decided, starts=cs))
    @printf("    trial %.6g -> chi2 %s %s\n", mid,
            isnan(chi) ? "NaN" : @sprintf("%.4f",chi),
            decided ? (chi<=thr ? "admissible" : "rejected") : "UNDECIDED")
    flush(stdout)
    if isnan(chi) || !decided; undecided=true; break; end
    chi<=thr ? (Ca=mid) : (Cr=mid)
  end
  edge = 10^((log10(Ca)+log10(Cr))/2)
  (edge=sig2(edge), edge_raw=edge, bracket=(Ca,Cr), undecided=undecided, starts=starts)
end

res=Dict{String,Any}()
for (endname, key) in (("upper","upper"), ("lower","lower"))
  spec=arg(endname,""); isempty(spec) && continue
  a,r = parse.(Float64, split(spec,","))
  ra = bisect(a,r,THR_ABS,"$endname absolute")
  rr = bisect(a,r,THR_REL,"$endname relative")
  res["edge_abs_$key"]=ra.edge; res["edge_rel_$key"]=rr.edge
  res["bracket_$key"]=collect(ra.bracket)
  res["starts_$key"]=[(C=s.C,chi2=s.chi2,decided=s.decided,starts=s.starts) for s in vcat(ra.starts,rr.starts)]
  res["undecided_$key"]=ra.undecided || rr.undecided
  @printf("  => %s: absolute %.6g S, relative %.6g S%s\n\n", endname, ra.edge, rr.edge,
          (ra.undecided||rr.undecided) ? "  [bisection stopped UNDECIDED]" : "")
end

prof=Dict{Float64,Float64}()
for f in sh
  d=load(joinpath(OUT,f))
  for i in eachindex(d["C"])
    isnan(d["lambda2"][i]) && continue
    d["lambda2"][i]<0.011 && (prof[d["C"][i]]=d["chi2"][i])
  end
end
adm=[k for k in keys(prof) if prof[k]<=THR_ABS]
floorC=10.0^LO*L
res["end_lower"] = haskey(res,"edge_abs_lower") ? "data" :
    (isapprox(minimum(adm),floorC;rtol=1e-3) ? "stated limit" : "unresolved")
res["end_upper"] = haskey(res,"edge_abs_upper") ? "data" : "unresolved"
res["chiref"]=cref; res["threshold_abs"]=THR_ABS; res["threshold_rel"]=THR_REL
res["col"]=COL; res["L"]=L; res["window"]=collect(WIN); res["models"]=MODELS

# ceilings at the upper edge, porosity from max(0.25, SNMR p05) to 0.40
const SNMR_P05=Dict("3P"=>0.2760,"3S"=>0.2601,"2Pv2"=>0.2466,"2Sv2"=>0.2413)
if haskey(res,"edge_abs_upper")
  C=res["edge_abs_upper"]; plo=max(0.25,SNMR_P05[COL]); mid=0.5*(plo+0.40)
  ceil_(φ,m)=(C/L)/φ^m*6400.0
  res["ceilings"]=[ceil_(0.40,1.50), ceil_(mid,1.70), ceil_(plo,1.90)]
  res["ceiling_porosity"]=[plo,0.40]
  @printf("  ceilings at %.6g S (porosity %.4f-0.40): %.0f / %.0f / %.0f mg/L\n",
          C, plo, res["ceilings"]...)
end

f=joinpath(OUT,"bisect_$(COL).jld2")
# jldsave takes keyword args, so a String-keyed Dict must be converted
jldsave(f; (Symbol(k) => v for (k,v) in res)...)
open(joinpath(OUT,"bisect_$(COL).md"),"w") do io
  println(io,"# Uniform-window edges, $COL\n")
  @printf(io,"thresholds: absolute %.2f, relative %.3f (chi2ref %.4f)\n\n",THR_ABS,THR_REL,cref)
  for e in ("upper","lower")
    haskey(res,"edge_abs_$e") || continue
    @printf(io,"- %s edge: absolute %.6g S, relative %.6g S, bracket %.6g-%.6g%s\n",
      e,res["edge_abs_$e"],res["edge_rel_$e"],res["bracket_$e"]...,
      res["undecided_$e"] ? "  [UNDECIDED]" : "")
  end
  @printf(io,"\nends: lower %s, upper %s\n",res["end_lower"],res["end_upper"])
  haskey(res,"ceilings") && @printf(io,"ceilings: %.0f / %.0f / %.0f mg/L\n",res["ceilings"]...)
  @printf(io,"models saved: %d\n",length(MODELS))
end
println("wrote ",f)
