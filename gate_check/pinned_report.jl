# pinned_report.jl -- M1.1 close-out item 2. From stored models only.
#
#   julia --project=. gate_check/pinned_report.jl --profile <edges_X.jld2> \
#         --truth <truth.jld2> --top 7.5 --base 23.5 --tag 3P
#
# Reports, per accepted node: how many cells sit on a bound, and each unit's
# conductance against truth. Then the largest forced conductance reached with NO
# cell on a bound, and the ceiling that one gives -- because that is the part of
# the band that owes nothing to where the clamps were set.
using JLD2, Statistics, Printf
const HERE=@__DIR__
const ARG=Dict{String,String}()
let i=1; while i<=length(ARGS); startswith(ARGS[i],"--") && i<length(ARGS) ?
    (ARG[ARGS[i][3:end]]=ARGS[i+1]; i+=2) : (i+=1); end; end
arg(k,d)=get(ARG,k,d)
const MGL=6400.0

PR=load(arg("profile","")); T=load(arg("truth",""))
lt=PR["lt"]; thick=PR["thick"]; zall=PR["zall"]; M=PR["models"]
ndata=PR["ndata"]; mb=PR["mbounds"]
chi = haskey(PR,"min_chi2") ? PR["min_chi2"] : PR["chi2"]
chiref = PR["chiref"]
wtop=parse(Float64,arg("top","7.5")); wbot=parse(Float64,arg("base","23.5"))
TAG=arg("tag","col"); tabs=ndata+3.84; trel=chiref+3.84
L=wbot-wtop

celltop=vcat(0.0,collect(Float64,zall)[1:end-1]); cellbot=collect(Float64,zall)
zt=T["z"]; sbt=T["σ_bulk"]; dzt=T["dz"]
roleof(a,b) = b<=wtop ? "cover" : (a>=wbot ? (a<wbot+20 ? "below" : "base") : "WINDOW")
units=[(Float64(u[1]),Float64(u[2]),roleof(Float64(u[1]),Float64(u[2]))) for u in T["units"]]
truthC(a,b)=sum(sbt[i]*dzt for i in eachindex(zt) if a<=zt[i]<b)
function modelC(m,a,b)
    σ=10 .^ m; s=0.0
    for j in eachindex(σ)
        ov=min(cellbot[j],b)-max(celltop[j],a); ov>0 && (s+=σ[j]*ov)
    end
    s
end
pinned(m)=count(x -> x<=mb[1]+1e-9 || x>=mb[2]-1e-9, m)
solved(m)= !any(isnan,m) && !all(x->x==-2.0,m)

@printf("pinned_report  tag=%s  window %.2f-%.2f m  bounds %s\n", TAG, wtop, wbot, string(mb))
@printf("absolute threshold %.2f ; relative %.2f (chiref %.4f)\n\n", tabs, trel, chiref)
@printf("%-10s %9s %7s %8s %-44s\n","C (S)","min chi2","pinned","rule","unit conductance / truth")
rows=NamedTuple[]
for i in sortperm(lt)
    m=M[i]; solved(m) || continue
    chi[i] <= tabs || continue
    rule = chi[i] <= trel ? "both" : "abs"
    parts=String[]
    for (a,b,r) in units
        tc=truthC(a,b); tc<=0 && continue
        push!(parts, @sprintf("%s %.2fx", r, modelC(m,a,b)/tc))
    end
    np=pinned(m)
    @printf("%-10.4f %9.3f %7d %8s %-44s\n", 10^lt[i]*thick, chi[i], np, rule, join(parts,"  "))
    push!(rows,(C=10^lt[i]*thick, chi=chi[i], pin=np, rule=rule))
end

clean=filter(r -> r.pin==0, rows)
println()
if isempty(clean)
    println("NO accepted node has every cell off its bounds.")
else
    best=clean[argmax([r.C for r in clean])]
    @printf("largest forced C with NO cell on a bound: %.4f S  (chi2 %.3f, rule %s)\n", best.C, best.chi, best.rule)
    @printf("  against the absolute edge of the full band, that is the part of the band\n")
    @printf("  that owes nothing to where the clamps were set.\n\n")
    tds(C,θ,m)=(C/L)/θ^m*MGL
    @printf("  ceiling at that C: widest %.0f, central %.0f, tightest %.0f mg/L\n",
            tds(best.C,0.25,1.9), tds(best.C,0.325,1.7), tds(best.C,0.40,1.5))
    @printf("  with the SNMR limit theta>=0.265: widest %.0f mg/L\n", tds(best.C,0.265,1.9))
end
# the relative-rule edge node
rel=filter(r -> r.rule=="both", rows)
isempty(rel) || @printf("\nlargest C accepted under the RELATIVE rule: %.4f S (chi2 %.3f, %d pinned)\n",
                        maximum(r.C for r in rel), rel[argmax([r.C for r in rel])].chi,
                        rel[argmax([r.C for r in rel])].pin)
