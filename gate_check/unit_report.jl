# unit_report.jl -- what the rest of the profile has to become for the data to
# accept a given window conductance. Reads stored models only; no solves.
#
#   julia --project=. gate_check/unit_report.jl --profile gate_check/results/profile_X.jld2 \
#         --truth data/truth_murray2P.jld2 --top 17.5 --base 33.5 --tag 2P
#
# For every ACCEPTED node (min chi2 <= N + 3.84) this reports, from that node's
# WINNING model and over the TRUE unit boundaries:
#   mean conductivity and conductance of each unit outside the window
#   min and max cell conductivity in each unit
#   how many cells sit on the parameter bounds
# with truth alongside. REPORT ONLY -- nothing here excludes anything.
using JLD2, Statistics, Printf
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE,".."))
const ARG = Dict{String,String}()
let i=1; while i<=length(ARGS); startswith(ARGS[i],"--") && i<length(ARGS) ?
    (ARG[ARGS[i][3:end]]=ARGS[i+1]; i+=2) : (i+=1); end; end
arg(k,d)=get(ARG,k,d)

PR = load(arg("profile","")); T = load(arg("truth",""))
haskey(PR,"models") || error("""
  this profile artefact stores no models.
  Re-run profile_minchi2.jl at or after the commit that added `models`.""")
M = PR["models"]; lt = PR["lt"]; thick = PR["thick"]; zall = PR["zall"]
ndata = PR["ndata"]; chi = PR["min_chi2"]; mb = PR["mbounds"]
wtop = parse(Float64, arg("top","17.5")); wbot = parse(Float64, arg("base","33.5"))
TAG = arg("tag","col"); tabs = ndata + 3.84

# cell edges and thicknesses on the inversion grid
celltop = vcat(0.0, collect(Float64, zall)[1:end-1]); cellbot = collect(Float64, zall)
dzc = cellbot .- celltop

# TRUE unit boundaries, labelled generically for the paper
zt = T["z"]; sbt = T["σ_bulk"]; dzt = T["dz"]
units = T["units"]
"generic role for a unit, from where it sits relative to the window"
function role(a, b)
    b <= wtop && return "cover/cap"
    a >= wbot && return a < wbot + 20 ? "below window" : "base"
    "WINDOW"
end
roles = [(Float64(u[1]), Float64(u[2]), role(Float64(u[1]), Float64(u[2]))) for u in units]

@printf("unit_report  tag=%s  window %.2f-%.2f m  bounds log10 sigma %s\n", TAG, wtop, wbot, string(mb))
@printf("accepted = min chi2 <= N + 3.84 = %.2f\n\n", tabs)
println("TRUTH, over the true unit boundaries:")
@printf("  %-14s %8s %10s %12s\n", "role", "top-base", "mean sigma", "conductance")
for (a,b,r) in roles
    k = findall(q -> a <= q < b, zt); isempty(k) && continue
    @printf("  %-14s %8s %10.5f %12.4f\n", r, @sprintf("%.0f-%.0f", a, b),
            mean(sbt[k]), sum(sbt[k].*dzt))
end

acc = findall(c -> isfinite(c) && c <= tabs, chi)
@printf("\n%d of %d nodes accepted\n", length(acc), length(lt))
for i in acc
    m = M[i]; any(isnan, m) && continue
    σ = 10 .^ m
    nb_lo = count(x -> x <= mb[1] + 1e-9, m); nb_hi = count(x -> x >= mb[2] - 1e-9, m)
    @printf("\n--- node %d: window C = %.4f S, min chi2 = %.3f ; cells on bounds: %d low, %d high\n",
            i, 10^lt[i]*thick, chi[i], nb_lo, nb_hi)
    @printf("  %-14s %8s %10s %12s %10s %10s %7s\n",
            "role","top-base","mean sigma","conductance","min sigma","max sigma","on bnd")
    for (a,b,r) in roles
        k = findall(j -> cellbot[j] > a && celltop[j] < b, eachindex(dzc))
        isempty(k) && continue
        ovl = [min(cellbot[j],b) - max(celltop[j],a) for j in k]
        @printf("  %-14s %8s %10.5f %12.4f %10.5f %10.5f %7d\n", r,
                @sprintf("%.0f-%.0f", a, b),
                sum(σ[k].*ovl)/sum(ovl), sum(σ[k].*ovl),
                minimum(σ[k]), maximum(σ[k]),
                count(j -> m[j] <= mb[1]+1e-9 || m[j] >= mb[2]-1e-9, k))
    end
end
