# sensitivity_depth.jl -- the depth below which the AEM data have negligible
# sensitivity, from the noise-weighted Jacobian. No inversion; one forward
# solve with the Jacobian computed.
#
#   julia --project=. gate_check/sensitivity_depth.jl --stage0 <path> \
#         --model <profile.jld2> --frac 0.01 --tag 3P
#
# DEFINITION. For a reference model, s_j = || W * J[:, j] ||_2, the noise-weighted
# sensitivity of all the data to cell j. Normalise by max(s). The sensitivity
# depth is the shallowest depth below which the normalised sensitivity stays
# under `frac` for every deeper cell. W is included because an unweighted
# Jacobian measures response, not information: a gate with large noise
# contributes little regardless of how strongly it responds.
#
# Cells below this depth are EXCLUDED FROM PINNED COUNTS and not fixed -- see
# DECISIONS.md for why. The count above and below is reported separately so the
# two can be recombined.
using JLD2, Statistics, Printf, LinearAlgebra, SparseArrays
const HERE=@__DIR__; const ROOT=normpath(joinpath(HERE,".."))
include(joinpath(ROOT,"src","aem_setup.jl"))
include(joinpath(ROOT,"src","occam_core.jl")); using .OccamCore
const SK = HiQGA.transD_GP.SkyTEM1DInversion
const ARG=Dict{String,String}()
let i=1; while i<=length(ARGS); startswith(ARGS[i],"--") && i<length(ARGS) ?
    (ARG[ARGS[i][3:end]]=ARGS[i+1]; i+=2) : (i+=1); end; end
arg(k,d)=get(ARG,k,d)

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

A = load(arg("stage0","")); aem = operator_from_stage0(A)
zall = collect(Float64, A["zall"]); n = length(aem.ρ) - 1
celltop = vcat(0.0, zall[1:end-1]); cellbot = zall
FRAC = parse(Float64, arg("frac","0.01")); TAG = arg("tag","col")

# reference model: a stored one if given, else a uniform halfspace
m = if haskey(ARG,"model") && isfile(ARG["model"])
    PR = load(ARG["model"]); ch = haskey(PR,"min_chi2") ? PR["min_chi2"] : PR["chi2"]
    ok = findall(i -> isfinite(ch[i]) && !any(isnan, PR["models"][i]) &&
                      !all(x -> x == -2.0, PR["models"][i]), eachindex(ch))
    isempty(ok) ? fill(-2.0, n) : PR["models"][ok[argmin(ch[ok])]]
else
    fill(-2.0, n)
end
@printf("sensitivity_depth  tag=%s  reference = %s  frac = %.3f\n", TAG,
        haskey(ARG,"model") ? "best stored model" : "uniform halfspace", FRAC)

chi2(aem, m; computeJ = true)
WJ = aem.W * aem.J
s = [norm(view(WJ, :, j)) for j in 1:n]
sn = s ./ maximum(s)

# shallowest depth below which EVERY deeper cell is under frac. In a function:
# a top-level `while` that decrements a counter is the soft-scope trap again.
function sensdepth(sn, cellbot, frac)
    n = length(sn); idx = n
    while idx >= 1 && sn[idx] < frac; idx -= 1; end
    idx >= n ? Inf : cellbot[idx]
end
zsens = sensdepth(sn, cellbot, FRAC)
@printf("  %d cells; peak sensitivity at %.1f m\n", n, cellbot[argmax(sn)])
@printf("  SENSITIVITY DEPTH = %.2f m   (every cell below it is under %.1f%% of peak)\n",
        zsens, 100FRAC)
@printf("  cells below it: %d of %d\n\n", count(>(zsens), cellbot), n)
@printf("%-8s %10s %12s %10s\n","cell","top-bot m","norm. sens","below?")
for j in 1:n
    (j <= 6 || j > n-8 || (sn[j] >= FRAC && sn[min(j+1,n)] < FRAC)) || continue
    @printf("%-8d %10s %12.5f %10s\n", j,
            @sprintf("%.1f-%.1f", celltop[j], cellbot[j]), sn[j],
            cellbot[j] > zsens ? "yes" : "")
end
jldsave(joinpath(HERE,"results","sensdepth_$(TAG).jld2");
    zall=zall, celltop=celltop, cellbot=cellbot, sens=s, sens_norm=sn,
    frac=FRAC, zsens=zsens, ncells_below=count(>(zsens), cellbot))
println("\nwrote results/sensdepth_$(TAG).jld2")
