# smoothest.jl -- the SMOOTHEST ADMISSIBLE model at each accepted node.
#
#   julia --project=. gate_check/smoothest.jl --profile <edges_or_profile.jld2> \
#         --stage0 <path> --top 7.5 --base 23.5 --mbounds -3.0,0.30103 \
#         --lams 10,1,0.1,0.01 --tag 3P_sed
#
# WHY. The band is min chi2 over the lambda2 set, so the model STORED at a node
# is the one from the lightest smoothing that won -- and at the light end the
# profile is spiky. For the paper's figure the useful model is the smoothest one
# that STILL PASSES: it shows what the data require rather than what the weakest
# regularisation permits.
#
# Per accepted node, warm-started from that node's stored model, lambda2 is
# walked DOWN from the heaviest. The first value that both meets the conductance
# constraint and stays under the threshold wins; if none does, the node falls
# back to its stored min-chi2 model and is marked.
#
# Output shape is the one the gate stores use, so a reader can pick it up
# without special-casing: `models`, `conductance`, `lt`, `mbounds`, `ndata`,
# `thick`, plus `lam_used` and `fellback` per node.
using JLD2, Statistics, Printf, LinearAlgebra, Random, SparseArrays
const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "occam_core.jl")); using .OccamCore
const SK = HiQGA.transD_GP.SkyTEM1DInversion

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) ? (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (i += 1)
    end
end
arg(k, d) = get(ARG, k, d)

function operator_from_stage0(A)
    ρ0 = fill(100.0, length(A["z"])); ρ0[A["z"] .< 0] .= 1e12
    aem = SK.dBzdt(; timeslow = A["times_low"], ramplow = A["ramp_low"],
        zRxlow = A["z_rx"], zTxlow = A["z_tx"],
        timeshigh = A["times_high"], ramphigh = A["ramp_high"],
        zRxhigh = A["z_rx"], zTxhigh = A["z_tx"],
        rRx = A["r_rx"], rTx = A["r_tx"], z = copy(A["z"]), ρ = ρ0,
        lowpassfcs = A["lowpass_fc"], calcjacobian = true)
    aem.dlow = copy(A["d_low"]);        aem.dhigh = copy(A["d_high"])
    aem.σlow = copy(A["sigma_low"]);    aem.σhigh = copy(A["sigma_high"])
    aem.selectlow = copy(A["select_low"]); aem.selecthigh = copy(A["select_high"])
    aem.ndatalow = Int(A["ndata_low"]); aem.ndatahigh = Int(A["ndata_high"])
    sg = vcat(aem.σlow[aem.selectlow], aem.σhigh[aem.selecthigh])
    aem.W = sparse(Diagonal(1 ./ sg)); aem.res = zeros(length(sg))
    aem.J = zeros(length(sg), length(aem.ρ) - 1)
    aem
end

function main()
    P = load(arg("profile", ""))
    A = load(arg("stage0", ""))
    aem = operator_from_stage0(A)
    zall = A["zall"]; n = length(aem.ρ) - 1; ndata = Int(P["ndata"])
    dzc = diff(vcat(0.0, collect(Float64, zall)))
    celltop = vcat(0.0, collect(Float64, zall)[1:end-1]); cellbot = collect(Float64, zall)
    R = makeR1(n); m0 = fill(-2.0, n)

    wtop = parse(Float64, arg("top", "7.5")); wbot = parse(Float64, arg("base", "23.5"))
    ov = [max(0.0, min(cellbot[i], wbot) - max(celltop[i], wtop)) for i in 1:n]
    fixidx = findall(>(1e-9), ov)
    dzeff = copy(dzc); dzeff[fixidx] .= ov[fixidx]
    thick = sum(dzeff[fixidx])

    mb = collect(Float64, P["mbounds"]); setbounds!(mb[1], mb[2])
    LAMS = parse.(Float64, split(arg("lams", "10,1,0.1,0.01"), ","))
    TAG = arg("tag", "col")

    lt = collect(Float64, haskey(P, "lt") ? P["lt"] : log10.(P["conductance"] ./ thick))
    chi = collect(Float64, haskey(P, "min_chi2") ? P["min_chi2"] : P["chi2"])
    M = P["models"]
    tabs = ndata + 3.84
    solved(m) = !isempty(m) && !any(isnan, m) && !all(x -> x == -2.0, m)
    acc = [i for i in eachindex(chi) if isfinite(chi[i]) && chi[i] <= tabs && solved(collect(Float64, M[i]))]

    @printf("smoothest  tag=%s  %d accepted of %d nodes  threshold %.2f  lams %s\n",
            TAG, length(acc), length(chi), tabs, string(LAMS))
    flush(stdout)

    out_m = [Float64[] for _ in eachindex(chi)]
    lam_used = fill(NaN, length(chi)); fell = falses(length(chi))
    out_chi = fill(NaN, length(chi)); nsolve = 0
    for i in acc
        start = collect(Float64, M[i]); C = 10^lt[i] * thick
        got = false
        for l2 in LAMS                       # heaviest first
            m, c, _, nf = occam(aem, start, m0, R, ndata; fixidx = fixidx, fixval = C,
                                mode = :map, ctype = :cond, dz = dzeff, lambda2 = l2,
                                beta2 = 0.01, withH = false, info = true)
            nsolve += 1
            if nf.conres < 1e-8 && nf.settled && c <= tabs
                out_m[i] = m; lam_used[i] = l2; out_chi[i] = c; got = true; break
            end
        end
        if !got
            out_m[i] = start; fell[i] = true; out_chi[i] = chi[i]
        end
        @printf("  node %2d  C %8.4f S  lambda2 %-6s chi2 %8.3f %s\n", i, C,
                got ? string(lam_used[i]) : "-", out_chi[i], got ? "" : "(fell back)")
        flush(stdout)
    end

    f = joinpath(ROOT, "gate_check", "results", "smooth_$(TAG).jld2")
    jldsave(f; models = out_m, conductance = (10 .^ lt) .* thick, lt = lt,
            mbounds = mb, ndata = ndata, thick = thick, chi2 = out_chi,
            lam_used = lam_used, fellback = fell, lams = LAMS,
            source = basename(arg("profile", "")), nsolves = nsolve,
            covered = count(i -> !isempty(out_m[i]), eachindex(out_m)), naccepted = length(acc))
    @printf("\nwrote %s  (%d of %d accepted nodes covered, %d solves)\n",
            relpath(f, ROOT), count(i -> !isempty(out_m[i]), acc), length(acc), nsolve)
end

main()
