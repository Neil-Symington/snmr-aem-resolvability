# profile_minchi2.jl -- M1.1 item 3. The forced-conductance profile under the
# agreed definition:
#
#   the band is MIN CHI2 subject to the forced window conductance and the
#   stated parameter bounds -- not chi2 at the MAP.
#
# At each node the lowest chi2 is kept over every start, both sweep directions
# and the whole lambda2 set. chi2 at the MAP for lambda2 = 1 is recorded
# alongside so the difference is visible. The verdict uses min chi2.
#
#   julia --project=. gate_check/profile_minchi2.jl --stage0 <path> \
#         --top 7.5 --base 23.5 --lo -3.0 --hi -1.0 --step 0.2 \
#         --lam 1,10 --truth data/truth_murray3P.jld2 --tag 3P
#
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
filesha(p) = bytes2hex(reinterpret(UInt8, [hash(read(p))]))

# ── copied VERBATIM from stages/stage4_gate_v2.jl:62-86 ────────────────────
# Copied rather than included because that file ends in main(). When the fix
# lands in the stage, this copy goes away.
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

A = load(arg("stage0", ""))
aem = operator_from_stage0(A)
zall = A["zall"]; n = length(aem.ρ) - 1; ndata = A["ndata"]
dzc = diff(vcat(0.0, collect(Float64, zall)))
R = makeR1(n); m0 = fill(-2.0, n); RtR = Matrix(R'R)
const BETA2 = 0.01

wtop = parse(Float64, arg("top","7.5")); wbot = parse(Float64, arg("base","23.5"))
celltop = vcat(0.0, collect(Float64, zall)[1:end-1]); cellbot = collect(Float64, zall)
ov = [max(0.0, min(cellbot[i], wbot) - max(celltop[i], wtop)) for i in 1:n]
fixidx = findall(>(1e-9), ov)
dzeff = copy(dzc); dzeff[fixidx] .= ov[fixidx]
thick = sum(dzeff[fixidx])

LAMS = parse.(Float64, split(arg("lam","1,10"), ","))
lo = parse(Float64, arg("lo","-3.0")); hi = parse(Float64, arg("hi","-1.0"))
step = parse(Float64, arg("step","0.2"))
lt = collect(lo:step:hi); NP = length(lt)
TAG = arg("tag","col")

# PARAMETER BOUNDS define admissibility. occam_core clamps m to fbounds; report them.
@printf("profile_minchi2  tag=%s  stage0=%s\n", TAG, basename(arg("stage0","")))
@printf("  N = %d data ; window %.2f-%.2f m (support %.4f m, %d of %d cells)\n",
        ndata, wtop, wbot, thick, length(fixidx), n)
@printf("  grid: log10 sigma_eff %.2f to %.2f step %.2f -> %d nodes ; lambda2 = %s\n",
        lo, hi, step, NP, string(LAMS))
# Parameter bounds define admissibility. occam_core holds them in Refs, so
# item 2 (extend the search below 1e-3 S/m) is a setbounds! call, not a config
# edit: the per-cell clamp is what limits how low the window conductance can be
# forced in the first place.
if haskey(ARG, "mbounds")
    a, b = parse.(Float64, split(ARG["mbounds"], ","))
    setbounds!(a, b)
end
@printf("  parameter bounds, log10 sigma per cell: %s  -> sigma %.1e to %.1e S/m\n",
        string(getbounds()), 10.0^getbounds()[1], 10.0^getbounds()[2])
flush(stdout)

"One constrained solve. Returns (chi2, phi, nbound, settled, model)."
function solve_at(C, l2, mstart)
    m, chi, _, nf = occam(aem, mstart, m0, R, ndata; fixidx = fixidx, fixval = C,
                          mode = :map, ctype = :cond, dz = dzeff, lambda2 = l2,
                          beta2 = BETA2, withH = false, info = true)
    (chi = chi, phi = nf.phi, nbound = nf.nbound, settled = nf.settled,
     conres = nf.conres, m = m)
end

# ── the unconstrained reference, multi-start, lowest chi2 ─────────────────
t0 = time()
refs = Any[]
push!(refs, (:m0, occam(aem, m0, m0, R, ndata; mode = :minimise)))
Random.seed!(20261114)
for k in 1:2
    ms = m0 .+ 0.5 .* randn(n)
    push!(refs, (Symbol("rand$k"), occam(aem, ms, m0, R, ndata; mode = :minimise)))
end
chirefs = [r[2][2] for r in refs]
chiref = minimum(chirefs)
@printf("\n  chiref multi-start: %s -> chiref = %.4f  (%s)\n",
        join([@sprintf("%.3f", c) for c in chirefs], ", "), chiref,
        string(refs[argmin(chirefs)][1]))
@printf("  unconstrained reference took %.1f s\n", time()-t0)
flush(stdout)

# ── the profile: both directions, cold starts, all lambda2 ────────────────
best_chi = fill(Inf, NP); best_src = fill("", NP); best_nb = zeros(Int, NP)
best_phi = fill(NaN, NP); map_chi_l1 = fill(NaN, NP)
nsolves = 0
for l2 in LAMS
    for (dirname, order) in (("fwd", 1:NP), ("bwd", NP:-1:1))
        warm = nothing
        for i in order
            C = 10^lt[i] * thick
            starts = Any[("cold_mref", nothing), ("cold_m0", m0)]
            warm !== nothing && push!(starts, ("warm_$dirname", warm))
            for (sname, ms) in starts
                st = ms === nothing ? refs[argmin(chirefs)][2][1] : ms
                r = solve_at(C, l2, st); nsolves += 1
                if r.conres < 1e-8 && r.chi < best_chi[i]
                    best_chi[i] = r.chi; best_phi[i] = r.phi
                    best_nb[i] = r.nbound
                    best_src[i] = "$(sname)@l2=$(l2)"
                end
                sname == "warm_$dirname" && (warm = r.m)
                ms === nothing && warm === nothing && (warm = r.m)
            end
        end
    end
    @printf("  lambda2 %-6g done, %d solves, %.1f min elapsed\n", l2, nsolves, (time()-t0)/60)
    flush(stdout)
end

# chi2 at the MAP for lambda2 = 1, for comparison
if 1.0 in LAMS
    warm = refs[argmin(chirefs)][2][1]
    for i in 1:NP
        r = solve_at(10^lt[i]*thick, 1.0, warm); nsolves += 1
        map_chi_l1[i] = r.chi; warm = r.m
    end
end

# ── report ────────────────────────────────────────────────────────────────
mkpath(joinpath(HERE,"results"))
jldsave(joinpath(HERE,"results","profile_$(TAG).jld2");
    lt = lt, thick = thick, conductance = 10 .^ lt .* thick,
    min_chi2 = best_chi, min_chi2_source = best_src, min_chi2_nbound = best_nb,
    min_phi = best_phi, map_chi2_lambda1 = map_chi_l1,
    chiref = chiref, chirefs = chirefs, ndata = ndata, lams = LAMS,
    window = (top = wtop, base = wbot), fbounds = collect(OccamCore.FBOUNDS),
    nsolves = nsolves, stage0 = arg("stage0",""))

tabs = ndata + 3.84; trel = chiref + 3.84
@printf("\n%-5s %10s %10s %10s %7s %-24s\n","i","C (S)","min chi2","MAP l2=1","nbound","winning start")
for i in 1:NP
    @printf("%-5d %10.4f %10.3f %10.3f %7d %-24s\n", i, 10^lt[i]*thick, best_chi[i],
            map_chi_l1[i], best_nb[i], best_src[i])
end
# consistency check
mn = minimum(best_chi)
@printf("\nCONSISTENCY: chiref %.4f  <=  min over nodes %.4f  (+1e-6) : %s\n",
        chiref, mn, chiref <= mn + 1e-6 ? "PASS" : "**FAIL**")
chiref <= mn + 1e-6 || @warn "chiref exceeds the constrained minimum: the unconstrained solve is not the minimum"

band(thr) = (ok = best_chi .<= thr; any(ok) ?
    (10^lt[findfirst(ok)]*thick, 10^lt[findlast(ok)]*thick) : (NaN, NaN))
ba = band(tabs); br = band(trel)
@printf("absolute band (N+3.84 = %.2f): %.4f - %.4f S\n", tabs, ba...)
@printf("relative band (chiref+3.84 = %.2f): %.4f - %.4f S\n", trel, br...)
tf = arg("truth","")
if !isempty(tf) && isfile(tf)
    T = load(tf); z = T["z"]; dz = T["dz"]; sb = T["σ_bulk"]
    tw = sum(sb[i]*dz for i in eachindex(z) if wtop <= z[i] < wbot)
    @printf("truth window C = %.4f S (truth hash %s)\n", tw, filesha(tf)[1:12])
    @printf("  absolute: %s   relative: %s\n",
            (tw>=ba[1] && tw<=ba[2]) ? "INSIDE" : "**OUTSIDE**",
            (isnan(br[1]) ? "empty band" : ((tw>=br[1] && tw<=br[2]) ? "INSIDE" : "**OUTSIDE**")))
end
@printf("\n%d solves, %.1f min total\n", nsolves, (time()-t0)/60)
