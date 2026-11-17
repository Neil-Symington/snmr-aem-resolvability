# truecheck.jl -- the true-model check for the gate.
#
#   julia --project=. gate_check/truecheck.jl --stage0 <path> --truth <path> \
#         --top 7.5 --base 23.5 --mbounds -3.0,0.30103 --tag 3P
#
# Two numbers, under the release (SEDIMENTARY) parameter bounds:
#
#   (a) chi2 of the DISCRETISED TRUTH. The truth lives on a 0.25 m grid; the
#       inversion grid has 65 cells. The truth is averaged onto the inversion
#       grid ARITHMETICALLY in conductivity, which is the average that preserves
#       conductance -- the quantity the gate forces. The geometric average is
#       reported alongside because the two differ by about a factor of 1.6 on
#       these columns, and that difference is not negligible against the
#       threshold.
#
#   (b) MIN CHI2 AT C_true, seeded from the discretised truth. This separates
#       two failures that look identical in the profile: a truth the data
#       genuinely reject, and a truth the solver cannot reach from a cold start.
#       If (a) is above the threshold but (b) is below it, the forced-conductance
#       profile at C_true was start-limited, not data-limited.
#
# TRAP, hit once here on 17 Nov: occam's `fixval` under ctype = :cond is the
# CONDUCTANCE IN S (occam_core.jl:172), not the window conductivity. Passing
# C_true/thick forces 1/16 of the intended conductance, which on 3P fell below
# the floor-limited minimum of 0.0160 S; the solve then pinned all nine window
# cells to the 0.001 S/m floor, landed at 0.0160 S, reported chi2 = 25.998 and
# settled = true, and looked like a clean acceptance at C_true. Nothing in the
# solver's return said the constraint had been missed by a factor of 13. So this
# script now recomputes the achieved conductance from the RETURNED model and
# discards any solve that is more than 1% off, rather than trusting settled.
#
# No tuning toward the truth: nothing here feeds back into the band. The band
# stays what profile_minchi2.jl computed from starts that never see the truth.
# This script only says whether the truth is reachable and at what cost.
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

# copied from gate_check/profile_minchi2.jl so the operator is built identically
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

"Average the truth onto the inversion grid. kind = :arith (conductance
preserving) or :geom. Cells deeper than the truth grid take its deepest value."
function discretise(ztruth, sb, celltop, cellbot, kind)
    n = length(cellbot); out = zeros(n)
    zmaxt = maximum(ztruth)
    for i in 1:n
        idx = findall(j -> celltop[i] <= ztruth[j] < cellbot[i], eachindex(ztruth))
        if isempty(idx)
            # below the truth grid, or a cell thinner than the truth spacing
            out[i] = celltop[i] >= zmaxt ? sb[end] :
                     sb[argmin(abs.(ztruth .- 0.5*(celltop[i]+cellbot[i])))]
        else
            out[i] = kind === :arith ? mean(sb[idx]) : exp(mean(log.(sb[idx])))
        end
    end
    out
end

function main()
    A = load(arg("stage0", ""))
    aem = operator_from_stage0(A)
    zall = A["zall"]; n = length(aem.ρ) - 1; ndata = A["ndata"]
    dzc = diff(vcat(0.0, collect(Float64, zall)))
    celltop = vcat(0.0, collect(Float64, zall)[1:end-1]); cellbot = collect(Float64, zall)
    R = makeR1(n); m0 = fill(-2.0, n)

    wtop = parse(Float64, arg("top","7.5")); wbot = parse(Float64, arg("base","23.5"))
    ov = [max(0.0, min(cellbot[i], wbot) - max(celltop[i], wtop)) for i in 1:n]
    fixidx = findall(>(1e-9), ov)
    dzeff = copy(dzc); dzeff[fixidx] .= ov[fixidx]
    thick = sum(dzeff[fixidx])

    if haskey(ARG, "mbounds")
        a, b = parse.(Float64, split(ARG["mbounds"], ","))
        setbounds!(a, b)
    end
    TAG = arg("tag","col")
    tf = arg("truth","")
    T = load(tf); zt = T["z"]; dzt = T["dz"]; sb = T["σ_bulk"]
    Ctrue = sum(sb[i]*dzt for i in eachindex(zt) if wtop <= zt[i] < wbot)

    @printf("truecheck  tag=%s  stage0=%s\n", TAG, basename(arg("stage0","")))
    @printf("  truth %s  hash %s\n", basename(tf), filesha(tf)[1:12])
    @printf("  N = %d data ; window %.2f-%.2f m (support %.4f m, %d of %d cells)\n",
            ndata, wtop, wbot, thick, length(fixidx), n)
    @printf("  parameter bounds, log10 sigma per cell: %s -> sigma %.1e to %.1e S/m\n",
            string(getbounds()), 10.0^getbounds()[1], 10.0^getbounds()[2])
    @printf("  truth window conductance C_true = %.4f S  (sigma_eff %.5f S/m)\n",
            Ctrue, Ctrue/thick)
    tabs = ndata + 3.84
    @printf("  absolute threshold N+3.84 = %.3f\n", tabs)
    flush(stdout)

    # ── (a) chi2 of the discretised truth ─────────────────────────────────
    # The discretisation that MATTERS is the one Stage 0 used to generate the
    # data: build_aem_jac does piecewise_interp(grid.z, truth z, truth σ_bulk)
    # and forward-models THAT, then checks it fits its own noisy data at
    # chi2/N ~ 1. So :interp is the generating model, and its chi2 is the floor
    # any model on this grid is trying to reach -- not a discretisation penalty.
    # :arith and :geom are reported as a sensitivity: they say how much the
    # answer to "does the truth fit?" depends on how the truth is put on the
    # grid, which on 3P turns out to be far more than the 3.84 allowance.
    res = Dict{Symbol,Any}()
    for kind in (:interp, :arith, :geom)
        σd = kind === :interp ? piecewise_interp(collect(Float64, A["z"]), zt, sb)[2:end] :
                                discretise(zt, sb, celltop, cellbot, kind)
        mt = log10.(σd)
        nclip = count(x -> x < getbounds()[1] || x > getbounds()[2], mt)
        mtc = copy(mt); clampm!(mtc)
        c = chi2(aem, mtc)
        Cd = sum(σd[i]*dzeff[i] for i in fixidx)
        res[kind] = (chi = c, m = mtc, C = Cd, nclip = nclip)
        @printf("\n(a) discretised truth, %s\n", kind === :interp ?
                "piecewise_interp onto the grid nodes -- THE GENERATING MODEL, as Stage 0 built it" :
                string(kind) * " average in conductivity (sensitivity only)")
        @printf("    window conductance of the discretised truth: %.4f S (truth %.4f S, factor %.3f)\n",
                Cd, Ctrue, Cd/Ctrue)
        @printf("    cells clipped by the parameter bounds: %d of %d\n", nclip, n)
        @printf("    chi2 = %.3f   -> %s the absolute threshold\n", c,
                c <= tabs ? "BELOW (accepted)" : "ABOVE (rejected)")
        flush(stdout)
    end

    # ── (b) min chi2 at C_true, seeded from the discretised truth ─────────
    # Seeded from the generating model -- the strongest start available, so a
    # rejection here is not a start-limited one.
    @printf("\n(b) min chi2 at C_true = %.4f S, seeded from the GENERATING model\n", Ctrue)
    best = (chi = Inf, nbound = -1, lam = NaN, m = Float64[], C = NaN)
    @printf("    floor-limited minimum reachable conductance (all cells at %.1e S/m over %.2f m) = %.4f S\n",
            10.0^getbounds()[1], thick, 10.0^getbounds()[1]*thick)
    @printf("    ceiling-limited maximum = %.4f S\n", 10.0^getbounds()[2]*thick)
    for l2 in (1.0, 10.0)
        t0 = time()
        m, chi, _, nf = occam(aem, res[:interp].m, m0, R, ndata;
                              fixidx = fixidx, fixval = Ctrue, mode = :map,
                              ctype = :cond, dz = dzeff, lambda2 = l2,
                              beta2 = 0.01, withH = false, info = true)
        Cach = sum(10.0^m[i]*dzeff[i] for i in fixidx)
        rel  = abs(Cach - Ctrue)/Ctrue
        @printf("    lambda2 = %-5.1f chi2 = %10.3f  nbound %3d  settled %s  C achieved %.4f S (%.2f%% off)  conres %.2e  (%.0f s)\n",
                l2, chi, nf.nbound, string(nf.settled), Cach, 100*rel, nf.conres, time()-t0)
        if rel > 0.01
            @printf("      CONSTRAINT NOT MET -- this solve is discarded, not reported as a chi2 at C_true.\n")
        else
            chi < best.chi && (best = (chi = chi, nbound = nf.nbound, lam = l2, m = m, C = Cach))
        end
        flush(stdout)
    end
    if !isfinite(best.chi)
        @printf("    NO SOLVE MET THE CONSTRAINT at C_true. Nothing is reported as accepted or\n")
        @printf("    rejected at C_true: the question is unanswered, not answered negatively.\n")
        @printf("    floor-limited minimum reachable conductance = %.4f S; C_true = %.4f S\n",
                10.0^getbounds()[1]*thick, Ctrue)
    end
    isfinite(best.chi) && @printf("    min chi2 at C_true = %.3f (lambda2 %.1f, %d cells on a bound, C achieved %.4f S) -> %s\n",
            best.chi, best.lam, best.nbound, best.C,
            best.chi <= tabs ? "BELOW the absolute threshold: C_true is ACCEPTED" :
                               "ABOVE the absolute threshold: C_true is REJECTED")
    isfinite(best.chi) && @printf("    cost of forcing the truth, over the discretised truth itself: %+.3f\n",
            best.chi - res[:interp].chi)

    out = joinpath("gate_check", "results", "truecheck_$(TAG).jld2")
    jldsave(out; tag = TAG, truth = tf, truth_hash = filesha(tf)[1:12],
            stage0 = arg("stage0",""), ndata = ndata, thr_abs = tabs,
            window = [wtop, wbot], thick = thick, mbounds = collect(getbounds()),
            Ctrue = Ctrue, C_discretised_arith = res[:arith].C,
            C_discretised_geom = res[:geom].C,
            chi2_truth_interp = res[:interp].chi, C_discretised_interp = res[:interp].C,
            chi2_truth_arith = res[:arith].chi, chi2_truth_geom = res[:geom].chi,
            nclip_interp = res[:interp].nclip,
            nclip_arith = res[:arith].nclip, nclip_geom = res[:geom].nclip,
            chi2_forced_at_Ctrue = best.chi, forced_nbound = best.nbound,
            forced_lambda2 = best.lam, forced_C_achieved = best.C,
            C_floor_limit = 10.0^getbounds()[1]*thick,
            C_ceiling_limit = 10.0^getbounds()[2]*thick,
            m_truth_interp = res[:interp].m, m_truth_arith = res[:arith].m,
            m_forced = best.m)
    println("\nwrote $out")
end

main()
