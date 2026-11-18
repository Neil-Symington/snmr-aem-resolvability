# stage4_gate_v2.jl -- Stage 4 of the Section 3 workflow, config-driven.
#
#   julia --project=. stage4_gate_v2.jl --config ../configs/synth_3P.toml
#
# Supersedes stage4_gate.jl, which is kept as the record of the superseded
# numbers. Four things differ.
#
# 1. NO TRUTH IS READ. stage4_gate.jl found the aquifer by scanning the truth
#    file's unit labels for the string "AQUIFER" and trimming 3 m off each end.
#    That is the truth firewall broken in the most consequential place in the
#    workflow: the gate was told where the water is. The window now comes from
#    STAGE 3, which derived it from the SNMR draws, and the AEM data comes from
#    the sealed Stage 0 artefact. Neither reads the column.
#
# 2. PARTIAL CELLS ARE OVERLAP-TRUNCATED, so the forced quantity is exactly the
#    conductance of the Stage 3 window. The AEM grid is log-spaced and its cells
#    do not align with any window: forcing whole cells made a nominal 10.0 m
#    unit carry 8.82 m of support. The constraint row uses an EFFECTIVE dz --
#    the overlap of each cell with the window -- so the support equals the
#    window to machine precision. Asserted below.
#
# 3. THE GRID RUNS TO THE MODEL BOUNDS and admissibility is enforced per node,
#    as the brief specifies, rather than by narrowing the grid. Nodes near the
#    bounds cannot satisfy the constraint (measured: 11 cells pinned and the
#    achieved conductance out by a factor of 9 at the extreme) and are excluded,
#    not normalised over.
#
# 4. TWO OUTPUTS, NOT ONE.
#    (a) PROFILE BAND. chi2 <= N + 3.84 referenced to the TARGET misfit N, edges
#        located by bisection. The upper edge converts to a clay-free salinity
#        through theta (SNMR draws) and m (prior): the DATA-ONLY BOUND, the one
#        number here that owes nothing to the clay prior.
#    (b) BAYESIAN CURVE. Laplace-corrected integrated likelihood,
#          log L(C) = -Phi*/2 - 1/2 logdet(Z'HZ) + 1/2 logdet(Z'PZ),
#        lambda-hat by ML-II, p(C) under a flat-in-log-C prior, and
#        P(TDS on the wrong side of X) integrating theta over the draws and m
#        over its prior. Delta is the MAXIMUM change in that P across
#        lambda-hat/10, lambda-hat, lambda-hat x10 AND across dropping one node
#        at each end of the grid. Delta <= tol -> data-carried, else
#        prior-carried.
using JLD2, Statistics, Printf, LinearAlgebra, Random, TOML, SparseArrays
const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "config.jl")); using .SiteConfig
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "occam_core.jl")); using .OccamCore
# The validated edge search, unit-tested against an injected mock solver in
# test/test_edgefinder.jl. Stage 4 wraps it rather than carrying a second
# implementation; the one it used to carry had both of the bugs that module
# exists to prevent.
include(joinpath(ROOT, "gate_check", "edgefinder.jl")); using .EdgeFinder
const SK = HiQGA.transD_GP.SkyTEM1DInversion

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) ? (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (i += 1)
    end
end
const CFGPATH = get(ARG, "config", "")
isempty(CFGPATH) && error("usage: stage4_gate_v2.jl --config <path>")
const NPTS_OVERRIDE = parse(Int, get(ARG, "npts", "0"))
# Tier F reuses lambda-hat from tier P: the ladder is 7 unconstrained solves per
# realisation and it is selecting a prior weight for the SAME column under a
# different NOISE realisation, which is not a thing the noise should move.
const LHAT_OVERRIDE = parse(Float64, get(ARG, "lhat", "0"))

"""Rebuild the SkyTEM operator from the SEALED artefact. No truth, no column."""
function operator_from_stage0(A)
    ρ0 = fill(100.0, length(A["z"]))        # a neutral halfspace; the data, not
    ρ0[A["z"] .< 0] .= 1e12                 # the start, carries the information
    aem = SK.dBzdt(; timeslow = A["times_low"], ramplow = A["ramp_low"],
        zRxlow = A["z_rx"], zTxlow = A["z_tx"],
        timeshigh = A["times_high"], ramphigh = A["ramp_high"],
        zRxhigh = A["z_rx"], zTxhigh = A["z_tx"],
        rRx = A["r_rx"], rTx = A["r_tx"], z = copy(A["z"]), ρ = ρ0,
        lowpassfcs = A["lowpass_fc"], calcjacobian = true)
    # Install the observation. A fresh operator has ZERO-LENGTH data fields --
    # makenoisydata! is what sizes them -- so these are field ASSIGNMENTS, not
    # broadcasts, and W has to be rebuilt by hand. W is Diagonal(1/sigma) over
    # the SELECTED gates; leaving it stale would weight the new data by the old
    # noise and every chi2 downstream would be quietly wrong.
    aem.dlow = copy(A["d_low"]);        aem.dhigh = copy(A["d_high"])
    aem.σlow = copy(A["sigma_low"]);    aem.σhigh = copy(A["sigma_high"])
    aem.selectlow = copy(A["select_low"]); aem.selecthigh = copy(A["select_high"])
    aem.ndatalow = Int(A["ndata_low"]); aem.ndatahigh = Int(A["ndata_high"])
    sg = vcat(aem.σlow[aem.selectlow], aem.σhigh[aem.selecthigh])
    aem.W = sparse(Diagonal(1 ./ sg))
    aem.res = zeros(length(sg))
    aem.J = zeros(length(sg), length(aem.ρ) - 1)
    aem
end

function main()
    cfg  = SiteConfig.resolve(CFGPATH)
    info = SiteConfig.assertvalid(cfg; root = ROOT)
    d, g, pr, hy, out = cfg["data"], cfg["gate"], cfg["priors"], cfg["hypothesis"], cfg["output"]
    odir = joinpath(ROOT, out["dir"])

    @printf("Stage 4 -- gate.  site %s  (%s)\n", cfg["name"], cfg["label"])
    @printf("  lineage: column %s, parent %s, role %s\n",
            d["column"], isempty(info.parent) ? "-" : info.parent, info.role)

    A  = load(joinpath(odir, "stage0_$(cfg["name"]).jld2"))
    S3 = load(joinpath(odir, "stage3_$(cfg["name"]).jld2"))
    println("  reads: stage0 artefact + stage3 window. NO truth file is opened.")

    aem = operator_from_stage0(A)
    zall = A["zall"]; n = length(aem.ρ) - 1; ndata = A["ndata"]
    dzc = diff(vcat(0.0, collect(Float64, zall)))
    R = makeR1(n); m0 = fill(-2.0, n); RtR = Matrix(R'R)
    BETA2 = Float64(g["bayes"]["beta2"])
    Pof(l2) = (M = l2 .* RtR; M[diagind(M)] .+= l2 * BETA2; Symmetric(M))

    # ---- the window, from STAGE 3 -------------------------------------
    wtop = Float64(S3["fixed_window"].top); wbot = Float64(S3["fixed_window"].base)
    celltop = vcat(0.0, collect(Float64, zall)[1:end-1]); cellbot = collect(Float64, zall)
    ov = [max(0.0, min(cellbot[i], wbot) - max(celltop[i], wtop)) for i in 1:n]
    fixidx = findall(>(1e-9), ov)
    dzeff = copy(dzc); dzeff[fixidx] .= ov[fixidx]      # overlap-truncated
    thick = sum(dzeff[fixidx])
    @printf("\n  Stage 3 window %.2f-%.2f m (%.2f m)\n", wtop, wbot, wbot - wtop)
    @printf("  constraint acts on %d of %d cells; cells span %.2f-%.2f m\n",
            length(fixidx), n, celltop[fixidx[1]], cellbot[fixidx[end]])
    @printf("  EFFECTIVE support (overlap-truncated) = %.6f m vs window %.6f m  [gap %.2e m]\n",
            thick, wbot - wtop, abs(thick - (wbot - wtop)))
    abs(thick - (wbot - wtop)) < 1e-9 ||
        error("overlap truncation failed: forced support $thick != window $(wbot-wtop)")
    println("  forced quantity: CONDUCTANCE  sum(10^m_i * dz_eff_i) = C  [ctype=:cond]")

    # ---- reference fit -----------------------------------------------
    mref, chiref = occam(aem, m0, m0, R, ndata; mode = :minimise)
    @printf("  :minimise reference fit: chi2 %.3f (/N %.3f)\n", chiref, chiref/ndata)
    flush(stdout)

    # ---- ML-II for lambda-hat ----------------------------------------
    local lhat, ladder_lhat, unimodal, mls, mlok, lamgrid
    lamgrid = Float64.(g["bayes"]["lambda2_ladder"])
    # With --lhat given, the ladder is SKIPPED, not run and discarded: the brief
    # asks tier F to reuse lambda-hat rather than re-select it, and the saving is
    # the ladder's 7 unconstrained solves per realisation.
    if LHAT_OVERRIDE > 0
        @printf("\nML-II ladder SKIPPED: lambda-hat = %.4g reused from tier P\n", LHAT_OVERRIDE)
        lhat = LHAT_OVERRIDE; ladder_lhat = NaN; unimodal = false
        mls = Float64[]; mlok = Bool[]; lamgrid = Float64[]
        @goto after_ladder
    end
    println("\nML-II for lambda2")
    @printf("  %-10s %12s %10s %12s %6s  %s\n", "lambda2", "chi2", "chi2/N", "log ML", "bound", "status")
    mls = Float64[]; mlok = Bool[]
    for l2 in lamgrid
        m, chi, H, nf = occam(aem, mref, m0, R, ndata; mode = :map, lambda2 = l2,
                              beta2 = BETA2, withH = true, info = true)
        P = Pof(l2)
        ld(X) = try; 2*sum(log, diag(cholesky(Symmetric(Matrix(X))).U)); catch; logabsdet(Matrix(X))[1]; end
        v = -chi/2 - ld(H)/2 + ld(P)/2
        ok = nf.nbound == 0 && nf.settled
        push!(mls, v); push!(mlok, ok)
        @printf("  %-10.3g %12.3f %10.3f %12.2f %6d  %s\n", l2, chi, chi/ndata, v, nf.nbound,
                ok ? "usable" : "REJECTED (bound/unsettled)")
        flush(stdout)
    end
    any(mlok) || error("no usable lambda2: every ladder fit was on a bound or unsettled")
    lhat = lamgrid[findmax(i -> mlok[i] ? mls[i] : -Inf, eachindex(mls))[2]]
    jb = findmax(i -> mlok[i] ? mls[i] : -Inf, eachindex(mls))[2]
    unimodal = all(diff(mls[1:jb]) .>= 0) && all(diff(mls[jb:end]) .<= 0)
    @printf("  lambda-hat = %.4g   log ML unimodal: %s\n", lhat, unimodal ? "YES" : "NO")
    ladder_lhat = lhat
    @label after_ladder

    # ---- the conductance grid ------------------------------------------
    sb = Float64.(g["sigma_bounds"])
    # --grid overrides the config, for the launcher's smoke only: a three-node
    # single pass exercises every sweep, both edge searches and every save path
    # in a couple of minutes. Production takes the config's value.
    mode = get(ARG, "grid", get(g, "grid_mode", "single_pass"))

    function logdets(H, P, m)
        a = zeros(n); a[fixidx] .= log(10) .* (10 .^ m[fixidx]) .* dzeff[fixidx]
        Z = nullspace(reshape(a, 1, n))
        ld(X) = (Sy = Symmetric(Z' * Matrix(X) * Z);
                 try; 2*sum(log, diag(cholesky(Sy).U)); catch; logabsdet(Matrix(Sy))[1]; end)
        (ld(H), ld(P))
    end

    function sweep(l2, ltv, Cgv; laplace = true, tag = "")
        NP = length(Cgv)
        chis = fill(NaN, NP); phis = fill(NaN, NP)
        ldH = fill(NaN, NP); ldP = fill(NaN, NP)
        nb = zeros(Int, NP); good = falses(NP); cres = fill(NaN, NP)
        warm = copy(mref); t0 = time()
        for (i, C) in enumerate(Cgv)
            m, chi, H, nf = occam(aem, warm, m0, R, ndata; fixidx = fixidx, fixval = C,
                                  mode = :map, ctype = :cond, dz = dzeff, lambda2 = l2,
                                  beta2 = BETA2, withH = true, info = true)
            if nf.nbound > 0 || nf.conres >= 1e-8 || !nf.settled
                m2, c2, H2, nf2 = occam(aem, mref, m0, R, ndata; fixidx = fixidx, fixval = C,
                                        mode = :map, ctype = :cond, dz = dzeff, lambda2 = l2,
                                        beta2 = BETA2, withH = true, info = true)
                nf2.phi < nf.phi && ((m, chi, H, nf) = (m2, c2, H2, nf2))
            end
            # Independent check of the forced quantity. nf.conres is occam own
            # |gcond(m) - fixval|/fixval, so it is already computed from the
            # returned model; recomputing it HERE, from the caller own dzeff
            # and fixidx, means a future change to gcond, or a dz mismatch
            # between caller and solver, cannot pass silently. 17 Nov: a driver
            # that checked neither reported a chi2 at a conductance 13x off the
            # one it named, with settled = true and a plausible-looking misfit.
            Cach = sum(10.0^m[j] * dzeff[j] for j in fixidx)
            relc = abs(Cach - C) / C
            (relc < 1e-8) == (nf.conres < 1e-8) || error("""
                stage 4 gate: occam constraint residual and an independent
                recomputation disagree at forced conductance C = $C S.
                  occam conres  = $(nf.conres)
                  recomputed    = $relc  (achieved $Cach S)
                The forced quantity is not what it is reported to be. Stop.""")

            chis[i] = chi; phis[i] = nf.phi; nb[i] = nf.nbound; cres[i] = nf.conres
            if laplace
                ldH[i], ldP[i] = logdets(H, Pof(l2), m)
            end
            # RULED 17 Nov: admissibility is the CONSTRAINT and CONVERGENCE, not
            # bound occupancy. A cell on a stated bound is admissible; the bounds
            # are inclusive. nbound is recorded per node and never gates.
            #
            # The retired `nf.nbound == 0` term was not merely reported: it NaN'd
            # nodes out of the band and, in two-pass mode, set the bracket of the
            # fine grid. On 3P it moved the band from 0.0506-0.4320 S to
            # 0.1426-0.1795 S and put the truth at 0.2063 S outside it, while the
            # nodes it dropped sat at chi2 25.3-28.9 against a threshold of 45.84.
            good[i] = nf.conres < 1e-8 && nf.settled
            warm = copy(m)
            print("\r    $(tag)lambda2 $(round(l2,sigdigits=3)): $i/$NP good $(count(good))   "); flush(stdout)
        end
        println()
        @printf("    %.0f s, %d/%d admissible, %d on a bound\n",
                time()-t0, count(good), NP, count(nb .> 0))
        (chi2 = chis, phi = phis, ldH = ldH, ldP = ldP, nbound = nb, good = good, conres = cres)
    end
    local lt, Cg, coarse_step, coarse_info
    if mode == "two_pass"
        NC = Int(g["n_nodes_coarse"]); NF = Int(g["n_nodes_fine"])
        ltc = collect(range(log10(sb[1]), log10(sb[2]), length = NC))
        Cgc = (10 .^ ltc) .* thick
        coarse_step = ltc[2] - ltc[1]
        @printf("\nPASS 1 -- coarse, %d nodes bound to bound (sigma_eff %.3g to %.3g S/m, %.3f dex apart)\n",
                NC, sb[1], sb[2], coarse_step)
        println("  locating the REACHABLE range from bound occupancy, at lambda-hat")
        SC = sweep(lhat, ltc, Cgc; laplace = false, tag = "coarse ")
        for i in 1:NC
            @printf("    node %d  sigma_eff %-10.4g chi2 %-11.3f nbound %-3d conres %-9.2e %s\n",
                    i, 10^ltc[i], SC.chi2[i], SC.nbound[i], SC.conres[i],
                    SC.good[i] ? "REACHABLE" : "unreachable")
        end
        reach = findall(SC.good)
        length(reach) >= 2 || error("""
            the coarse pass found $(length(reach)) reachable node(s) of $NC. The fine
            grid cannot be bracketed. Either this column's conductance lies outside
            the model bounds, or the coarse grid is too sparse to land inside the
            reachable range -- raise [gate].n_nodes_coarse and re-run.""")
        lo_i, hi_i = first(reach), last(reach)
        lt = collect(range(ltc[lo_i], ltc[hi_i], length = NF))
        Cg = (10 .^ lt) .* thick
        coarse_info = (ltc = ltc, chi2 = SC.chi2, nbound = SC.nbound,
                       good = SC.good, reach_lo = ltc[lo_i], reach_hi = ltc[hi_i])
        @printf("\n  reachable range: sigma_eff %.4g to %.4g S/m (%.3f dex), coarse nodes %d-%d of %d\n",
                10^ltc[lo_i], 10^ltc[hi_i], ltc[hi_i]-ltc[lo_i], lo_i, hi_i, NC)
        @printf("PASS 2 -- fine, %d nodes INSIDE it, %.4f dex apart (single-pass was %.3f dex)\n",
                NF, lt[2]-lt[1], (log10(sb[2])-log10(sb[1]))/(Int(g["n_nodes"])-1))
    else
        NPTS = NPTS_OVERRIDE > 0 ? NPTS_OVERRIDE : Int(g["n_nodes"])
        lt = collect(range(log10(sb[1]), log10(sb[2]), length = NPTS))
        Cg = (10 .^ lt) .* thick
        coarse_step = NaN; coarse_info = nothing
        @printf("\n  SINGLE-PASS grid: sigma_eff %.3g to %.3g S/m (MODEL BOUNDS), %d nodes\n",
                sb[1], sb[2], NPTS)
    end
    @printf("  C from %.4g to %.4g S; inadmissible nodes are EXCLUDED, not normalised over\n",
            first(Cg), last(Cg))

    println("\nsweeping conductance")
    SW = Dict{Float64,Any}()
    for l2 in (lhat/10, lhat, lhat*10); SW[l2] = sweep(l2, lt, Cg); end

    # ---- band edges by TRUE BISECTION ----------------------------------
    # Each step solves the constrained problem at a new C rather than
    # interpolating between grid nodes. That is what lets the fine grid drop to
    # 15 nodes: the grid finds the shape of the curve, the bisection finds where
    # it crosses N + delta_chi2. An edge that bisection cannot bracket inside
    # the fine grid is OPEN, and under the two-pass grid it is open by at least
    # one coarse step -- which is the error tier F exists to catch.
    l2prof = minimum(collect(keys(SW)))
    target = ndata + Float64(g["band"]["delta_chi2"])
    etol = Float64(g["band"]["edge_tol"]); emax = Int(g["band"]["edge_max_iter"])
    function chi2_at(C, l2)
        m, chi, _H, nf = occam(aem, mref, m0, R, ndata; fixidx = fixidx, fixval = C,
                               mode = :map, ctype = :cond, dz = dzeff, lambda2 = l2,
                               beta2 = BETA2, withH = false, info = true)
        ok = nf.conres < 1e-8 && nf.settled      # see the ruling above
        (chi, ok)
    end
    # THE SEARCH IS THE VALIDATED MODULE, not a second implementation.
    #
    # What stood here was Stage 4's own bisection, and it carried both of the
    # bugs gate_check/edgefinder.jl was written to prevent: it moved the bracket
    # from the midpoint's status alone and never re-established it from the
    # solved set, and it had no chi2 convergence test at all -- it ran to a
    # tolerance on position or to an iteration cap. EdgeFinder re-establishes
    # the bracket every refinement, asserts it valid before computing an edge,
    # requires a node to settle against ITSELF rather than against whichever
    # node currently bounds the bracket, and is unit-tested against an injected
    # mock solver in test/test_edgefinder.jl. It is now the only search; the old
    # one is deleted rather than left selectable.
    function edge_via_module(l2, upper::Bool, thr::Float64)
        SWp = SW[l2]
        # the already-solved grid, as EdgeFinder's S: abscissa -> (chi, m)
        S = Dict{Float64,Any}()
        for i in eachindex(lt)
            SWp.good[i] && isfinite(SWp.chi2[i]) &&
                (S[lt[i]] = (chi = SWp.chi2[i], m = mref))   # thr applied by EdgeFinder
        end
        length(S) >= 2 || return (edge = NaN, converged = false, state = "no solved nodes",
                                  nsolves = 0)
        solvefn = (x, from) -> begin
            chi, ok = chi2_at(10^x * thick, l2)
            (chi = ok ? chi : Inf, m = from)
        end
        if upper
            r = EdgeFinder.bisect_edge(S, solvefn, thr;
                                       dex_tol = etol, maxref = emax)
        else
            # EdgeFinder locates an UPPER crossing, so the lower edge is found by
            # reflecting the abscissa: a decreasing crossing becomes increasing.
            Sr = Dict(-k => v for (k, v) in S)
            solvefnr = (x, from) -> solvefn(-x, from)
            rr = EdgeFinder.bisect_edge(Sr, solvefnr, thr;
                                        dex_tol = etol, maxref = emax)
            r = (edge = isfinite(rr.edge) ? -rr.edge : NaN, converged = rr.converged,
                 state = rr.state, nsolves = rr.nsolves)
        end
        (edge = r.edge, converged = r.converged, state = r.state, nsolves = r.nsolves)
    end
    @printf("\nband edges by bisection at lambda2 = %.4g\n", l2prof)
    @printf("  absolute rule: chi2 <= N + %.2f = %.2f ; relative: chiref + %.2f = %.2f\n",
            Float64(g["band"]["delta_chi2"]), target,
            Float64(g["band"]["delta_chi2"]), chiref + Float64(g["band"]["delta_chi2"]))
    # BOTH RULES, from the same solved grid. The absolute rule (N + delta) drives
    # the verdict; the relative rule (chiref + delta) is reported beside it, as
    # agreed, so the two can be compared without a second run.
    target_rel = chiref + Float64(g["band"]["delta_chi2"])
    rlo     = edge_via_module(l2prof, false, target)
    rhi     = edge_via_module(l2prof, true,  target)
    rlo_rel = edge_via_module(l2prof, false, target_rel)
    rhi_rel = edge_via_module(l2prof, true,  target_rel)
    elo, ehi = rlo.edge, rhi.edge
    open_lo, open_hi = !rlo.converged, !rhi.converged
    nlo, nhi = rlo.nsolves, rhi.nsolves
    @printf("  lower edge sigma_eff %.5g S/m  C %.4f S  %s [%s]  (%d solves)\n",
            isfinite(elo) ? 10^elo : NaN, isfinite(elo) ? 10^elo*thick : NaN,
            rlo.converged ? "converged" : "NOT CONVERGED", rlo.state, nlo)
    @printf("  upper edge sigma_eff %.5g S/m  C %.4f S  %s [%s]  (%d solves)\n",
            isfinite(ehi) ? 10^ehi : NaN, isfinite(ehi) ? 10^ehi*thick : NaN,
            rhi.converged ? "converged" : "NOT CONVERGED", rhi.state, nhi)

    # ---- assert the fine grid's end nodes are interior -----------------
    # The fine grid was bracketed by coarse nodes already shown to be
    # interior, so its ends must be too. If they are not, the bracket is wrong
    # and the band computed inside it cannot be trusted.
    if mode == "two_pass"
        SWh = SW[lhat]
        ends_ok = SWh.good[1] && SWh.good[end]
        @printf("\n  fine-grid end nodes interior at lambda-hat: first %s, last %s  -> %s\n",
                SWh.good[1] ? "yes" : "NO", SWh.good[end] ? "yes" : "NO",
                ends_ok ? "ASSERTION PASSES" : "ASSERTION FAILS")
        ends_ok || @warn """The fine grid's end nodes are not interior although they
            sit on coarse nodes that were. The bracket is unreliable for this column;
            treat its band edges as OPEN and widen n_nodes_coarse."""
    end

    jldsave(joinpath(odir, "stage4_$(cfg["name"]).jld2");
        column = d["column"], parent = get(d,"parent",""), role = info.role,
        window = (top = wtop, base = wbot), thick = thick,
        fixidx = fixidx, dzeff = dzeff[fixidx], lt = lt, Cg = Cg,
        lambda_hat = lhat, lamgrid = lamgrid, logML = mls, ml_usable = mlok,
        unimodal = unimodal, ndata = ndata, chiref = chiref,
        sweeps = Dict(string(k) => v for (k,v) in SW),
        m_prior = collect(info.m_prior), sigma_bounds = sb,
        grid_mode = mode, coarse_step_dex = coarse_step,
        edge_lo_dex = elo, edge_hi_dex = ehi, edge_lo_open = open_lo, edge_hi_open = open_hi,
        edge_lo_C = isfinite(elo) ? 10^elo*thick : NaN,
        edge_hi_C = isfinite(ehi) ? 10^ehi*thick : NaN,
        edge_lo_state = rlo.state, edge_hi_state = rhi.state,
        edge_lo_converged = rlo.converged, edge_hi_converged = rhi.converged,
        target_absolute = target, target_relative = target_rel,
        edge_lo_C_rel = isfinite(rlo_rel.edge) ? 10^rlo_rel.edge*thick : NaN,
        edge_hi_C_rel = isfinite(rhi_rel.edge) ? 10^rhi_rel.edge*thick : NaN,
        edge_lo_rel_converged = rlo_rel.converged,
        edge_hi_rel_converged = rhi_rel.converged,
        search = "gate_check/edgefinder.jl EdgeFinder.bisect_edge",
        edge_solves = nlo + nhi, edge_profile_lambda2 = l2prof,
        ladder_lambda_hat = ladder_lhat,
        coarse = coarse_info === nothing ? nothing :
                 Dict("ltc"=>coarse_info.ltc, "chi2"=>coarse_info.chi2,
                      "nbound"=>coarse_info.nbound, "good"=>coarse_info.good,
                      "reach_lo"=>coarse_info.reach_lo, "reach_hi"=>coarse_info.reach_hi))
    println("\nwrote ", relpath(joinpath(odir, "stage4_$(cfg["name"]).jld2"), ROOT))
end
main()
