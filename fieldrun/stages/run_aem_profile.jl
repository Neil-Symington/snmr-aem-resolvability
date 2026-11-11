# Ported from udf_realdata/ into the repository. The include paths and the
# project are the only things changed from the working version: COMMON is
# fieldrun/common/ and the project is the repository root, not ../ek_prerun.
# The control file first: support modules read it at load time.
include(joinpath(@__DIR__, "..", "common", "stage_prelude.jl"))
# run_aem_profile.jl — the admissible conductivity band for each site's target
# interval, by constrained profile likelihood on the AEM data.
#
#   julia --project=. run_aem_profile.jl
#
# This is synthetics/occam_coverage.jl reduced to the part that does not need a
# sampler. For a grid of forced values of the interval's geometric-mean
# conductivity, a constrained Gauss-Newton fit finds the best misfit achievable
# with the interval held there. Where that misfit satisfies
#
#     chi2(forced) <= ndata + 3.84
#
# the data cannot reject that conductivity. The band edges are then bisected off
# the grid, because the grid steps ~0.14 dex while the misfit can jump tens of
# chi2 across one step, so a grid-read edge is located only to within a step.
# Bisection only ever moves the upper edge UP, i.e. it only ever loosens the
# certification -- the same direction as in the report.
#
# The misfit reference is ndata + 3.84, NOT the achieved minimum + 3.84. That
# fix is recorded in the project: referencing to an unconstrained smooth minimum
# that over-fits measures every Dchi2 from an over-fitted floor and excluded the
# truth on the synthetics. Both are printed; the band uses ndata.
#
# NOTHING HERE USES A POSTERIOR OR A PRIOR. That is the point of the gate: every
# objection that can be raised about the sampler lands somewhere else.

using JLD2, Printf, Statistics, DelimitedFiles

include(joinpath(COMMON, "aem_lib.jl"))


const NPTS   = 16
const PAD    = 20.0      # grid spans sigma_unc / PAD to sigma_unc * PAD
const NPASS  = 2
const DCHI   = [1.0, 3.84, 6.63]   # 68 / 95 / 99 % on one parameter
const RCHI   = 1.0       # stop bisecting when the bracket is worth < this chi2
const NBISECT = 14

# --wide: extend the forced grid DOWN to the model bound instead of stopping at
# sigma_unc/PAD, and classify every node three ways instead of two.
#
# WHY. The published band's lower edge is open at 19 of 19 soundings, and that
# was being read as physics. It is not: the grid runs sigma_unc/20 to
# sigma_unc*20, and at 16 of those 19 the grid floor sits 1.2 to 2.8 decades
# ABOVE the model bound of 1e-4 S/m. So "open at the floor" means the sweep
# stopped, not that the data ran out. At the other three (EK18_30, OK18_37,
# bon_mid_4) the grid floor falls BELOW the bound, so the lowest nodes could
# only be reached with cells clamped -- the defect the synthetics gate was just
# fixed for.
#
# THREE STATES, not two. A node whose model has cells on the bound has a chi2
# that is an UPPER bound on the true constrained minimum, because the optimiser
# could not go where it wanted. Calling such a node "rejected" is therefore
# unsafe in the direction that matters: it narrows the band and makes
# certification too easy. It is not admissible either. It is INDETERMINATE, and
# a band edge that abuts an indeterminate region is not a chi2 crossing.
const WIDE = any(==("--wide"), ARGS)
const SUF  = WIDE ? "_wide" : ""
const BD   = Bounds()                 # the field default, (-4.0, 0.5) in log10 sigma
const BTOL = 1e-9
nbound(m) = count(x -> x <= BD.lo + BTOL || x >= BD.hi - BTOL, m)

const FORCE = any(==("--force"), ARGS)
# --repolish: load an existing profile and keep polishing until the passes stop
# recovering misfit, then re-bisect. A profile that is still improving is an
# UPPER bound on chi2_min at the points that stalled, which overstates Dchi2,
# narrows the admissible band and therefore makes certification too EASY --
# the unsafe direction. The initial run uses NPASS passes; this drives it to
# convergence without redoing the sweep.
const REPOLISH = any(==("--repolish"), ARGS)
const MAXPASS = let i = findfirst(==("--maxpass"), ARGS)
    i === nothing ? 8 : parse(Int, ARGS[i+1])
end

function profile_site(field_id; grid = build_grid())
    # Idempotent: the profile is the expensive stage, so an existing result is
    # reused unless --force. That lets the chain be re-entered after a single
    # site is re-run without redoing the other four.
    done = joinpath(OUT, "aem_profile$(SUF)_$(field_id).jld2")
    if isfile(done) && !FORCE && !REPOLISH
        P = load(done)
        @printf("%s: profile already present (%.1f-%.1f m), skipping; --force to redo\n",
                field_id, P["win"][1], P["win"][2])
        return (field_id = field_id, win = (P["win"][1], P["win"][2]),
                σ_unc = P["sigma_unc"], lo = P["edge_lo"], hi = P["edge_hi"],
                ndata = P["ndata"], chi2n_unc = P["chi2_unc"] / P["ndata"])
    end
    A = load(joinpath(OUT, "aem_occam_$(field_id).jld2"))
    df = joinpath(OUT, "detect_$(field_id).jld2")
    isfile(df) || (println("  no detection file for $field_id, skipping"); return nothing)
    Dt = load(df)
    za, zb = Dt["win"]
    (isfinite(za) && isfinite(zb)) || (println("  $field_id has no target interval"); return nothing)

    g = load_gates(field_id)
    scale = A["noise_scale"]
    aem = build_operator(g.d_lm, g.σ_lm, g.d_hm, g.σ_hm, A["height_m"];
                         grid = grid, calcjacobian = true, σ_scale = scale)
    n = length(aem.ρ) - grid.nfixed
    ndata = length(aem.res)
    zall = grid.zall
    R = makeR1(n); m0 = fill(-2.0, n)
    fixidx = findall(z -> za <= z <= zb, zall)
    isempty(fixidx) && (println("  window has no AEM cells"); return nothing)

    # The stored model IS this solver's unconstrained minimum on this operator,
    # at the same noise scale, reference and regulariser, so it is refined from
    # itself rather than re-solved from cold starts: occam() exits on the first
    # step that cannot improve, which costs one Jacobian instead of a few
    # hundred forwards. The stored chi2 is checked against the re-evaluated one.
    χstored = chi2(aem, A["m_log10sigma"])
    abs(χstored - A["chi2"]) > 1e-6 * max(1.0, A["chi2"]) &&
        @printf("  NOTE: stored chi2 %.4f, re-evaluated %.4f — operator differs\n",
                A["chi2"], χstored)
    mu, χu = occam_multi(aem, [A["m_log10sigma"]], m0, R, ndata; mode = :minimise)
    σ_unc = 10^mean(mu[fixidx])

    # Node spacing is preserved at the published value so a --wide profile is
    # directly comparable with the one it replaces; only the extent changes.
    const_step = 2 * log10(PAD) / (NPTS - 1)
    lt, npts = if WIDE
        # Cap the TOP at the model ceiling too. sigma_unc * PAD exceeds it at
        # the conductive sites -- UD24's published grid ran to 15.0 S/m against
        # a ceiling of 3.16 -- so those nodes were unreachable in the same way
        # the floor nodes were, and for the same reason. It never mattered at
        # UDF because every band closed well below, but a grid should not
        # contain points the parameterisation cannot represent.
        hi = min(log10(σ_unc * PAD), BD.hi)
        k  = max(NPTS, ceil(Int, (hi - BD.lo) / const_step) + 1)
        (LinRange(BD.lo, hi, k), k)
    else
        (LinRange(log10(σ_unc / PAD), log10(σ_unc * PAD), NPTS), NPTS)
    end
    targets = 10 .^ collect(lt)

    @printf("\n%s  interval %.1f-%.1f m (%d of %d AEM cells, %s)\n", field_id,
            za, zb, length(fixidx), n, Dt["win_source"])
    @printf("  ndata %d, noise x %.3f, unconstrained chi2/N %.4f, interval sigma %.4f S/m\n",
            ndata, scale, χu / ndata, σ_unc)
    @printf("  forced grid %.5f to %.5f S/m, %d points, step %.3f dex%s\n",
            targets[1], targets[end], npts, lt[2] - lt[1],
            WIDE ? "  [--wide: floor is the model bound]" : "")
    flush(stdout)

    χprof = fill(NaN, npts); models = Vector{Vector{Float64}}(undef, npts)
    nbd = zeros(Int, npts)
    i0 = argmin(abs.(lt .- log10(σ_unc)))
    fitval(lv, starts) = occam_multi(aem, starts, m0, R, ndata; fixidx = fixidx,
                                     fixval = lv, mode = :minimise, ctype = :mean)
    t0 = time()
    if REPOLISH && isfile(done)
        P0 = load(done)
        lt = log10.(P0["targets"]); targets = P0["targets"]
        χprof = copy(P0["chi2_prof"]); M0 = P0["models"]
        for i in 1:length(targets); models[i] = M0[:, i]; end
        @printf("  repolishing a stored profile (chi2 min %.2f), up to %d passes\n",
                minimum(filter(isfinite, χprof)), MAXPASS); flush(stdout)
    else
        for i in vcat(i0:npts, reverse(1:i0-1))
            prev = i == i0 ? mu : (i > i0 ? models[i-1] : models[i+1])
            m, χ = fitval(lt[i], [prev, mu])
            χprof[i] = χ; models[i] = m; nbd[i] = nbound(m)
        end
    end
    # Convergence is judged on the points that can MOVE THE BAND, not on the
    # whole grid. The far, strongly-rejected points sit at chi2 in the
    # thousands and keep yielding large absolute improvements that are
    # irrelevant: recovering 556 chi2 at a point whose chi2 is 291,714 cannot
    # change an edge referenced to ndata + 3.84. "Relevant" is anything within
    # a factor of ten of the level.
    relevant(χv) = findall(c -> isfinite(c) && c <= 10 * (ndata + DCHI[2]), χv)
    for pass in 1:(REPOLISH ? MAXPASS : NPASS)
        gained = 0.0; gained_rel = 0.0
        rel = Set(relevant(χprof))
        ibest = argmin(χprof)
        for i in 1:npts
            starts = Vector{Vector{Float64}}([models[i]])
            i > 1 && push!(starts, models[i-1])
            i < npts && push!(starts, models[i+1])
            i != ibest && push!(starts, models[ibest])
            m, χ = fitval(lt[i], starts)
            if χ < χprof[i] - 1e-6
                d = χprof[i] - χ
                gained += d
                (i in rel) && (gained_rel += d)
                χprof[i] = χ; models[i] = m; nbd[i] = nbound(m)
            end
        end
        @printf("  polish pass %d recovered %.3f chi2 (%.3f near the level)\n",
                pass, gained, gained_rel); flush(stdout)
        gained_rel < 0.1 && break
    end
    @printf("  profile done in %.0f s\n", time() - t0); flush(stdout)

    # A constrained minimum cannot beat the unconstrained one; if it does, the
    # unconstrained run stalled and is not the reference.
    χref = min(χu, minimum(χprof))
    χref < χu - 1e-9 &&
        @printf("  NOTE: constrained sweep beat the unconstrained fit by %.3f chi2; re-referencing\n",
                χu - χref)

    rlevel = ndata + DCHI[2]
    # THREE STATES. indet = the fit put cells on the model bound, so its chi2 is
    # an upper bound and admissibility cannot be decided; adm = admissible and
    # clean; the rest are cleanly rejected.
    indet = findall(nbd .> 0)
    adm   = findall(i -> isfinite(χprof[i]) && χprof[i] <= rlevel && nbd[i] == 0, 1:npts)
    if !isempty(indet)
        @printf("  %d of %d nodes INDETERMINATE (cells on the model bound), sigma %.5g to %.5g S/m\n",
                length(indet), npts, targets[minimum(indet)], targets[maximum(indet)])
        nadm_ind = count(i -> isfinite(χprof[i]) && χprof[i] <= rlevel, indet)
        nadm_ind > 0 && @printf("  of those, %d would have passed the chi2 level had they been trusted\n", nadm_ind)
    end
    edge_lo = edge_hi = NaN; lo_open = hi_open = false
    lo_dchi = hi_dchi = NaN
    # initialised here, not in the else-branch below: if nothing is admissible
    # the branch never runs and jldsave would hit an undefined name
    lo_reason = hi_reason = "no admissible node on the grid"
    if isempty(adm)
        @printf("  NOTHING admissible at chi2 <= %.2f (min %.2f) — the interval cannot be\n",
                rlevel, minimum(χprof))
        @printf("  fixed at ANY conductivity on this grid without exceeding the level\n")
    else
        ilo, ihi = minimum(adm), maximum(adm)
        function bisect(i_in, i_out)
            lin, lout = lt[i_in], lt[i_out]
            m_in, m_out = models[i_in], models[i_out]
            χin, χout = χprof[i_in], χprof[i_out]
            for _ in 1:NBISECT
                abs(χout - χin) < RCHI && break
                lmid = (lin + lout) / 2
                m, χ = fitval(lmid, [m_in, m_out, mu])
                if χ <= rlevel
                    lin, m_in, χin = lmid, m, χ
                else
                    lout, m_out, χout = lmid, m, χ
                end
                abs(lout - lin) < 1e-4 && break
            end
            10^lin, abs(χout - χin)
        end
        # An edge is only a chi2 crossing if the node OUTSIDE it is cleanly
        # rejected. If it is indeterminate, or there is no node outside at all,
        # the edge is a limit of the sweep or of the parameterisation and must
        # be reported as open rather than bisected to a spurious precision.
        lo_reason = hi_reason = "chi2 crossing"
        if ilo == 1
            edge_lo = targets[1]; lo_open = true
            lo_reason = WIDE ? "grid floor = model bound" : "grid floor (sweep limit)"
        elseif nbd[ilo - 1] > 0
            edge_lo = targets[ilo]; lo_open = true
            lo_reason = "node below is indeterminate (on the model bound)"
        else
            edge_lo, lo_dchi = bisect(ilo, ilo - 1)
        end
        if ihi == npts
            edge_hi = targets[npts]; hi_open = true
            hi_reason = "grid ceiling (sweep limit)"
        elseif nbd[ihi + 1] > 0
            edge_hi = targets[ihi]; hi_open = true
            hi_reason = "node above is indeterminate (on the model bound)"
        else
            edge_hi, hi_dchi = bisect(ihi, ihi + 1)
        end
        @printf("  lower edge: %s\n  upper edge: %s\n", lo_reason, hi_reason)
        # headroom to the parameterisation, which the field table should carry
        @printf("  headroom: lower edge %.2f dex above the model floor, upper edge %.2f dex below the ceiling\n",
                log10(edge_lo) - BD.lo, BD.hi - log10(edge_hi))
        @printf("  95%% ADMISSIBLE BAND [%.5f, %.5f] S/m%s%s\n", edge_lo, edge_hi,
                lo_open ? " [OPEN low]" : "", hi_open ? " [OPEN high]" : "")
        @printf("     = [%.0f, %.0f] ohm.m, %.2f decades wide; unconstrained fit at %.4f\n",
                1 / edge_hi, 1 / edge_lo, log10(edge_hi / edge_lo), σ_unc)
    end
    for (lbl, Δ) in zip(("68%", "95%", "99%"), DCHI)
        ok = findall(i -> isfinite(χprof[i]) && χprof[i] <= ndata + Δ && nbd[i] == 0, 1:npts)
        isempty(ok) ? @printf("    %-4s nothing on the grid\n", lbl) :
            @printf("    %-4s [%.5f, %.5f] S/m\n", lbl, targets[minimum(ok)], targets[maximum(ok)])
    end
    flush(stdout)

    jldsave(joinpath(OUT, "aem_profile$(SUF)_$(field_id).jld2");
        field_id = field_id, win = [za, zb], win_source = Dt["win_source"],
        targets = targets, chi2_prof = χprof, chi2_ref = χref, chi2_unc = χu,
        ndata = ndata, noise_scale = scale, rlevel = rlevel, thresholds = DCHI,
        sigma_unc = σ_unc, zall = zall, fixidx = fixidx,
        models = hcat(models...), m_unconstrained = mu,
        edge_lo = edge_lo, edge_hi = edge_hi,
        edge_lo_open = lo_open, edge_hi_open = hi_open,
        edge_lo_dchi = lo_dchi, edge_hi_dchi = hi_dchi,
        npts = npts, pad = PAD, npass = NPASS,
        nbound = nbd, wide = WIDE, bound_lo = BD.lo, bound_hi = BD.hi,
        edge_lo_reason = lo_reason, edge_hi_reason = hi_reason)
    (field_id = field_id, win = (za, zb), σ_unc = σ_unc,
     lo = edge_lo, hi = edge_hi, ndata = ndata, chi2n_unc = χu / ndata)
end

function main()
    d = readdlm(joinpath(FieldConfig.cfgpath(FieldConfig.CFG[], "aem_sounds"), "manifest.csv"), ',', String; header = true)
    ci = Dict(strip(h) => i for (i, h) in enumerate(vec(d[2])))
    rows = []
    nattempt = 0
    for s in d[1][:, ci["field_id"]]
        (!isempty(ONLY) && s != ONLY) && continue
        nattempt += 1
        r = try profile_site(s) catch e
            @warn "profile failed for $s" exception = (e, catch_backtrace()); nothing
        end
        r !== nothing && push!(rows, r)
    end
    println("\n", "="^92)
    @printf("%-14s %-16s %9s %11s %11s %7s\n", "site", "interval m", "chi2/N",
            "sigma_lo", "sigma_hi", "decades")
    for r in rows
        @printf("%-14s %6.1f-%-9.1f %9.3f %11.5f %11.5f %7.2f\n",
                r.field_id, r.win[1], r.win[2], r.chi2n_unc, r.lo, r.hi,
                log10(r.hi / r.lo))
    end
    println("="^92)
    # A stage that produced NOTHING is a failed stage, not an empty table.
    # All three sweeps printed a bare header and exited 0 when every site was
    # skipped for missing inputs, so run_batch.sh reported "ok" for a sweep
    # that did no work at all. Say so, and say it in the exit code.
    if isempty(rows)
        @printf("\nFAILED: no site produced a result (%d attempted).\n", nattempt)
        println("A skipped site above means its input file is missing; check that")
        println("the per-site stages actually wrote their outputs.")
        exit(1)
    end
    nskip = nattempt - length(rows)
    nskip > 0 && @printf("NOTE: %d of %d sites skipped or failed.\n", nskip, nattempt)
end

main()
