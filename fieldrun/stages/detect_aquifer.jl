# Ported from udf_realdata/ into the repository. The include paths and the
# project are the only things changed from the working version: COMMON is
# fieldrun/common/ and the project is the repository root, not ../ek_prerun.
# detect_aquifer.jl — "is there an aquifer?", by the paper's Stage 3 criteria,
# applied to the Upper Darling Floodplain SNMR posteriors.
#
# Identical to ek_realdata/detect_aquifer.jl except for `kernel_for`, which
# rebuilds the SMR kernel from the UDF job assembly rather than the East
# Kimberley inversion spreadsheet, and uses the SAME AEM conductivity draw the
# inversion was sampled under rather than the smooth model.
#
#   julia --project=. detect_aquifer.jl
#
# The criteria are transcribed from the synthetics, not reinvented:
#
#   * constraint depth zlim          stage1_target.jl constraint_depth()
#                                    the deepest cell whose minimum detectable
#                                    water-content change is <= 0.05 v/v
#   * extractable unit, per draw     stage1_target.jl allruns()
#                                    contiguous, theta > 0.15, >= 10 m, above zlim
#   * consensus interval             consensus_unit.jl
#                                    depths inside a unit in >= 50% of draws
#   * W = int theta dz               stage1_integral.jl
#                                    over the consensus interval, with P(W >= w)
#   * E = max_interval int (theta - floor) dz   stage1_integral.jl
#                                    maximum-subarray, floor swept
#   * best L-metre screen            interval_from_integral.jl
#                                    no theta threshold anywhere; L = 12 m headline
#
# Every statistic is computed per draw and then summarised. Nothing is computed
# from q50: the integral of a percentile profile is not a percentile of the
# integral, and the report's own numbers were corrected for exactly that error.
#
# The sensitivity DOI (doi_study/doi_sensitivity.jl) is computed alongside and
# reported, because the posterior-width DOI that the pre-run table carries is
# partly a readout of model complexity and should not be quoted on its own.

using JLD2, Statistics, Printf, LinearAlgebra, DelimitedFiles

# The control file first: support modules read it at load time.
include(joinpath(@__DIR__, "..", "common", "stage_prelude.jl"))

# PRERUN is gone: prerun_data.jl and prerun_lib.jl are in common/
include(joinpath(COMMON, "prerun_data.jl"))
include(joinpath(COMMON, "udf_smr.jl"))

const HERE = @__DIR__
const OUT  = joinpath(HERE, "outputs")

const THETA_MIN   = 0.15      # stage1_target.jl default
const THICK_MIN   = 10.0      # stage1_target.jl default
const THETA_DETECT = 0.05     # stage1_target.jl detectability floor
const FLOORS      = (0.05, 0.10, 0.15, 0.20)
const FLOOR_MAIN  = 0.15
const SCREENS     = (8.0, 10.0, 12.0, 16.0, 20.0)
const SCREEN_MAIN = 10.0   # the paper uses L = 10 m
const WGRID       = collect(0.0:0.05:8.0)
# a screen whose base sits within this fraction of the constraint depth is
# PINNED: resting on the limit rather than located by the data
const PIN_FRAC    = 0.90
# A depth is RESPONDING while the posterior there is narrower than this
# fraction of the prior width. Paper Stage 2 asks for "the depth below which the
# profile stops responding to the data, and the integrated water bound below
# that depth" -- this is that depth, measured on the posterior itself.
#
# It is needed because the sensitivity-based constraint depth structurally
# OVERSTATES resolution: it asks whether a change at one depth alone would be
# detectable with everything else known, which conditions away the trade-offs
# that actually limit a layered-earth inversion (the paper's own objection to
# single-depth DOI, section 2.5). At bon_mid_4 the detectability limit said
# 115.4 m and the sensitivity DOI 123.9 m on an 80 m loop, while the posterior
# retained 65% of the prior p5-p95 width at 99.7 m -- and a spurious "aquifer"
# was selected there.
#
# Measured as (posterior p95 - p5) / (prior p95 - p5), and read as the TERMINAL
# run of uninformed cells at the base of the column; see the comment at its use
# for why the first crossing is the wrong rule.
const RESPOND_FRAC = 0.50

"Cell tops, bases and thicknesses from node depths (stage1_integral.jl cellsof)."
function cellsof(z)
    mid = (z[1:end-1] .+ z[2:end]) ./ 2
    top = vcat(0.0, mid)
    bot = vcat(mid, z[end] + (z[end] - z[end-1]) / 2)
    (top = top, bot = bot, dz = bot .- top)
end

"Every contiguous run above zlim with theta > tmin and thickness >= Lmin."
function allruns(θ, g, zlim; θmin = THETA_MIN, tmin = THICK_MIN)
    out = NTuple{4,Float64}[]
    run = 0.0; cur0 = 0.0; acc = 0.0
    for i in eachindex(θ)
        if g.bot[i] <= zlim && θ[i] > θmin
            run == 0 && (cur0 = g.top[i]; acc = 0.0)
            run += g.dz[i]; acc += θ[i] * g.dz[i]
        else
            run >= tmin && push!(out, (cur0, g.top[i], run, acc / run))
            run = 0.0
        end
    end
    run >= tmin && push!(out, (cur0, g.bot[end], run, acc / run))
    out
end

"Contiguous index runs with p >= thr."
function runs_above(p, thr)
    out = UnitRange{Int}[]; i = 1
    while i <= length(p)
        if p[i] >= thr
            j = i; while j < length(p) && p[j+1] >= thr; j += 1; end
            push!(out, i:j); i = j + 1
        else
            i += 1
        end
    end
    out
end

"Maximum-subarray: the contiguous interval holding the most water above a floor."
function maxsub(x)
    best = -Inf; bi0 = bi1 = 1; cur = 0.0; c0 = 1
    for i in eachindex(x)
        if cur <= 0
            cur = x[i]; c0 = i
        else
            cur += x[i]
        end
        if cur > best
            best = cur; bi0 = c0; bi1 = i
        end
    end
    (best, bi0, bi1)
end

"Index ranges spanning about L metres and ending above zlim."
function windows(g, zlim, L)
    n = length(g.dz); out = UnitRange{Int}[]
    for i in 1:n
        g.top[i] >= zlim && break
        j = i
        while j < n && g.bot[j] < zlim && (g.bot[j] - g.top[i]) < L
            j += 1
        end
        (g.bot[j] <= zlim && (g.bot[j] - g.top[i]) >= 0.6L) && push!(out, i:j)
    end
    out
end

"""
Minimum detectable water-content change per cell, and the constraint depth.

`K̃` is the complex amplitude kernel (nq x nz) with cell thickness already
folded in, so V0 = |K̃ w|. A change dtheta in cell i alone moves datum q by
|K̃[q,i]| dtheta, so the change detectable at one sigma overall is
dtheta_min[i] = 1 / sqrt(sum_q (|K̃[q,i]|/sigma_q)^2).
"""
function detectability(K̃, σ, z; θdetect = THETA_DETECT)
    nz = size(K̃, 2)
    dmin = [1 / sqrt(sum(abs2(abs(K̃[q, i]) / σ[q]) for q in axes(K̃, 1))) for i in 1:nz]
    i = findlast(<=(θdetect), dmin)
    (dmin = dmin, zlim = i === nothing ? 0.0 : Float64(z[i]))
end

"""
Cumulative noise-normalised sensitivity from the bottom up, and its DOI.

doi_study/doi_sensitivity.jl, in full: C(j) is how many sigma of data change
the whole section at or below cell j could produce for a dtheta perturbation.
"""
function sens_doi(K̃, σ, z, w_ref; Δw = 0.05)
    nz = size(K̃, 2)
    w = fill(w_ref, nz)
    zc = K̃ * w
    u = conj.(zc) ./ max.(abs.(zc), eps())
    J = real.(u .* K̃)
    s2 = vec(sum((J .* Δw ./ σ) .^ 2, dims = 1))
    C = sqrt.(reverse(cumsum(reverse(s2))))
    idx = findfirst(<(1.0), C)
    (C = C, doi = idx === nothing ? Float64(z[end]) : (idx == 1 ? 0.0 : Float64(z[idx-1])))
end

"""
Rebuild the SMR kernel for a site from the stored sounding inputs.

The kernel must use the SAME conductivity the inversion was sampled under.
Each UDF ensemble drew its own Stage-1 AEM profile, so the stored
`aem_sigma`/`aem_ztop` (ensemble 1's draw) is used rather than the smooth
Occam model. That matters much more here than at East Kimberley: pointwise AEM
conductivity is uncertain by 4-64x over this floodplain, so rebuilding on a
different earth would evaluate detection against a kernel the posterior never
saw.
"""
function kernel_for(D; ens::Int = 1)
    fid = D["field_id"]
    r = udf_meta(fid)
    b = build_udf_sounding(merge(r, (matpath = udf_matfile(fid)[1],)))
    b.ok || error("$fid: $(b.why)")
    td = joinpath(OUT, "aem_td_$(fid).jld2")
    cm = if ens > 1 && isfile(td)
        A = load(td); Dr = A["sigma_draws"]
        SNMRForward.ConductivityModel(Dr[:, mod1(ens, size(Dr, 2))], diff(A["zbnd"]))
    elseif haskey(D, "aem_sigma")
        SNMRForward.ConductivityModel(D["aem_sigma"], diff(D["aem_ztop"]))
    else
        A = load(joinpath(OUT, "aem_occam_$(fid).jld2"))
        SNMRForward.ConductivityModel(A["sigma"], diff(A["zbnd"]))
    end
    F = SMRPInversion.MRSForward_square(r.loop_side / 2, b.zbnd, b.q,
            deg2rad(r.inc_deg), 0.0, b.Be, cm)
    SMRPInversion.SMRDetOperator(sounding_fixed(merge(b, (F = F,)))).K̃
end

function analyse(field_id)
    f = joinpath(OUT, "smr_transd_$(field_id).jld2")
    isfile(f) || (println("  no trans-D file for $field_id, skipping"); return nothing)
    D = load(f)
    Θ = D["theta"]; z = D["zbnd"]; g = cellsof(z); nd = size(Θ, 2)
    σ̂ = D["sigma_hat"]
    s50 = quantile(exp.(D["log_s"]), 0.5)
    σ_eff = sqrt(max(s50, 0.0)) .* σ̂

    K̃ = kernel_for(D)
    # Posterior predictive, EVALUATED WITHIN KERNEL.
    #
    # Each ensemble sampled under its own Stage-1 AEM draw, so there is no
    # single kernel the pooled posterior fits. Scoring every draw against one
    # of them gives chi2/N 2.5-39 here and says nothing: measured per ensemble
    # against its own kernel the same draws sit at 0.98-1.05, while ensemble 2-4
    # draws scored against ensemble 1's kernel reach 5.8-81. Over this
    # floodplain pointwise AEM conductivity is uncertain by 4-64x, so the
    # kernels genuinely differ.
    #
    # This is the third pooled diagnostic that per-draw kernels break, after
    # R-hat and the split-half drift. Any statistic that compares draws to data
    # must be computed within kernel.
    npp = min(400, nd)
    nens = Int(get(D, "nensembles", 1)); per = nd ÷ max(nens, 1)
    Ks = Vector{Any}(undef, nens)
    for e in 1:nens
        Ks[e] = e == 1 ? K̃ : kernel_for(D; ens = e)
    end
    kof(j) = Ks[clamp(fld(j - 1, per) + 1, 1, nens)]
    jj = round.(Int, range(1, nd, length = npp))
    pp = hcat([abs.(kof(j) * Θ[:, j]) for j in jj]...)
    ppq = [quantile(view(pp, i, :), p) for i in axes(pp, 1), p in (0.05, 0.5, 0.95)]
    chi2n_draws = [mean(((D["V0"] .- abs.(kof(j) * Θ[:, j])) ./
                         (sqrt(exp(D["log_s"][j])) .* σ̂)) .^ 2) for j in jj]
    det = detectability(K̃, σ_eff, z)
    # the limit is the SHALLOWER of the sensitivity-based depth and the depth at
    # which the posterior stops responding
    #
    # Measure "responding" as the posterior p5-p95 width as a fraction of the
    # PRIOR's own p5-p95 width (0.9 of the uniform range, not the full range):
    # a cell the data say nothing about sits at wfrac ~ 1, not ~ 0.9.
    #
    # And take the TERMINAL run, not the first crossing. A first-crossing rule
    # cuts the column at the first wide cell, and the widest cell in a
    # well-resolved sounding is the aquifer itself -- theta is large there, so
    # its absolute uncertainty is large too. At bon_mid_4 first-crossing put the
    # limit at 51.7 m, i.e. through the middle of the water the bore confirms.
    # Loss of resolution is a property of the bottom of the column: scan up from
    # the base while the posterior is still uninformed and stop at the first
    # cell that is not.
    pw = 0.9 * (D["fbounds"][2] - D["fbounds"][1])
    wfrac = (D["q95"] .- D["q05"]) ./ pw
    it = length(z) + 1
    while it - 1 >= 1 && wfrac[it-1] > RESPOND_FRAC; it -= 1; end
    z_respond = it > length(z) ? Float64(z[end]) : (it == 1 ? 0.0 : Float64(z[it-1]))
    zlim = min(det.zlim, z_respond)
    # The synthetics compute the constraint depth at s = 1, i.e. from sigma-hat
    # alone. On real records sigma-hat is decay-model mismatch and understates
    # the error, so the honest constraint depth uses the INFERRED noise. Both
    # are reported, because the difference is large and it is the reason the
    # criteria see less column here than they do on a synthetic.
    det_hat = detectability(K̃, σ̂, z)
    sd = [sens_doi(K̃, σ_eff, z, wr).doi for wr in (0.02, 0.05, 0.10)]

    println("="^100)
    @printf("%s — %d draws, %d pulse moments, sigma_eff median %.2f nV (s p50 %.2f)\n",
            field_id, nd, length(σ̂), median(σ_eff) * 1e9, s50)
    @printf("  depth limit %.1f m = min(detectability %.1f m, posterior-responding %.1f m)\n",
            zlim, det.zlim, z_respond)
    @printf("  constraint depth (detectable dtheta <= %.2f v/v): %.1f m at the inferred\n",
            THETA_DETECT, det.zlim)
    @printf("      noise, %.1f m at sigma-hat alone (s = 1, the synthetics' convention)\n",
            det_hat.zlim)
    @printf("  sensitivity DOI at w_ref 0.02/0.05/0.10: %.1f / %.1f / %.1f m\n", sd...)
    @printf("  posterior-width DOI (NOT quoted alone): %.1f m\n", D["doi_width"])
    @printf("  convergence: rhat_max %.3f, ess_bulk_min %.0f, gates %s\n",
            D["rhat_max"], D["ess_bulk_min"], D["gates_pass"] ? "PASS" : "FAIL")
    @printf("  posterior predictive chi2/N under the inferred noise: p05 %.2f p50 %.2f p95 %.2f\n",
            quantile(chi2n_draws, .05), quantile(chi2n_draws, .5), quantile(chi2n_draws, .95))

    # ── the threshold rule, per draw ────────────────────────────────────────
    units = [allruns(view(Θ, :, j), g, zlim) for j in 1:nd]
    nun = [length(u) for u in units]
    P1 = mean(nun .>= 1); P2 = mean(nun .>= 2)
    occ = zeros(length(z))
    for u in units, (a, b, _, _) in u, i in eachindex(z)
        (g.top[i] >= a - 1e-9 && g.bot[i] <= b + 1e-9) && (occ[i] += 1)
    end
    occ ./= nd
    @printf("\n  THRESHOLD RULE (theta > %.2f over >= %.0f m, above %.1f m)\n",
            THETA_MIN, THICK_MIN, zlim)
    @printf("    P(>= 1 extractable unit) %.3f    P(>= 2) %.3f\n", P1, P2)

    cons = runs_above(occ, 0.5)
    cons = [r for r in cons if g.bot[last(r)] - g.top[first(r)] > 2.0]
    have_cons = !isempty(cons)
    if have_cons
        for r in cons
            @printf("    consensus interval (>= 50%% of draws): %.1f-%.1f m (%.1f m)\n",
                    g.top[first(r)], g.bot[last(r)], g.bot[last(r)] - g.top[first(r)])
        end
    else
        @printf("    NO consensus interval at 50%% — the threshold rule gives no unit here\n")
    end

    # ── W over a fixed interval ─────────────────────────────────────────────
    # Interval preference: the consensus unit when the threshold rule yields
    # one; otherwise the modal best-screen window, which needs no threshold.
    # ── Stage 3, with a determinedness check ────────────────────────────────
    # "Slide a window and keep the best" ranks on W alone, which gives no
    # defence against a deep, broad, poorly-determined feature outranking a
    # shallow, sharp, well-determined one. At bon_mid_4 that is exactly what
    # happened: the deep candidate won on W (3.42 m vs 2.04 m) while being
    # twice as uncertain (CV 0.155 vs 0.078) and less agreed-upon
    # (P(modal) 0.35 vs 0.56), and the bore 62 m away shows it does not exist.
    #
    # The separator is PINNING, the diagnostic already validated above: the
    # deep candidate's base sat at 113.4 m against a 115.4 m limit. So the
    # best screen is now found twice -- over everything, and over windows whose
    # base clears the limit by 10% -- and both are reported. The UNPINNED one
    # is primary, because a window resting on the limit is a lower bound on a
    # quantity the sounding has stopped resolving, not an estimate.
    ws_main = windows(g, zlim, SCREEN_MAIN)
    "Best screen over a given window set: per-draw W, modal window, diagnostics."
    function pick(ws)
        isempty(ws) && return nothing
        bi = Vector{Int}(undef, nd); bv = Vector{Float64}(undef, nd)
        for j in 1:nd
            bw = -Inf; k0 = 1
            for (k, r) in enumerate(ws)
                v = 0.0
                @inbounds for i in r; v += Θ[i, j] * g.dz[i]; end
                v > bw && (bw = v; k0 = k)
            end
            bi[j] = k0; bv[j] = bw
        end
        cnt = [count(==(k), bi) / nd for k in eachindex(ws)]
        km = argmax(cnt)
        (W = bv, idx = ws[km], top = g.top[first(ws[km])], bot = g.bot[last(ws[km])],
         P = cnt[km], cv = std(bv) / mean(bv), cnt = cnt, which = bi)
    end
    all_pick = pick(ws_main)
    # windows clear of the limit: the base must leave PIN_FRAC of headroom
    ws_free  = [r for r in ws_main if g.bot[last(r)] <= PIN_FRAC * zlim]
    free_pick = pick(ws_free)
    if all_pick !== nothing && free_pick !== nothing
        @printf("\n  STAGE 3 CANDIDATES (L = %.0f m)\n", SCREEN_MAIN)
        @printf("    %-22s %7s %7s %7s %7s %8s\n", "", "top", "base", "W p50", "CV", "P(modal)")
        @printf("    %-22s %7.1f %7.1f %7.2f %7.3f %8.2f\n", "best W, any depth",
                all_pick.top, all_pick.bot, quantile(all_pick.W, .5), all_pick.cv, all_pick.P)
        @printf("    %-22s %7.1f %7.1f %7.2f %7.3f %8.2f\n", "best W, clear of limit",
                free_pick.top, free_pick.bot, quantile(free_pick.W, .5), free_pick.cv, free_pick.P)
        if abs(all_pick.top - free_pick.top) > 1e-6
            if all_pick.bot > PIN_FRAC * zlim
                @printf("    DIFFER, and the best-W screen rests on the constraint depth (%.1f of %.1f m).\n",
                        all_pick.bot, zlim)
                @printf("    Ranking on W alone cannot see this; the unpinned candidate is primary.\n")
            else
                # both candidates clear the limit: restricting the window set has
                # reallocated draws and moved the MODAL window, not the pinning
                @printf("    DIFFER on the modal window only (%.1f m base, limit %.1f m): neither is pinned.\n",
                        all_pick.bot, zlim)
            end
        end
    end
    # the unpinned candidate is primary where one exists
    chosen = free_pick === nothing ? all_pick : free_pick
    best_idx = chosen === nothing ? Int[] : chosen.which
    best_val = chosen === nothing ? Float64[] : chosen.W
    modal = chosen === nothing ? nothing : chosen.idx

    if have_cons
        ix = collect(cons[1]); wsrc = "consensus interval"
    elseif modal !== nothing
        ix = collect(modal); wsrc = "modal best-$(Int(SCREEN_MAIN))m screen"
    else
        ix = Int[]; wsrc = "none"
    end
    za = isempty(ix) ? NaN : g.top[first(ix)]
    zb = isempty(ix) ? NaN : g.bot[last(ix)]
    W = isempty(ix) ? Float64[] : [sum(Θ[ix, j] .* g.dz[ix]) for j in 1:nd]

    if !isempty(W)
        @printf("\n  W = int theta dz over %.1f-%.1f m (%s)\n", za, zb, wsrc)
        @printf("    p05 %.2f  p50 %.2f  p95 %.2f m of water   CV %.3f\n",
                quantile(W, .05), quantile(W, .5), quantile(W, .95), std(W) / mean(W))
        @printf("    P(W >= 1 m) %.3f   P(W >= 2 m) %.3f   P(W >= 3 m) %.3f\n",
                mean(W .>= 1), mean(W .>= 2), mean(W .>= 3))
        cvcell = median([std(Θ[i, :]) / mean(Θ[i, :]) for i in ix])
        @printf("    single-cell CV over the same interval: %.3f\n", cvcell)
    end

    # ── E, the threshold-free scan ──────────────────────────────────────────
    above = findall(i -> g.bot[i] <= zlim, 1:length(z))
    function excess(θc)
        e = Vector{Float64}(undef, nd); a = similar(e); b = similar(e)
        for j in 1:nd
            x = (Θ[above, j] .- θc) .* g.dz[above]
            v, i0, i1 = maxsub(x)
            e[j] = v; a[j] = g.top[above[i0]]; b[j] = g.bot[above[i1]]
        end
        (E = e, top = a, bot = b)
    end
    println("\n  E = max over intervals of int (theta - floor) dz, above the constraint depth")
    @printf("    %7s %8s %8s %8s  %-22s %s\n", "floor", "p05", "p50", "p95",
            "p50 interval", "CV")
    Emain = nothing
    for θc in FLOORS
        r = excess(θc)
        @printf("    %7.2f %8.2f %8.2f %8.2f  %7.1f - %-11.1f %.3f\n",
                θc, quantile(r.E, .05), quantile(r.E, .5), quantile(r.E, .95),
                median(r.top), median(r.bot), std(r.E) / mean(r.E))
        θc == FLOOR_MAIN && (Emain = r)
    end

    # ── the screen-length rule ──────────────────────────────────────────────
    println("\n  BEST L-METRE SCREEN (no theta threshold)")
    @printf("    %5s | %-24s %7s | %7s %7s %7s | %s\n", "L (m)", "p50 top-base",
            "modal", "W p05", "W p50", "W p95", "P(W>=2m)")
    screen = Dict{Float64,Any}()
    for L in SCREENS
        ws = windows(g, zlim, L); isempty(ws) && continue
        tops = Float64[]; bots = Float64[]; best = Float64[]; which = Int[]
        for j in 1:nd
            bw = -Inf; bi = 1
            for (k, r) in enumerate(ws)
                v = 0.0
                @inbounds for i in r; v += Θ[i, j] * g.dz[i]; end
                v > bw && (bw = v; bi = k)
            end
            push!(best, bw); push!(which, bi)
            push!(tops, g.top[first(ws[bi])]); push!(bots, g.bot[last(ws[bi])])
        end
        cnt = [count(==(k), which) / nd for k in eachindex(ws)]
        km = argmax(cnt)
        @printf("    %5.0f | %9.1f - %-12.1f %6.2f | %7.2f %7.2f %7.2f | %.3f\n",
                L, median(tops), median(bots), cnt[km],
                quantile(best, .05), quantile(best, .5), quantile(best, .95),
                mean(best .>= 2))
        screen[L] = (best = best, tops = tops, bots = bots, which = which,
                     cnt = cnt, modal_top = g.top[first(ws[km])],
                     modal_bot = g.bot[last(ws[km])])
    end

    if !haskey(screen, SCREEN_MAIN)
        # Can only happen if the constraint depth is shallower than 0.6*L, i.e.
        # the sounding does not constrain even one screen length of section.
        @printf("\n  no %.0f m screen fits above the constraint depth (%.1f m) — ",
                SCREEN_MAIN, zlim)
        println("nothing to report for this site")
        return nothing
    end
    S = screen[SCREEN_MAIN]
    @printf("\n    where to put a %.0f m screen — top candidates:\n", SCREEN_MAIN)
    ws = windows(g, zlim, SCREEN_MAIN)
    for k in sortperm(S.cnt, rev = true)[1:min(4, length(S.cnt))]
        S.cnt[k] < 0.01 && break
        @printf("      %.1f - %.1f m   P %.3f\n",
                g.top[first(ws[k])], g.bot[last(ws[k])], S.cnt[k])
    end

    # ── how much of the answer is the depth limit? ─────────────────────────
    # Every criterion here is defined "above the constraint depth", so if the
    # water sits at or below that depth the rules return a lower bound rather
    # than the aquifer. Rerunning the screen rule at three limits measures how
    # much of the answer is the limit rather than the data.
    nolimit_best = Float64[]; nolimit_top = NaN; nolimit_bot = NaN
    println("\n  DEPTH-LIMIT SENSITIVITY — the best 12 m screen at three limits")
    @printf("    %-34s %8s | %-22s %8s\n", "limit", "depth", "modal screen", "W p50")
    zsens = Any[]
    for (lbl, zl) in (("detectability, inferred noise", zlim),
                      ("detectability, sigma-hat only", det_hat.zlim),
                      ("sensitivity DOI", sd[2]),
                      ("none (full model column)", Float64(g.bot[end])))
        ws = windows(g, zl, SCREEN_MAIN)
        if isempty(ws)
            @printf("    %-34s %8.1f | %-22s %8s\n", lbl, zl, "none fits", "—")
            continue
        end
        bv = Float64[]; wh = Int[]
        for j in 1:nd
            bw = -Inf; bi = 1
            for (k, r) in enumerate(ws)
                v = 0.0
                @inbounds for i in r; v += Θ[i, j] * g.dz[i]; end
                v > bw && (bw = v; bi = k)
            end
            push!(bv, bw); push!(wh, bi)
        end
        cnt = [count(==(k), wh) / nd for k in eachindex(ws)]
        km = argmax(cnt)
        @printf("    %-34s %8.1f | %6.1f - %-6.1f P %.2f %8.2f m\n", lbl, zl,
                g.top[first(ws[km])], g.bot[last(ws[km])], cnt[km], quantile(bv, .5))
        push!(zsens, (label = lbl, zlim = zl, top = g.top[first(ws[km])],
                      bot = g.bot[last(ws[km])], P = cnt[km], W = bv))
        lbl == "none (full model column)" && (nolimit_best = bv;
            nolimit_top = g.top[first(ws[km])]; nolimit_bot = g.bot[last(ws[km])])
    end
    pinned = !isempty(zsens) && (zsens[1].bot >= zlim - 1.5)
    pinned && println("    NOTE: the best screen sits hard against the constraint depth, so the")
    pinned && println("    water is at or below what this sounding resolves — W is a LOWER BOUND.")

    # The verdict is stated AT the constraint depth, and says so, because every
    # statistic behind it is truncated there. Where the best screen is pinned
    # against that depth the water column is a lower bound, so a negative
    # verdict is "not resolved by this sounding", not "not present".
    base = if !isempty(W) && mean(W .>= 2) >= 0.9
        "AQUIFER (W >= 2 m, P >= 0.9)"
    elseif !isempty(W) && mean(W .>= 1) >= 0.9
        "MARGINAL (W >= 1 m, P >= 0.9; not 2 m)"
    else
        "NOT RESOLVED (P < 0.9 for 1 m)"
    end
    verdict = base * @sprintf(" above %.0f m", zlim) *
              (pinned ? "; screen pinned at that depth, W is a LOWER BOUND" : "")
    @printf("\n  VERDICT: %s\n", verdict)

    jldsave(joinpath(OUT, "detect_$(field_id).jld2");
        field_id = field_id, z = z, dz = g.dz, ztop = g.top, zbot = g.bot,
        zlim = zlim, dmin = det.dmin, doi_sens = sd[2],
        doi_sens_lo = sd[1], doi_sens_hi = sd[3], doi_width = D["doi_width"],
        occ = occ, P1 = P1, P2 = P2, nunits = nun,
        have_consensus = have_cons,
        win = [za, zb], win_source = wsrc, win_idx = ix,
        W = W, wgrid = WGRID, PW = [mean(W .>= w) for w in WGRID],
        E = Emain.E, Etop = Emain.top, Ebot = Emain.bot, floor_main = FLOOR_MAIN,
        PE = [mean(Emain.E .>= w) for w in WGRID],
        screen_L = collect(SCREENS),
        screen_best = [screen[L].best for L in SCREENS if haskey(screen, L)],
        screen_top = [screen[L].modal_top for L in SCREENS if haskey(screen, L)],
        screen_bot = [screen[L].modal_bot for L in SCREENS if haskey(screen, L)],
        screen_main_best = S.best, screen_main_tops = S.tops,
        screen_main_bots = S.bots, screen_main_cnt = S.cnt,
        screen_main_top = S.modal_top, screen_main_bot = S.modal_bot,
        pred_V0 = ppq, V0 = D["V0"], q = D["q"], chi2n_draws = chi2n_draws,
        zlim_sigmahat = det_hat.zlim, dmin_sigmahat = det_hat.dmin,
        zlim_detectability = det.zlim, z_respond = z_respond,
        width_frac = wfrac, respond_frac = RESPOND_FRAC,
        zsens_labels = [x.label for x in zsens], zsens_zlim = [x.zlim for x in zsens],
        zsens_top = [x.top for x in zsens], zsens_bot = [x.bot for x in zsens],
        zsens_P = [x.P for x in zsens],
        zsens_W50 = [quantile(x.W, .5) for x in zsens],
        zsens_W05 = [quantile(x.W, .05) for x in zsens],
        zsens_W95 = [quantile(x.W, .95) for x in zsens],
        screen_pinned = pinned,
        screen_nolimit_best = nolimit_best, screen_nolimit_top = nolimit_top,
        screen_nolimit_bot = nolimit_bot,
        theta_min = THETA_MIN, thick_min = THICK_MIN, verdict = verdict,
        ndraws = nd, s_p50 = s50, sigma_eff = σ_eff)
    (field_id = field_id, zlim = zlim, doi_sens = sd[2], P1 = P1,
     win = (za, zb), wsrc = wsrc, W = W, verdict = verdict,
     screen_main = S.best, modal = (S.modal_top, S.modal_bot))
end

function main()
    d = readdlm(joinpath(FieldConfig.cfgpath(FieldConfig.CFG[], "aem_sounds"), "manifest.csv"), ',', String; header = true)
    ci = Dict(strip(h) => i for (i, h) in enumerate(vec(d[2])))
    sites = d[1][:, ci["field_id"]]
    res = []
    for s in sites
        r = try analyse(s) catch e
            @warn "detection failed for $s" exception = (e, catch_backtrace()); nothing
        end
        r !== nothing && push!(res, r)
    end
    println("\n", "="^100)
    @printf("%-14s %7s %7s %7s  %-22s %8s %8s %8s  %s\n", "site", "zlim", "DOIsens",
            "P1unit", "interval (source)", "W p50", "P(W>=2)", "scr p50", "verdict")
    for r in res
        @printf("%-14s %7.1f %7.1f %7.3f  %6.1f-%-6.1f %-8s %8.2f %8.3f %8.2f  %s\n",
                r.field_id, r.zlim, r.doi_sens, r.P1, r.win[1], r.win[2],
                r.wsrc == "consensus interval" ? "(cons)" : "(screen)",
                isempty(r.W) ? NaN : quantile(r.W, .5),
                isempty(r.W) ? NaN : mean(r.W .>= 2),
                quantile(r.screen_main, .5), r.verdict)
    end
    println("="^100)
end

main()
