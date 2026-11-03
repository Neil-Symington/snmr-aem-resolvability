# occam_core.jl — the constrained Gauss-Newton ("constrained Occam") solver,
# shared by occam_constrained.jl (the profile) and occam_coverage.jl (the
# posterior-coverage audit). Lifted verbatim out of occam_constrained.jl so the
# two steps cannot drift apart; the constants below are its original ones and
# must not be changed casually — stored profiles were computed with them.
#
# Gauss-Newton in log10-conductivity with R1 smoothness, as HiQGA's gradientinv,
# except that (a) cells inside a depth window are held FIXED (or their MEAN is
# fixed) and only the free cells are solved for, and (b) :minimise drives chi2 as
# low as it will go rather than stopping at a target.
#
# The constraint is in the normal equations, not a clamp: with Delta_fixed = 0 the
# system restricts to H[free,free], g[free]. Smoothness across the contact is
# still penalised (R'R*m uses the FULL model), so the fit cannot buy misfit with
# an artificial jump at the boundary.
module OccamCore

using Statistics, Printf, LinearAlgebra, SparseArrays
using HiQGA.transD_GP
const SK = transD_GP.SkyTEM1DInversion

export chi2, makeR1, occam, occam_multi, setbounds!, getbounds

# The solver's clamp on log10 sigma. Held in Refs rather than as consts so the
# open-edge check can lower the floor without editing the file: the question
# "is the band's upper edge set by the data or by sigma_min?" can only be asked
# by moving sigma_min. setbounds!(-4.0, 1.0) gives sigma_min = 1e-4 S/m.
const LO_REF = Ref(-3.0)
const HI_REF = Ref(1.0)
setbounds!(lo, hi) = (LO_REF[] = lo; HI_REF[] = hi; nothing)
getbounds() = (LO_REF[], HI_REF[])
const Λ2 = 10 .^ reverse(LinRange(-2.0, 7.0, 8))
const ALPHAS = (1.0, 0.4)
const Β2 = 0.1
const NSTEPS = 25
# Levenberg damping for :map, adapted per step rather than swept. Trying every
# value on a ladder at every iteration costs one forward model per entry per
# alpha and was ~12x slower than it needed to be; the standard LM update
# (shrink mu on an accepted step, grow it on a rejected one) converges in the
# same number of steps at a fraction of the evaluations.
const MU0, MU_DOWN, MU_UP, MU_MIN, MU_MAX = 1e-2, 1/3, 3.0, 1e-10, 1e8
# RELATIVE convergence tolerance. An absolute test (phi_new < phi - 1e-9) never
# fires on phi ~ 1e4: a 1e-10 improvement still counts as progress, so the fit
# runs to the step cap long after it has converged. Measured: phi identical to
# three decimals from step 13 to step 28, at ~2 s per step.
const PHI_RTOL = 1e-9
const NSTEPS_MAP = 120
# backtracking line search along the Gauss-Newton direction: the cheap fix for
# a good direction with a bad length, which is most rejected LM steps.
const NLS_TRIES = 6
const LS_SHRINK = 0.35
# Convergence and cost control for :map. PHI_RTOL (1e-9) is the ACCEPTANCE
# threshold -- it decides whether a trial step is downhill at all. It is far too
# tight to use as a stopping rule: a step improving Phi ~ 1e3 by 1e-6 passes it,
# so without a separate test the loop always runs to the step cap making
# meaningless progress. PHI_CONV is the stopping rule, and NFEV_MAX bounds the
# work when neither test fires.
const PHI_CONV = 1e-7
const NFEV_MAX = 1500
const MU_TRIES = 8

chi2(aem, m; computeJ=false) = (SK.getresidual(aem, m; computeJ=computeJ); norm(aem.W*aem.res)^2)

function makeR1(n)
    I_, J_, V_ = Int[], Int[], Float64[]
    for i in 2:n
        push!(I_,i); push!(J_,i);   push!(V_, 1.0)
        push!(I_,i); push!(J_,i-1); push!(V_,-1.0)
    end
    sparse(I_, J_, V_, n, n)
end

"""
    occam(aem, mstart, m0, R, ndata; fixidx, fixval, mode, ctype)

`mode = :minimise` drives chi2 as low as it will go (profile likelihood).
`mode = :target`   stops at the smoothest model reaching chi2 <= ndata (Occam).
`mode = :map`      MAP under a FIXED Gaussian prior, lambda2 supplied by the
                   caller and identical at every grid point.

WHY :map EXISTS. :minimise and :target both pick lambda2 off the ladder at every
iteration, keeping whichever candidate lowers chi2 most. That makes lambda2 a
damping device, not a prior: it differs between grid points and drifts small at
convergence, where J'WJ is rank-deficient (fewer data than cells). An integrated
likelihood

    log L = -chi2/2 - 1/2 logdet(J'WJ + P),     P = lambda2*(R'R + beta2*I)

is a Laplace approximation to the integral over the free cells, and it is only
coherent if P is the SAME prior precision everywhere -- otherwise the dropped
constants 1/2 log|P| and (p/2) log 2pi do not cancel on normalisation, and the
determinant term inherits the optimiser's path rather than the data's
information. :map fixes lambda2 so that formula means what it says. lambda2 then
becomes a declared prior to be varied in the robustness sweep, not a hidden
tuning knob.
`ctype = :uniform` forces every cell in `fixidx` to `fixval` (strict: removes
shape freedom too, so its Dchi2 OVERSTATES the cost and the interval it gives is
conservative). `ctype = :mean` fixes only the window's mean log10-conductivity
and lets the shape vary — the proper profile for "this unit's conductivity".
"""
function occam(aem, mstart, m0, R, ndata; fixidx=Int[], fixval=0.0, mode=:minimise,
               ctype=:uniform, nsteps=NSTEPS, verbose=false,
               lambda2=nothing, beta2=Β2, withH=false, dz=nothing, info=false)
    mode === :map && (nsteps = max(nsteps, NSTEPS_MAP))
    n = length(m0)
    meanmode = (ctype === :mean) && !isempty(fixidx)
    condmode = (ctype === :cond) && !isempty(fixidx)
    if condmode
        mode === :map || throw(ArgumentError("ctype=:cond is only implemented for mode=:map"))
        dz === nothing && throw(ArgumentError("ctype=:cond needs dz, the cell thicknesses"))
        fixval > 0 || throw(ArgumentError("ctype=:cond takes fixval = conductance in S, > 0"))
    end
    free = (meanmode || condmode) ? collect(1:n) : setdiff(1:n, fixidx)
    # c is the constraint row: c'm = fixval. For :mean it averages the interval.
    c = zeros(n); meanmode && (c[fixidx] .= 1/length(fixidx))
    # ---- :cond forces CONDUCTANCE, sum 10^m_i dz_i over fixidx = fixval.
    # That is nonlinear in m, so unlike :mean it needs (a) the constraint row
    # relinearised at every iterate, a_i = ln10 * 10^m_i * dz_i, (b) the step
    # projected onto its tangent null space, and (c) the model re-projected back
    # onto the constraint surface after the step. (c) is exact and closed-form:
    # scaling every cell in the unit by a common factor multiplies the
    # conductance by that factor, so m[fixidx] .+= log10(C/g(m)) lands on the
    # surface in one shot. clamp! can knock it off again, hence the short loop.
    dzf = condmode ? Float64[dz[i] for i in fixidx] : Float64[]
    gcond(mm) = sum(10 .^ mm[fixidx] .* dzf)
    function projcond!(mm)
        for _ in 1:12
            gc = gcond(mm)
            (isfinite(gc) && gc > 0) || break
            δ = log10(fixval / gc)
            abs(δ) < 1e-13 && break
            mm[fixidx] .+= δ
            clamp!(mm, LO_REF[], HI_REF[])
        end
        mm
    end
    # constraint row and the linearised residual a'Δ = rc, at the current m
    function conlin(mm)
        if condmode
            a = zeros(n)
            a[fixidx] .= log(10) .* (10 .^ mm[fixidx]) .* dzf
            return a, fixval - gcond(mm)
        else
            return c, fixval - dot(c, mm)
        end
    end
    # put the start on the constraint surface
    m = copy(mstart)
    if condmode
        clamp!(m, LO_REF[], HI_REF[]); projcond!(m)
    elseif meanmode
        m[fixidx] .+= fixval - mean(m[fixidx])        # put the start on the constraint
    elseif !isempty(fixidx)
        m[fixidx] .= fixval
    end
    clamp!(m, LO_REF[], HI_REF[])
    condmode && projcond!(m)
    RtR = Matrix(R'R); target = float(ndata)
    if mode === :map
        lambda2 === nothing &&
            throw(ArgumentError("mode=:map needs an explicit lambda2 -- it is the prior, not a knob"))
    end
    # In :minimise/:target the ladder IS lambda2: it is simultaneously the
    # smoothness prior and the Levenberg damping, re-chosen every step. :map
    # cannot do that -- the prior has to stay put -- so the ladder becomes a
    # damping term mu added to the diagonal of the STEP Hessian only. The prior
    # stays at lambda2, and the Hessian reported by withH excludes mu, which is
    # zero at a converged point anyway.
    ladder = Λ2
    Β2eff = float(beta2)
    l2fix = mode === :map ? float(lambda2) : NaN
    # MAP minimises the PENALISED objective, not the misfit. Using argmin(chi2)
    # here (as :minimise does) walks toward the maximum-likelihood model, and the
    # Hessian would then be evaluated somewhere that is not the MAP.
    pen(mm) = l2fix * (sum(abs2, R*mm) + Β2eff * sum(abs2, mm .- m0))
    obj(χ, mm) = mode === :map ? χ + pen(mm) : χ
    χcur = chi2(aem, m)

    if mode === :map
        # Levenberg-Marquardt on the penalised objective, fixed prior.
        #
        # THE ABORT BUG, fixed here. The previous version ran the mu ladder
        # MU_TRIES = 8 times inside one step and then did `accepted || break`,
        # abandoning the whole optimisation. Starting from MU0 = 1e-2 with
        # MU_UP = 3, that reaches mu = 21.9 at best, so any first step needing
        # heavier damping than that returned mstart UNCHANGED -- measured at
        # tag 2 and tag 2S, where chi2 came back equal to chi2(m0) to 12
        # figures and max|m - m0| was exactly zero, for every lambda2 below the
        # one node that happened to escape. ML-II then "selected" that node.
        # Two changes: the mu ladder now runs to MU_MAX within a single step
        # (one Jacobian per step, not one per ladder entry), so `accepted ||
        # break` means damping was exhausted rather than merely rationed; and a
        # backtracking line search along the Gauss-Newton direction handles the
        # common case where the direction is fine and only its length is wrong,
        # which is far cheaper than re-solving with a larger mu.
        μ = MU0
        Φcur = χcur + pen(m)
        nacc = 0; gnorm = NaN; exhausted = false; nj = 0; nfev = 0
        converged = false; budget = false
        for step in 1:nsteps
            χcur = chi2(aem, m; computeJ = true); nj += 1; nfev += 1
            Φcur = χcur + pen(m)
            JtW = aem.J' * aem.W; Wr = aem.W * aem.res
            Hb = JtW*JtW' .+ l2fix.*RtR
            g  = JtW*Wr .+ l2fix.*(RtR*m) .+ l2fix*Β2eff .* (m .- m0)
            Hb[diagind(Hb)] .+= l2fix*Β2eff
            a, rc = conlin(m)
            # projected gradient: the stationarity measure that actually applies
            # on a constrained surface (the raw gradient need not vanish there)
            if meanmode || condmode
                na = norm(a); gnorm = na > 0 ? norm(g .- (dot(a, g)/na^2).*a) : norm(g)
            else
                gnorm = norm(g[free])
            end
            accepted = false
            while true
                H = copy(Hb); H[diagind(H)] .+= μ
                local Δf
                if meanmode || condmode
                    K = [H a; a' 0.0]
                    sol = try; K \ vcat(-g, rc); catch; (K + 1e-8I) \ vcat(-g, rc); end
                    Δf = sol[1:n]
                else
                    Hff, gf = Symmetric(H[free,free]), g[free]
                    Δf = try; -(cholesky(Hff) \ gf); catch; -((Hff + 1e-8I) \ gf); end
                end
                if !any(isnan, Δf)
                    t = 1.0
                    for _ in 1:NLS_TRIES          # backtracking line search
                        mn = copy(m); mn[free] .+= t .* Δf
                        if condmode
                            clamp!(mn, LO_REF[], HI_REF[]); projcond!(mn)
                        elseif meanmode
                            mn[fixidx] .+= fixval - mean(mn[fixidx])
                            clamp!(mn, LO_REF[], HI_REF[])
                        else
                            isempty(fixidx) || (mn[fixidx] .= fixval)
                            clamp!(mn, LO_REF[], HI_REF[])
                        end
                        χn = chi2(aem, mn); Φn = χn + pen(mn); nfev += 1
                        if Φn < Φcur - max(PHI_RTOL*abs(Φcur), 1e-12)
                            rel = (Φcur - Φn) / max(abs(Φcur), 1e-12)
                            m = mn; χcur = χn; Φcur = Φn
                            μ = max(μ*MU_DOWN, MU_MIN); accepted = true; nacc += 1
                            rel < PHI_CONV && (converged = true)
                            verbose && @printf("   step %2d  chi2 %.3f  phi %.3f  mu %.2g  t %.3g  rel %.2g\n",
                                               step, χn/ndata, Φn, μ, t, rel)
                            break
                        end
                        t *= LS_SHRINK
                    end
                end
                accepted && break
                if nfev >= NFEV_MAX; budget = true; break; end
                if μ >= MU_MAX*(1 - 1e-12); exhausted = true; break; end
                μ = min(μ*MU_UP, MU_MAX)
            end
            # Stop on genuine convergence, on the work budget, or when damping
            # was exhausted without an accepted step. Only the last of those is
            # the old `accepted || break`, and it now means what it says.
            (accepted && !converged && !budget) || break
        end
        if withH || info
            chi2(aem, m; computeJ = true)
            JtW = aem.J' * aem.W
            H = JtW*JtW' .+ l2fix.*RtR
            H[diagind(H)] .+= l2fix*Β2eff
            if info
                nlo = count(x -> x <= LO_REF[] + 1e-9, m)
                nhi = count(x -> x >= HI_REF[] - 1e-9, m)
                cres = condmode ? abs(gcond(m) - fixval)/fixval :
                       meanmode ? abs(dot(c, m) - fixval) : 0.0
                nfo = (pen = pen(m), phi = χcur + pen(m), gnorm = gnorm,
                       nbound = nlo + nhi, nlo = nlo, nhi = nhi,
                       conres = cres, nacc = nacc, njac = nj, nfev = nfev,
                       mu = μ, exhausted = exhausted, converged = converged,
                       budget = budget, nsteps = nsteps, capped = (nj >= nsteps),
                       # A legitimate stop is either the relative-Phi test
                       # firing, or no descent direction being found at any mu
                       # up to MU_MAX after real progress. Running out of steps
                       # while Phi is still falling, or hitting the work budget,
                       # is not convergence -- and the chi2 such a fit reports
                       # is an upper bound, which is exactly the thing that must
                       # not be allowed to win an ML-II comparison.
                       settled = (nacc >= 1) && (converged || exhausted) &&
                                 !budget && (nj < nsteps),
                       stuck = (nacc == 0), interior = (nlo + nhi) == 0)
                return m, χcur, Symmetric(H), nfo
            end
            return m, χcur, Symmetric(H)
        end
        return m, χcur
    end

    for step in 1:nsteps
        χcur = chi2(aem, m; computeJ=true)
        JtW = aem.J' * aem.W; Wr = aem.W * aem.res
        cands = Tuple{Float64,Float64,Vector{Float64},Float64}[]   # (objective, l2/mu, model, chi2)
        for l² in ladder
            lp = mode === :map ? l2fix : l²        # prior weight
            mu = mode === :map ? l² : 0.0          # Levenberg damping, step only
            H = JtW*JtW' .+ lp.*RtR
            g = JtW*Wr .+ lp.*(RtR*m) .+ lp*Β2eff .* (m .- m0)
            H[diagind(H)] .+= lp*Β2eff + mu
            local Δf
            if meanmode
                # Equality-constrained Newton step via the KKT system
                #   [H  c] [Δ]   [-g]
                #   [c' 0] [λ] = [ 0]
                # Projecting the UNCONSTRAINED Newton step onto null(c') instead is
                # NOT the same thing and gives a badly scaled direction whenever H
                # couples the interval to its neighbours — which it does.
                K = [H c; c' 0.0]
                rhs = vcat(-g, 0.0)
                sol = try; K \ rhs; catch; (K + 1e-8I) \ rhs; end
                Δf = sol[1:n]
            else
                Hff, gf = Symmetric(H[free,free]), g[free]
                Δf = try; -(cholesky(Hff) \ gf); catch; -((Hff + 1e-8I) \ gf); end
            end
            any(isnan, Δf) && continue
            for α in ALPHAS
                mn = copy(m); mn[free] .+= α .* Δf
                if meanmode
                    mn[fixidx] .+= fixval - mean(mn[fixidx])   # clamp() can drift it
                elseif !isempty(fixidx)
                    mn[fixidx] .= fixval
                end
                clamp!(mn, LO_REF[], HI_REF[])
                χn = chi2(aem, mn); push!(cands, (obj(χn, mn), l², mn, χn))
            end
        end
        isempty(cands) && break
        if mode === :target
            ok = filter(c -> c[1] <= target, cands)
            best = isempty(ok) ? cands[argmin(first.(cands))] : ok[argmax([c[2] for c in ok])]
            improved = best[1] < χcur - 1e-9
            m = best[3]; χcur = best[4]
            verbose && @printf("   step %2d  chi2 %.3f\n", step, χcur/ndata)
            (χcur <= target || !improved) && break
        else
            # :minimise and :map share this branch. With a one-entry ladder the
            # argmin is over the line-search alphas only, which is what :map wants.
            best = cands[argmin(first.(cands))]
            Φcur = mode === :map ? χcur + pen(m) : χcur
            improved = best[1] < Φcur*(1 - 1e-9) - 1e-9
            !improved && break
            m = best[3]; χcur = best[4]
            verbose && @printf("   step %2d  chi2 %.3f\n", step, χcur/ndata)
        end
    end
    if withH
        # Laplace Hessian at the MAP, in the SAME parameterisation the prior is
        # written in. Recomputed from the converged model rather than reused from
        # the last trial step, which belongs to a line-search candidate.
        chi2(aem, m; computeJ=true)
        l² = float(lambda2 === nothing ? last(Λ2) : lambda2)
        JtW = aem.J' * aem.W
        H = JtW*JtW' .+ l².*RtR
        H[diagind(H)] .+= l²*Β2eff
        return m, χcur, Symmetric(H)
    end
    m, χcur
end

"Best of several starts — guards against a local minimum changing the profile."
function occam_multi(aem, starts, m0, R, ndata; kw...)
    bm = copy(first(starts)); bχ = Inf
    for s in starts
        m, χ = occam(aem, s, m0, R, ndata; kw...)
        if χ < bχ
            bm = m
            bχ = χ
        end
    end
    bm, bχ
end

end # module
