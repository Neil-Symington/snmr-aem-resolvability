# occam_udf.jl — the constrained Gauss-Newton solver of synthetics/occam_core.jl,
# with the model bounds made arguments instead of module constants.
#
# WHY A COPY. occam_core.jl hard-codes LO, HI = -3.0, 1.0 in log10 conductivity,
# i.e. 1000 ohm.m at the resistive end. That is fine for the Murray-Basin
# synthetics it was written for, and wrong here: the East Kimberley resistivity
# profiles supplied with the SNMR data run to ~2700 ohm.m at surface, so the
# synthetics' bound would clip the model at the resistive end and the admissible
# band would be a statement about the bound rather than about the data.
# Everything else -- the KKT mean-constraint, the lambda ladder, the two-step
# line search, the :minimise / :target modes -- is unchanged, deliberately, so a
# profile computed here means the same thing as one computed there.
#
# Model vector is log10 CONDUCTIVITY, as in synthetics/aem_occam.jl.
module OccamUDF

using Statistics, Printf, LinearAlgebra, SparseArrays
using HiQGA.transD_GP
const SK = transD_GP.SkyTEM1DInversion

export chi2, makeR1, occam, occam_multi, Bounds

const Λ2 = 10 .^ reverse(LinRange(-2.0, 7.0, 8))
const ALPHAS = (1.0, 0.4)
const Β2 = 0.1
const NSTEPS = 25

"Model bounds in log10 conductivity. Default spans 10,000 to 0.3 ohm.m."
struct Bounds
    lo::Float64
    hi::Float64
end
Bounds() = Bounds(-4.0, 0.5)

chi2(aem, m; computeJ = false) =
    (SK.getresidual(aem, m; computeJ = computeJ); norm(aem.W * aem.res)^2)

function makeR1(n)
    I_, J_, V_ = Int[], Int[], Float64[]
    for i in 2:n
        push!(I_, i); push!(J_, i);     push!(V_, 1.0)
        push!(I_, i); push!(J_, i - 1); push!(V_, -1.0)
    end
    sparse(I_, J_, V_, n, n)
end

"""
    occam(aem, mstart, m0, R, ndata; fixidx, fixval, mode, ctype, bounds)

`mode = :minimise` drives chi2 as low as it will go (profile likelihood).
`mode = :target`   stops at the smoothest model reaching chi2 <= ndata.
`ctype = :mean`    fixes only the window's MEAN log10 conductivity, letting the
                   shape vary -- the proper profile for "this unit's conductivity".
`ctype = :uniform` forces every cell in `fixidx` to `fixval` (conservative).
"""
function occam(aem, mstart, m0, R, ndata; fixidx = Int[], fixval = 0.0,
               mode = :minimise, ctype = :uniform, nsteps = NSTEPS,
               bounds::Bounds = Bounds(), verbose = false)
    LO, HI = bounds.lo, bounds.hi
    n = length(m0)
    meanmode = (ctype === :mean) && !isempty(fixidx)
    free = meanmode ? collect(1:n) : setdiff(1:n, fixidx)
    c = zeros(n); meanmode && (c[fixidx] .= 1 / length(fixidx))
    m = copy(mstart)
    if meanmode
        m[fixidx] .+= fixval - mean(m[fixidx])
    elseif !isempty(fixidx)
        m[fixidx] .= fixval
    end
    clamp!(m, LO, HI)
    RtR = Matrix(R'R); target = float(ndata)
    χcur = chi2(aem, m)
    for step in 1:nsteps
        χcur = chi2(aem, m; computeJ = true)
        JtW = aem.J' * aem.W; Wr = aem.W * aem.res
        cands = Tuple{Float64,Float64,Vector{Float64}}[]
        for l² in Λ2
            H = JtW * JtW' .+ l² .* RtR
            g = JtW * Wr .+ l² .* (RtR * m) .+ l² * Β2 .* (m .- m0)
            H[diagind(H)] .+= l² * Β2
            local Δf
            if meanmode
                K = [H c; c' 0.0]
                rhs = vcat(-g, 0.0)
                sol = try; K \ rhs; catch; (K + 1e-8I) \ rhs; end
                Δf = sol[1:n]
            else
                Hff, gf = Symmetric(H[free, free]), g[free]
                Δf = try; -(cholesky(Hff) \ gf); catch; -((Hff + 1e-8I) \ gf); end
            end
            any(isnan, Δf) && continue
            for α in ALPHAS
                mn = copy(m); mn[free] .+= α .* Δf
                if meanmode
                    mn[fixidx] .+= fixval - mean(mn[fixidx])
                elseif !isempty(fixidx)
                    mn[fixidx] .= fixval
                end
                clamp!(mn, LO, HI)
                push!(cands, (chi2(aem, mn), l², mn))
            end
        end
        isempty(cands) && break
        if mode === :target
            ok = filter(c -> c[1] <= target, cands)
            best = isempty(ok) ? cands[argmin(first.(cands))] :
                                 ok[argmax([c[2] for c in ok])]
            improved = best[1] < χcur - 1e-9
            m = best[3]; χcur = best[1]
            verbose && @printf("   step %2d  chi2/N %.3f\n", step, χcur / ndata)
            (χcur <= target || !improved) && break
        else
            best = cands[argmin(first.(cands))]
            improved = best[1] < χcur * (1 - 1e-6) - 1e-9
            !improved && break
            m = best[3]; χcur = best[1]
            verbose && @printf("   step %2d  chi2/N %.3f\n", step, χcur / ndata)
        end
    end
    m, χcur
end

"Best of several starts — guards against a local minimum changing the profile."
function occam_multi(aem, starts, m0, R, ndata; kw...)
    bm = copy(first(starts)); bχ = Inf
    for s in starts
        m, χ = occam(aem, s, m0, R, ndata; kw...)
        if χ < bχ
            bm = m; bχ = χ
        end
    end
    bm, bχ
end

end # module
