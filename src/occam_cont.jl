# occam_cont.jl -- STAGE A of the solver work: continuation in the smoothing
# weight. NEW FILE; it calls OccamCore.occam and changes nothing inside it.
#
# WHY. At lambda2 = 0.01 a solve from a cold or distant start spends most of
# its budget travelling, and runs out: 33 of 40 verification solves stopped at
# ~180 accepted steps on the forward-call budget NFEV_MAX = 1500, not on the
# step cap, with chi2 still falling. A solve at lambda2 = 10 is cheap and
# lands near the smooth solution; using it to start lambda2 = 1, then 0.1,
# then 0.01, means each stage begins close to its own answer.
#
# THE STALL RULE. Each stage also stops once chi2 has moved by less than
# `chistall` over the last `chistall_n` accepted steps. Measured on one traced
# solve: 30 of 118 steps ran after chi2 had stopped moving by 0.01, 88.5% of
# the time being Jacobian evaluation. The rule lives in OccamCore.occam behind
# `chistall`, which DEFAULTS TO ZERO -- with it zero the solve is bit-for-bit
# what it was, and no record run passes it.
#
# This returns the FINAL stage's model and chi2, plus a per-stage trace, so a
# caller can see what each lambda2 bought.
module OccamCont
export occam_continuation, LAM_LADDER

using ..OccamCore

"The default ladder, coarse to fine. The last entry is the answer."
const LAM_LADDER = (10.0, 1.0, 0.1, 0.01)

"""
    occam_continuation(aem, mstart, m0, R, ndata; ladder, chistall, kw...)

Solve at each lambda2 in `ladder`, warm-starting each from the previous.
Returns `(m, chi, nf, trace)`; `trace` has one NamedTuple per stage with
lambda2, chi2, steps, nfev, njac, converged, stalled, conres and nbound.
Any keyword `occam` takes is passed through, so bounds, dz, fixidx, fixval
and ctype are the caller's.
"""
function occam_continuation(aem, mstart, m0, R, ndata;
                            ladder = LAM_LADDER, chistall = 0.01,
                            chistall_n = 3, nsteps = 400, kw...)
    m = copy(mstart)
    chi = Inf; nf = nothing
    trace = NamedTuple[]
    for l2 in ladder
        m, chi, _, nf = OccamCore.occam(aem, m, m0, R, ndata;
                                        lambda2 = l2, nsteps = nsteps,
                                        chistall = chistall,
                                        chistall_n = chistall_n,
                                        info = true, kw...)
        push!(trace, (lambda2 = l2, chi2 = chi, steps = nf.nacc, nfev = nf.nfev,
                      njac = nf.njac, converged = nf.converged,
                      stalled = get(nf, :stalled, false), conres = nf.conres,
                      nbound = nf.nbound, capped = nf.capped))
    end
    (m, chi, nf, trace)
end
end
