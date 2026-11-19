# edgefinder.jl -- warm-started bisection for a threshold crossing, with the
# solver INJECTED so the algorithm can be tested without the AEM code.
#
# `solvefn(x, from) -> (chi = Float64, m = Any)` solves at abscissa x (log10
# sigma) warm-started from model `from`, and is expected to return a LOWER chi2
# the closer `from` is to x's own answer. That continuation behaviour is the
# whole reason bisection is used rather than a fixed grid.
#
# Two bugs this structure exists to prevent, both found in the AEM run of
# 17 Nov 2026:
#
#  1. The bracket was updated from the midpoint's status alone. Re-solving the
#     bracketing nodes can push the "rejected" one BELOW the threshold, which
#     invalidates the bracket silently and under-reports the edge. The bracket
#     is now re-established from the full solved set every refinement and
#     asserted valid before any edge is computed.
#  2. The chi2 convergence test compared the bracket's lower node between
#     successive refinements -- but that is a DIFFERENT node each time, so the
#     quantity measured the profile's slope and could never settle. A node is
#     now compared only with itself.
#  3. `bracket` locates an UPPER crossing only -- highest accepted node, lowest
#     rejected node above it. A caller wanting the LOWER edge must reflect the
#     abscissa, or it silently re-refines the upper edge and returns it.
#     gate_check/edge_bisect.jl did exactly that until 19 Nov 2026: it printed
#     the correct lower bracket, then called this module on the unreflected set.
#     The reflection now lives here, behind `side = :lower`, so no caller has to
#     remember it. stages/stage4_gate_v2.jl:433-441 reflects inline and is
#     correct; it can be collapsed onto this path after the gate freeze.
module EdgeFinder
export bisect_edge

"Highest accepted node, and the lowest rejected node above it."
function bracket(S, thr)
    ks = sort(collect(keys(S)))
    acc = [k for k in ks if S[k].chi <= thr]
    isempty(acc) && return (nothing, nothing, "no accepted node")
    a = maximum(acc)
    above = [k for k in ks if k > a && S[k].chi > thr]
    isempty(above) && return (a, nothing, "bracket lost")
    (a, minimum(above), "valid")
end

nearest(S, x) = argmin(k -> abs(k - x), collect(keys(S)))

"""
    bisect_edge(S, solvefn, thr; ...)

`S` maps abscissa -> (chi, m) for already-solved nodes and is updated in place.
Returns `(edge, converged, state, hist, nsolves)`. CONVERGED requires, on two
successive refinements: a valid bracket; the edge moving under `dex_tol`; and
re-solving each bracket node warm from the newest midpoint changing THAT SAME
node's chi2 by under `chi_tol`.
"""
function bisect_edge(S, solvefn, thr; side = :upper, kwargs...)
    side === :upper && return bisect_upper!(S, solvefn, thr; kwargs...)
    side === :lower || error("side must be :upper or :lower, got $side")
    # reflect: a decreasing crossing becomes an increasing one
    Sr = Dict(-k => v for (k, v) in S)
    rr = bisect_upper!(Sr, (x, from) -> solvefn(-x, from), thr;
                       xlabel = x -> -x, kwargs...)
    for (k, v) in Sr                      # carry new solves back to the caller's set
        (!haskey(S, -k) || v.chi < S[-k].chi) && (S[-k] = v)
    end
    (edge = isfinite(rr.edge) ? -rr.edge : NaN, converged = rr.converged,
     state = rr.state, hist = [merge(h, (edge = -h.edge,)) for h in rr.hist],
     nsolves = rr.nsolves, raised = rr.raised)
end

"The upper crossing. `xlabel` maps an abscissa back to the caller's sign for
messages only; it never enters the arithmetic."
function bisect_upper!(S, solvefn, thr; maxref = 8, dex_tol = 0.02, chi_tol = 0.5,
                     extend_step = 0.1, max_extend = 6, xlabel = identity)
    hist = NamedTuple[]; nsolves = 0; streak = 0
    prev_edge = NaN; edge = NaN; state = "valid"; raised = String[]
    for it in 1:maxref
        a, b, st = bracket(S, thr)
        if st == "no accepted node"
            return (edge = NaN, converged = false, state = st, hist = hist, nsolves = nsolves, raised = raised)
        end
        if st == "bracket lost"
            # extend upward until a rejected node appears, and say so
            x = maximum(keys(S)); found = false
            for _ in 1:max_extend
                x += extend_step
                r = solvefn(x, S[nearest(S, x)].m); nsolves += 1
                S[x] = r
                if r.chi > thr; found = true; break; end
            end
            state = "bracket lost, extended"
            found || return (edge = NaN, converged = false,
                             state = "bracket lost, extension exhausted",
                             hist = hist, nsolves = nsolves, raised = raised)
            a, b, st = bracket(S, thr)
        end
        (S[a].chi <= thr && S[b].chi > thr) ||
            error("invalid bracket after re-establishment: chi2($a)=$(S[a].chi), chi2($b)=$(S[b].chi), thr=$thr")

        mid = (a + b) / 2
        rm = solvefn(mid, S[nearest(S, mid)].m); nsolves += 1
        if !haskey(S, mid) || rm.chi < S[mid].chi
            S[mid] = rm            # store the solver's own record, not a copy of two fields
        end
        # re-solve each bracket node warm from the new midpoint; compare each
        # node WITH ITSELF, which is criterion (c)
        dmax = 0.0
        for q in (a, b)
            rq = solvefn(q, S[mid].m); nsolves += 1
            d = abs(rq.chi - S[q].chi); dmax = max(dmax, d)
            if rq.chi < S[q].chi - 1e-12
                S[q] = rq
            elseif rq.chi > S[q].chi + chi_tol
                push!(raised, "node $(round(xlabel(q),digits=4)): $(round(S[q].chi,digits=3)) -> $(round(rq.chi,digits=3))")
            end
        end
        # re-establish and interpolate
        a2, b2, st2 = bracket(S, thr)
        if st2 != "valid"
            push!(hist, (step = abs(b - a), edge = NaN, chi_lo = NaN, chi_hi = NaN,
                         move = NaN, dchi = dmax, valid = false, nsolves = nsolves))
            streak = 0; state = st2; continue
        end
        ca, cb = S[a2].chi, S[b2].chi
        edge = cb > ca ? a2 + (thr - ca) * (b2 - a2) / (cb - ca) : (a2 + b2)/2
        move = isnan(prev_edge) ? NaN : abs(edge - prev_edge)
        push!(hist, (step = abs(b2 - a2), edge = edge, chi_lo = ca, chi_hi = cb,
                     move = move, dchi = dmax, valid = true, nsolves = nsolves))
        streak = (!isnan(move) && move < dex_tol && dmax < chi_tol) ? streak + 1 : 0
        prev_edge = edge
        streak >= 2 && return (edge = edge, converged = true, state = "located",
                               hist = hist, nsolves = nsolves, raised = raised)
    end
    (edge = edge, converged = false, state = state == "valid" ? "not converged" : state,
     hist = hist, nsolves = nsolves, raised = raised)
end
end
