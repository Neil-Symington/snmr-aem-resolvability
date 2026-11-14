# smoothness.jl -- two tests on a chi2 profile, because one is not enough.
#
# SPIKE test (the original). Flags a node sitting above the straight chord
# through its two neighbours. Catches an isolated stall, as on 3S where chi2
# jumped 28.5 -> 50.5 -> 221.7. BLIND ON A RISING LIMB: on a convex profile
# every interior node lies BELOW its chord, so the test can never fire there
# however badly the refits stall. That is why it reported "none" on 3P across a
# limb running 25.09 -> 29.08 -> 47.88.
#
# CURVATURE test (the replacement). Fits a least-squares quadratic in log10
# sigma to the neighbourhood EXCLUDING the node itself, then compares. A
# quadratic follows a limb, so a node that is too high stands out whether the
# profile is flat, rising or curving.
module Smoothness
using Statistics
export spike_flags, curvature_flags, both_flags

"Node i above the chord through i-1 and i+1 by more than tol."
function spike_flags(x, y; tol = 2.0)
    n = length(y); f = falses(n)
    for i in 2:n-1
        (isfinite(y[i-1]) && isfinite(y[i]) && isfinite(y[i+1])) || continue
        pred = y[i-1] + (y[i+1]-y[i-1]) * (x[i]-x[i-1]) / (x[i+1]-x[i-1])
        y[i] > pred + tol && (f[i] = true)
    end
    f
end

"""
Node i above a quadratic fitted to its neighbours (i-2,i-1,i+1,i+2), by more
than `tol`. `relative` also allows a proportional tolerance, so a node high up a
steep limb is not flagged for being a few per cent off a fit through points that
are themselves hundreds apart.
"""
function curvature_flags(x, y; dex_tol = 0.05)
    n = length(y); f = falses(n); pred = fill(NaN, n)
    # Fit log10(chi2), not chi2. A profile running 25 to 464 is smooth in the
    # log and wildly curved in the linear scale, where a quadratic through
    # neighbours extrapolates to nonsense -- on 3P it predicted 16.5 at a node
    # whose chi2 was 47.9, which is below chiref and cannot be a misfit at all.
    ly = [isfinite(v) && v > 0 ? log10(v) : NaN for v in y]
    for i in 1:n
        isfinite(ly[i]) || continue
        js = [j for j in (i-2, i-1, i+1, i+2) if 1 <= j <= n && isfinite(ly[j])]
        length(js) >= 3 || continue
        A = [ones(length(js)) [x[j] for j in js] [x[j]^2 for j in js]]
        c = A \ [ly[j] for j in js]
        lp = c[1] + c[2]*x[i] + c[3]*x[i]^2
        pred[i] = 10^lp
        # tolerance in dex: a node more than dex_tol above the fitted curve
        ly[i] > lp + dex_tol && (f[i] = true)
    end
    (flags = f, predicted = pred)
end

"Either test firing. Returns both so the reason is visible."
function both_flags(x, y; tol = 2.0, dex_tol = 0.05)
    s = spike_flags(x, y; tol = tol)
    c = curvature_flags(x, y; dex_tol = dex_tol)
    (spike = s, curvature = c.flags, predicted = c.predicted,
     any = s .| c.flags)
end
end
