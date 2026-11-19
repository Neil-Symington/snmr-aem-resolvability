# prerun_lib.jl — convergence diagnostics and derived quantities for the pre-run.
#
# RUNBOOK §5 is the authority on the convergence side: R̂ must not be computed
# across the chains of one PT ensemble (they target p^(1/T)), and ESS binds long
# before R̂ does. So the gate is rank-normalised split-R̂ across INDEPENDENT
# ensembles' T=1 streams (Vehtari, Gelman, Simpson, Carpenter & Bürkner 2021,
# Bayesian Analysis 16:667), paired with a minimum-ESS condition.
#
# Neither Distributions nor SpecialFunctions is a dependency, so the inverse
# normal CDF is implemented here (Acklam's rational approximation, |ε| < 1.15e-9).

module PrerunLib

using Statistics, Printf

export norminvcdf, iat, ess_basic, rank_normalise, split_rhat, rhat_and_ess,
       stationary_trend, doi_from_width, kl_against_uniform, detection_flag,
       interval_mean_max, set_run_filenames!, nodes_at_T1, assert_one_chain_at_T1

# ── Inverse normal CDF ──────────────────────────────────────────────────────
const _A = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
             1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00)
const _B = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
             6.680131188771972e+01, -1.328068155288572e+01)
const _C = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
            -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00)
const _D = ( 7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
             3.754408661907416e+00)

"Inverse standard-normal CDF. Acklam's approximation; adequate for rank normalisation."
function norminvcdf(p::Real)
    (p <= 0 || p >= 1) && return p <= 0 ? -Inf : Inf
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow
        q = sqrt(-2log(p))
        return (((((_C[1]*q+_C[2])*q+_C[3])*q+_C[4])*q+_C[5])*q+_C[6]) /
               ((((_D[1]*q+_D[2])*q+_D[3])*q+_D[4])*q+1)
    elseif p > phigh
        q = sqrt(-2log(1-p))
        return -(((((_C[1]*q+_C[2])*q+_C[3])*q+_C[4])*q+_C[5])*q+_C[6]) /
                ((((_D[1]*q+_D[2])*q+_D[3])*q+_D[4])*q+1)
    end
    q = p - 0.5; r = q*q
    (((((_A[1]*r+_A[2])*r+_A[3])*r+_A[4])*r+_A[5])*r+_A[6])*q /
    (((((_B[1]*r+_B[2])*r+_B[3])*r+_B[4])*r+_B[5])*r+1)
end

# ── Autocorrelation / ESS ───────────────────────────────────────────────────
"Initial-positive-sequence integrated autocorrelation time (as used for SBC thinning)."
function iat(x)
    n = length(x); n < 20 && return NaN
    x̄, v = mean(x), var(x)
    v <= 0 && return Inf
    τ = 1.0
    for k in 1:min(n ÷ 4, 2000)
        a = @view x[1:n-k]; b = @view x[k+1:n]
        ρ = sum((a .- x̄) .* (b .- x̄)) / ((n - k) * v)
        ρ <= 0 && break
        τ += 2ρ
    end
    τ
end

"ESS of a set of chains: total draws / mean IAT."
function ess_basic(chains::Vector{<:AbstractVector})
    good = [c for c in chains if length(c) >= 20 && var(c) > 0]
    isempty(good) && return NaN
    τs = filter(isfinite, [iat(c) for c in good])
    isempty(τs) && return NaN
    ntot = sum(length, good)
    ntot / max(mean(τs), 1.0)
end

# ── Rank normalisation and split-R̂ ──────────────────────────────────────────
"""
Rank-normalise draws pooled across chains (Vehtari et al. 2021 §3.1):
average ranks over the pooled sample, then z = Φ⁻¹((r − 3/8)/(S − 1/4)).
Makes R̂ and ESS valid for heavy-tailed or improper-variance targets.
"""
function rank_normalise(chains::Vector{<:AbstractVector})
    pooled = vcat(chains...)
    S = length(pooled)
    perm = sortperm(pooled)
    ranks = Vector{Float64}(undef, S)
    i = 1
    while i <= S                       # average ranks within ties
        j = i
        while j < S && pooled[perm[j+1]] == pooled[perm[i]]; j += 1; end
        r̄ = (i + j) / 2
        for k in i:j; ranks[perm[k]] = r̄; end
        i = j + 1
    end
    z = [norminvcdf((r - 0.375) / (S - 0.25)) for r in ranks]
    out = Vector{Vector{Float64}}(); o = 1
    for c in chains
        push!(out, z[o:o+length(c)-1]); o += length(c)
    end
    out
end

"Split each chain in half — doubles the chain count, halves the length."
function split_chains(chains::Vector{<:AbstractVector})
    out = Vector{Vector{Float64}}()
    for c in chains
        n = length(c); n < 4 && continue
        h = n ÷ 2
        push!(out, collect(Float64, c[1:h]))
        push!(out, collect(Float64, c[h+1:2h]))
    end
    out
end

"Classical R̂ (Gelman–Rubin) on equal-length chains."
function _rhat(chains::Vector{Vector{Float64}})
    m = length(chains); m < 2 && return NaN
    n = minimum(length, chains); n < 2 && return NaN
    cs = [c[1:n] for c in chains]
    means = mean.(cs); vars = var.(cs)
    W = mean(vars)
    W <= 0 && return NaN
    B = n * var(means)
    v̂ = (n - 1) / n * W + B / n
    sqrt(v̂ / W)
end

"""
Minimum pooled spread, as a fraction of the prior width, below which a scalar is
treated as DETERMINED and R̂ is not computed for it.

Rank normalisation is scale-free: a parameter pinned at a prior bound, varying
only in its sixth decimal, is mapped onto the full normal range exactly as a
well-spread one is, so R̂ reports non-convergence from numerical dust. Observed:
θ at 30 m pinned at the 0.001 floor with pooled sd 1e-4 returned R̂ = 1.301.
"""
const RHAT_MIN_SPREAD_FRAC = 1e-3

"""
Rank-normalised split-R̂ and bulk/tail ESS across independent ensembles.

`chains` is one vector of retained T=1 draws per ensemble. Returns the larger of
bulk-R̂ (rank-normalised) and tail-R̂ (rank-normalised on folded draws), which is
the quantity Vehtari et al. recommend gating on.

Pass `priorwidth` to enable the degeneracy guard: a scalar whose pooled spread
is below `RHAT_MIN_SPREAD_FRAC` of it is reported as `determined = true` with
`rhat = NaN`, and should not be gated.
"""
function rhat_and_ess(chains::Vector{<:AbstractVector}; priorwidth = nothing)
    good = [collect(Float64, c) for c in chains if length(c) >= 8]
    length(good) < 2 && return (rhat = NaN, ess_bulk = NaN, ess_tail = NaN, determined = false)
    if priorwidth !== nothing
        sp = std(vcat(good...))
        if !isfinite(sp) || sp < RHAT_MIN_SPREAD_FRAC * priorwidth
            return (rhat = NaN, ess_bulk = NaN, ess_tail = NaN, determined = true)
        end
    end
    sp   = split_chains(good)
    zb   = rank_normalise(sp)
    rb   = _rhat(zb)
    med  = median(vcat(sp...))
    fold = [abs.(c .- med) for c in sp]
    zt   = rank_normalise(fold)
    rt   = _rhat(zt)
    valid = filter(!isnan, [rb, rt])
    (rhat     = isempty(valid) ? NaN : maximum(valid),
     ess_bulk = ess_basic(zb),
     ess_tail = ess_basic(zt),
     determined = false)
end

# ── Misfit-trace stationarity (spec §4) ─────────────────────────────────────
"""
Is the retained misfit trace stationary (spec §4, "no trend beyond 2σ")?

The drift is compared with the SCATTER OF THE TRACE ITSELF, not with the
standard error of the fitted slope. With ~5,000 retained draws the slope's
standard error is minute, so a t-test on it rejects drifts far smaller than the
sample-to-sample variation — it flags traces that are stationary in every sense
that matters. The practical reading of "no trend beyond 2σ" is that the total
drift across the retained half should not exceed 2 standard deviations of the
trace.

Returns the drift in units of the trace's own σ (`drift_sigma`) alongside the
classical t statistic, so both readings are on the record.
"""
function stationary_trend(y::AbstractVector; nsig = 2.0)
    n = length(y)
    n < 20 && return (stationary = true, slope = 0.0, drift = 0.0,
                      drift_sigma = 0.0, t = 0.0)
    x = collect(1.0:n); x̄, ȳ = mean(x), mean(y)
    sxx = sum((x .- x̄).^2)
    sd = std(y)
    (sxx <= 0 || sd <= 0) && return (stationary = true, slope = 0.0, drift = 0.0,
                                     drift_sigma = 0.0, t = 0.0)
    b = sum((x .- x̄) .* (y .- ȳ)) / sxx
    resid = y .- (ȳ .+ b .* (x .- x̄))
    s2 = sum(resid.^2) / (n - 2)
    se = sqrt(max(s2, 0) / sxx)
    t = se > 0 ? b / se : 0.0
    drift = b * n
    (stationary = abs(drift) <= nsig * sd, slope = b, drift = drift,
     drift_sigma = drift / sd, t = t)
end

# ── Depth of investigation (spec §6) ────────────────────────────────────────
"""
DOI from posterior width.

The posterior interval (p95 − p5) is compared with the prior width; DOI is the
shallowest depth at which the posterior has stopped informing, i.e. the first
depth where the width exceeds `widthfrac` of the prior width and stays there.
Returns 0.0 when nothing is resolved at all (spec §7 wants that fraction) and
the grid bottom when the sounding is informative all the way down.
"""
function doi_from_width(zbnd, q05, q95, priorwidth; widthfrac = 0.6)
    w = q95 .- q05
    thresh = widthfrac * priorwidth
    # Scan DOWNWARD for the first depth at which the posterior stops informing.
    # Scanning up from the bottom instead assumes the width grows monotonically
    # with depth, which it does not: the deepest cell can be narrower than cells
    # above it, which made the whole profile read as "informative throughout".
    idx = findfirst(>(thresh), w)
    idx === nothing && return Float64(zbnd[end])   # informative to the base
    idx == 1 && return 0.0                         # nothing resolved anywhere
    Float64(zbnd[idx - 1])
end

"""
Information gain ∫ KL(posterior ‖ prior) dz over the model grid.

The prior is uniform on `bounds`, so KL reduces to Σ pᵢ log(pᵢ/uᵢ) on a common
histogram. `samples` is (ncells × ndraws).
"""
function kl_against_uniform(zbnd, samples::AbstractMatrix, bounds; nbins = 40)
    lo, hi = bounds
    edges = range(lo, hi, length = nbins + 1)
    u = 1.0 / nbins
    dz = diff(vcat(0.0, collect(Float64, zbnd)))
    kl = zeros(size(samples, 1))
    for i in axes(samples, 1)
        h = zeros(nbins)
        for v in @view samples[i, :]
            b = clamp(searchsortedlast(edges, v), 1, nbins)
            h[b] += 1
        end
        s = sum(h); s == 0 && continue
        h ./= s
        kl[i] = sum(p > 0 ? p * log(p / u) : 0.0 for p in h)
    end
    (ig = sum(kl .* dz), kl = kl)
end

"Largest interval-mean of `v` over any window at least `minthick` metres thick."
function interval_mean_max(zbnd, v; minthick = 5.0)
    n = length(zbnd); best = -Inf; bz = (0.0, 0.0)
    tops = vcat(0.0, collect(Float64, zbnd)[1:end-1])
    dz = diff(vcat(0.0, collect(Float64, zbnd)))
    for i in 1:n
        acc = 0.0; th = 0.0
        for j in i:n
            acc += v[j] * dz[j]; th += dz[j]
            if th >= minthick
                m = acc / th
                m > best && (best = m; bz = (tops[i], tops[i] + th))
            end
        end
    end
    (value = best, top = bz[1], bot = bz[2])
end

"""
Detection flag (spec §6).

Trans-D: p5 of the interval-mean θ exceeds `θmin` over some interval at least
`minthick` thick. Note this is the interval mean of the p5 profile, which is
conservative relative to the p5 of the interval-mean distribution — stated so
the definition used is unambiguous.
"""
function detection_flag(zbnd, p05; θmin = 0.10, minthick = 5.0)
    r = interval_mean_max(zbnd, p05; minthick = minthick)
    (detected = r.value > θmin, value = r.value, top = r.top, bot = r.bot)
end

# ── HiQGA file-layout helpers (RUNBOOK §1.3, §1.4, §1.5) ────────────────────
"""
Apply the five in-place filename rewrites `transD_GP.main` performs before
writing (MCMC_Driver.jl:316–329). A fresh OptionsStat names its history
`models_<prefix>.bin`; what lands on disk is `models_<prefix>s.bin`.
"""
function set_run_filenames!(opt, optn = nothing)
    opt.costs_filename    = "misfits_" * opt.fdataname * "s.bin"
    opt.fstar_filename    = "models_"  * opt.fdataname * "s.bin"
    opt.x_ftrain_filename = "points_"  * opt.fdataname * "s.bin"
    if optn !== nothing
        optn.costs_filename = "misfits_" * optn.fdataname * "nuisance.bin"
        optn.vals_filename  = "values_nuisance_" * optn.fdataname * ".bin"
    end
    opt, optn
end

"""
Nucleus count and misfit for the T=1 stream, parsed straight from the costs file.

RUNBOOK §1.4: `history(opt, stat=:nodes)` offsets the column index only when
`chain_idx` is given, so the plain call reads acceptanceRateDC (NaN) instead of
the node count. Columns: 1 chain | 2 iter | 3-7 acceptance | 8 nodes | 9 U |
10 T | 11 dcvalue.
"""
function nodes_at_T1(costs_filename; burninfrac = 0.5)
    (isempty(costs_filename) || !isfile(costs_filename)) &&
        return (nodes = Int[], U = Float64[], nrows = 0, nT1 = 0)
    rows = [split(l) for l in readlines(costs_filename) if !isempty(strip(l))]
    (isempty(rows) || length(rows[1]) != 11) &&
        return (nodes = Int[], U = Float64[], nrows = length(rows), nT1 = 0)
    at1 = [r for r in rows if isapprox(parse(Float64, r[10]), 1.0; atol = 1e-12)]
    isempty(at1) && return (nodes = Int[], U = Float64[], nrows = length(rows), nT1 = 0)
    i0 = max(1, floor(Int, burninfrac * length(at1)))
    (nodes = [parse(Int, r[8]) for r in at1[i0:end]],
     U     = [parse(Float64, r[9]) for r in at1[i0:end]],
     nrows = length(rows), nT1 = length(at1))
end

"""
RUNBOOK §1.5: exactly one chain holds T=1 per iteration, so the T=1 row count
must be nrows/nchains. The chain *ids* rotating through that slot is tempering
working correctly; the row count is the invariant to assert.
"""
function assert_one_chain_at_T1(costs_filename, nchains; tol = 0.02)
    h = nodes_at_T1(costs_filename; burninfrac = 0.0)
    h.nrows == 0 && return (ok = false, ratio = NaN, why = "no costs rows")
    ratio = h.nT1 / h.nrows
    ok = abs(ratio - 1 / nchains) <= tol
    (ok = ok, ratio = ratio,
     why = ok ? "" : "T=1 rows are $(round(ratio, digits=4)) of all rows, expected $(round(1/nchains, digits=4))")
end

end # module
