# test_gate_solver.jl -- Step 2's acceptance gate for the :map solver.
#
#   julia --project=. test_gate_solver.jl [--tags 3,2,3S,2S]
#
# These are assertions, not diagnostics. The point of the overnight failure was
# that a broken solver produced plausible-looking L(C) curves and nothing in the
# pipeline objected; reading monotonicity off a log by eye is how that happened.
#
# TEST A -- :map started from the converged :minimise model must not go
#   BACKWARDS. The brief specifies lambda2 -> 0 so that the MAP coincides with
#   the ML model and the test is an identity check; in practice lambda2 = 1e-6
#   leaves H with a condition number around 1e13 (42 data, 65 cells, so the
#   prior carries 23 directions on its own), the Gauss-Newton step is then
#   enormous, and the line search backtracks on every one of 120 steps. The test
#   costs ~1400 forward models and measures conditioning rather than the solver.
#   So it is run at lambda2 where H is well posed, with the identity check
#   replaced by the property it was standing in for: Phi(result) <= Phi(mref),
#   and at least one step accepted unless the start is already stationary.
#   That is what actually failed before -- zero steps accepted, mstart returned.
#
# TEST B -- unconstrained chi2 must be non-decreasing in lambda2 across the
#   ladder (more smoothing cannot improve a fit) AND must actually move between
#   adjacent nodes. chi2 identical to 12 figures across a 100x change in the
#   prior means the optimiser returned its starting model.
using JLD2, Statistics, Printf, LinearAlgebra
const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "occam_core.jl")); using .OccamCore
include(joinpath(ROOT, "src", "aem_jac.jl"))
const OUT = joinpath(HERE, "..", "data")

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) ? (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (i += 1)
    end
end
const TAGS = split(get(ARG, "tags", "3,2,3S,2S"), ',')
const LADDER = 10.0 .^ (0:6)
const RTOL_A = 1e-5          # drift tolerance on test A, in log10 sigma
const MOVE_B = 1e-6          # minimum relative chi2 change between ladder nodes

mutable struct Res; name::String; pass::Bool; msg::String; end
results = Res[]
note(n, p, m) = (push!(results, Res(n, p, m)); @printf("  [%s] %-28s %s\n", p ? "PASS" : "FAIL", n, m); flush(stdout))

for tag in TAGS
    @printf("\n%s\ntag %s\n%s\n", "="^70, tag, "="^70); flush(stdout)
    T = load(joinpath(OUT, "truth_murray$(tag).jld2"))
    B = build_aem_jac(T, 123); aem = B.aem
    n = length(aem.ρ) - 1; ndata = length(aem.res)
    R = makeR1(n); m0 = fill(-2.0, n)

    # ---- TEST A ----
    mref, χref = occam(aem, m0, m0, R, ndata; mode = :minimise)
    @printf("  :minimise reference  chi2 %.4f  (/N %.4f)\n", χref, χref/ndata); flush(stdout)
    pen(mm, l2) = l2 * (sum(abs2, R*mm) + 0.1 * sum(abs2, mm .- m0))
    for l2 in (1.0, 100.0)
        Φref = χref + pen(mref, l2)
        m, χ, _, nf = occam(aem, mref, m0, R, ndata; mode = :map, lambda2 = l2,
                            beta2 = 0.1, withH = true, info = true)
        Φ = χ + pen(m, l2)
        ok = Φ <= Φref * (1 + 1e-9) && (nf.nacc >= 1 || nf.gnorm < 1e-6)
        note("A lambda2=$(l2)", ok,
             @sprintf("phi %.4f vs %.4f, chi2 %.3f, %d accepted, |g_proj| %.2e%s",
                      Φ, Φref, χ, nf.nacc, nf.gnorm, nf.stuck ? "  <-- STUCK" : ""))
    end

    # ---- TEST B ----
    χs = Float64[]; stuck = 0; unsettled = 0; onbound = 0
    for l2 in LADDER
        _m, χ, _H, nf = occam(aem, mref, m0, R, ndata; mode = :map, lambda2 = l2,
                              beta2 = 0.1, withH = true, info = true)
        push!(χs, χ)
        nf.stuck && (stuck += 1); nf.settled || (unsettled += 1); nf.nbound > 0 && (onbound += 1)
        @printf("    lambda2 %-9.3g chi2 %12.4f  (/N %8.3f)  acc %3d  njac %3d  bound %2d  %s\n",
                l2, χ, χ/ndata, nf.nacc, nf.njac, nf.nbound,
                nf.settled ? "settled" : nf.stuck ? "STUCK" : nf.capped ? "step cap" :
                nf.budget ? "fev budget" : "unsettled"); flush(stdout)
    end
    note("B no stuck node", stuck == 0, "$stuck of $(length(LADDER)) returned the start model")
    # Reported, not asserted: a step-capped or boundary node is not a bug in the
    # solver, it is a node whose log ML must not be allowed to win (see §0.2).
    # stage4_gate.jl excludes them from lambda-hat; the count is what matters.
    note("B nodes settled", true,
         @sprintf("%d of %d settled, %d step-capped/unsettled, %d on the bound",
                  length(LADDER)-unsettled, length(LADDER), unsettled, onbound))
    d = diff(χs)
    note("B chi2 non-decreasing", all(d .>= -1e-6 .* max.(abs.(χs[1:end-1]), 1.0)),
         @sprintf("smallest step-to-step change %+.4g (negative = a decrease, which fails)",
                  minimum(d)))
    rel = abs.(d) ./ max.(abs.(χs[1:end-1]), 1.0)
    note("B chi2 moves per node", all(rel .> MOVE_B),
         @sprintf("smallest relative change %.3g", minimum(rel)))
end

@printf("\n%s\n%d of %d checks passed\n", "="^70, count(r -> r.pass, results), length(results))
for r in results; r.pass || @printf("  FAILED: %s -- %s\n", r.name, r.msg); end
exit(all(r -> r.pass, results) ? 0 : 1)
