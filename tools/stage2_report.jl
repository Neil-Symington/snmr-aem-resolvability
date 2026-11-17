# stage2_report.jl -- what Stage 2 achieved, per column, against the truth and
# against the TRUE noise. No inversion, no tuning; reads what is on disk.
#
#   julia --project=. tools/stage2_report.jl --tag 3S [--control] \
#         [--top 7.5] [--base 23.5]
#
# Four things, in the order Neil asked for them on 18 Nov:
#
#  1. THE NOISE THE INVERSION SETTLES ON, as a factor of the true noise, per
#     ensemble. The likelihood is hierarchical: s multiplies the stated sigma-hat
#     and is sampled, with a log-normal prior centred on s = 1. If s drifts above
#     1 the inversion is buying fit by widening the noise, which is where kernel
#     error would show up.
#
#  2. ACHIEVED DATA MISFIT AGAINST THE TRUE NOISE, chi2 per datum. The sampler's
#     own misfit is weighted by s*sigma-hat, so it cannot answer this; every
#     posterior draw is forwarded again here and weighted by sigma-hat itself.
#
#  3. KERNEL ERROR IN DATA SPACE. The data were generated through a kernel on the
#     TRUE conductivity; the inversion used one on the pooled Stage 1 median. The
#     true water content forwarded through the inversion's kernel, differenced
#     against the noise-free data, says how large that mismatch is in units of
#     the noise, before any inversion is involved. One forward call.
#
#  4. WINDOW WATER CONTENT AGAINST TRUTH.
#
# The control run (--control, kernel on the true conductivity) makes 3 identically
# zero by construction, which is the point of it: whatever 1, 2 and 4 then show is
# noise and non-uniqueness, not kernel error.
using JLD2, Statistics, Printf, DelimitedFiles
using SNMRForward, SMRPInversion
# The multi-chain rank-normalised split-Rhat and ESS of Vehtari et al. 2021,
# computed ACROSS the four independent ensembles' T=1 streams. The per-chain
# estimator this replaced was about 3x optimistic; see DECISIONS.md.
include(joinpath(@__DIR__, "..", "convergence_test", "prerun_lib.jl"))
using .PrerunLib: rhat_and_ess

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(ROOT, "data")
const RUN  = joinpath(OUT, "smr_runs")

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) && !startswith(ARGS[i+1], "--") ?
            (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (ARG[ARGS[i][3:end]] = "true"; i += 1)
    end
end
arg(k, d) = get(ARG, k, d)

const TAG     = arg("tag", "3S")
const CONTROL = haskey(ARG, "control")
const SUF     = CONTROL ? "_truekernel" : ""
# Window defaults are the columns' declared windows, the same ones the gate
# forces conductance over, so water and conductance are reported on one interval.
const WIN = Dict("3P" => (7.5, 23.5), "3S" => (7.5, 23.5),
                 "2P" => (17.5, 33.5), "2S" => (17.5, 33.5))
const WTOP = parse(Float64, arg("top",  string(WIN[TAG][1])))
const WBOT = parse(Float64, arg("base", string(WIN[TAG][2])))

"Integral of theta over [a,b] on the SNMR cell grid, by cell overlap."
function window_integral(θ, zsmr, a, b)
    top = vcat(0.0, zsmr[1:end-1]); bot = zsmr
    ov = [max(0.0, min(bot[i], b) - max(top[i], a)) for i in eachindex(zsmr)]
    (sum(θ .* ov), sum(ov))
end

"Mean theta over the window, per draw."
window_mean(Θ, zsmr, a, b) =
    (w = window_integral(ones(length(zsmr)), zsmr, a, b)[2];
     [window_integral(Θ[:, k], zsmr, a, b)[1] / w for k in axes(Θ, 2)])

"""
s samples for one ensemble at T = 1, as a FACTOR on sigma-hat.

THE TEMPERATURE SWAPS BETWEEN CHAINS. That is what parallel tempering does, and
it defeats every shortcut tried here first:

  - the nuisance VALUES file has three columns, CHAIN, SAMPLE, ln s -- reading
    the first two the wrong way round kept the hot chains and dropped T = 1
    entirely, and reported s around 2 to 4.6;
  - "chain 1 is the T = 1 chain" is false: every one of the four chains visits
    every one of the four temperatures (1.0, 1.357, 1.842, 2.5), and exactly a
    quarter of all rows are at T = 1;
  - "the T = 1 chain is the tightest in ln s" is also false, and the assertion
    that said so fired on two of eight ensembles.

The temperature is recorded per (chain, sample) in the nuisance MISFITS file,
whose five columns are CHAIN, SAMPLE, misfit, T, accept-rate. So the values are
joined to it on (chain, sample) and only T = 1 rows are kept. The first half of
the retained samples is dropped as burn-in, matching the stage's burninfrac.

The hot chains tolerate worse fits and so inflate s -- on 3S ensemble 1 the
chain-1 rows give 0.83 against 1.9 to 3.1 for chains 2 to 4 -- which is exactly
why the selection has to be right.
"""
function s_factor(prefix, e)
    fv = joinpath(RUN, "values_nuisance_$(prefix)_e$(e)_.bin")
    fm = joinpath(RUN, "misfits_$(prefix)_e$(e)_nuisance.bin")
    (isfile(fv) && isfile(fm)) || return Float64[]
    V = readdlm(fv); M = readdlm(fm)
    (size(V, 1) == 0 || size(M, 1) == 0) && return Float64[]
    size(M, 2) >= 4 || error("$(basename(fm)) has $(size(M,2)) columns, expected 5 " *
                             "(chain, sample, misfit, T, accept): cannot select T = 1")
    # (chain, sample) -> T
    Tof = Dict{Tuple{Int,Int},Float64}()
    for r in axes(M, 1)
        Tof[(Int(M[r, 1]), Int(M[r, 2]))] = Float64(M[r, 4])
    end
    ch = Int.(V[:, 1]); sm = Int.(V[:, 2]); val = Float64.(V[:, 3])
    at1 = [get(Tof, (ch[r], sm[r]), NaN) == 1.0 for r in eachindex(ch)]
    any(at1) || error("no nuisance sample at T = 1 in $(basename(fv))")
    smax = maximum(sm[at1])
    keep = at1 .& (sm .> smax / 2)
    exp.(val[keep])
end

function main()
    cf = joinpath(OUT, "cond_C2_typical$(TAG)T$(SUF).jld2")
    df = joinpath(OUT, "smr_data_murray$(TAG)T.jld2")
    for p in (cf, df)
        isfile(p) || error("missing $p")
    end
    C = load(cf); D = load(df); T = load(joinpath(OUT, "truth_murray$(TAG).jld2"))

    ens  = C["ens"]
    zsmr = collect(Float64, D["zsmr"])
    V    = collect(Float64, D["typical_V"])
    σ̂    = collect(Float64, D["typical_σ̂"])
    V0   = collect(Float64, D["V0_clean"])
    wtrue = collect(Float64, D["w_true"])
    N    = length(V)

    @printf("Stage 2 report  tag=%s%s\n", TAG, CONTROL ? "  CONTROL (true-conductivity kernel)" : "")
    @printf("  posterior : %s\n", basename(cf))
    @printf("  data      : %s\n", basename(df))
    @printf("  kernel    : %s\n", string(get(C, "kernel", "?")))
    haskey(C, "pooled_ndraws") && C["pooled_ndraws"] > 0 &&
        @printf("  pooled    : %d draws across %d folders\n",
                C["pooled_ndraws"], length(C["pooled_folders"]))
    @printf("  window    : %.2f-%.2f m ; %d data, sigma-hat median %.2f nV\n",
            WTOP, WBOT, N, 1e9 * median(σ̂))

    # ---- 3. kernel error in data space, before any inversion ---------------
    @printf("\n(3) KERNEL ERROR IN DATA SPACE  (true water content through the inversion's kernel)\n")
    if !haskey(C, "sigma_kernel_used")
        println("    sigma_kernel_used is absent from the posterior file: this run predates")
        println("    18 Nov. Cannot be computed without re-deriving the pooled median.")
    else
        σk = collect(Float64, C["sigma_kernel_used"]); zk = collect(Float64, C["z_kernel_used"])
        ctk = SNMRForward.ConductivityModel(σk, diff(zk))
        Fk  = SNMRForward.MRSForward_square(D["R"], zsmr, collect(Float64, D["qgrid"]),
                                            deg2rad(D["inc_deg"]), 0.0, D["Be"], ctk)
        Vk  = abs.(SNMRForward.forward(Fk, wtrue))
        dσ  = abs.(Vk .- V0) ./ σ̂
        @printf("    median %.2f sigma, max %.2f sigma, %.1f%% of %d data above 3 sigma\n",
                median(dσ), maximum(dσ), 100 * count(>(3.0), dσ) / N, N)
        @printf("    in amplitude: median %.2f%%, max %.2f%% of the noise-free data\n",
                100 * median(abs.(Vk .- V0) ./ V0), 100 * maximum(abs.(Vk .- V0) ./ V0))
    end

    # ---- 1, 2, 4 per ensemble ---------------------------------------------
    prefix = "C2_typical$(TAG)T"
    @printf("\n(1) INFERRED NOISE FACTOR s, and (2) chi2 PER DATUM against the TRUE noise\n")
    @printf("    %-4s %22s %28s\n", "ens", "s (factor on sigma-hat)", "chi2/N vs true noise")
    @printf("    %-4s %8s %6s %6s %10s %8s %8s\n", "", "p50", "p05", "p95", "p50", "p05", "p95")
    local Fi
    if haskey(C, "sigma_kernel_used")
        σk = collect(Float64, C["sigma_kernel_used"]); zk = collect(Float64, C["z_kernel_used"])
        ctk = SNMRForward.ConductivityModel(σk, diff(zk))
        Fi = SNMRForward.MRSForward_square(D["R"], zsmr, collect(Float64, D["qgrid"]),
                                           deg2rad(D["inc_deg"]), 0.0, D["Be"], ctk)
    else
        Fi = nothing
    end
    for e in eachindex(ens)
        sf = s_factor(prefix, e)
        Θ  = Float64.(ens[e])
        chi = Float64[]
        if Fi !== nothing && size(Θ, 2) > 0
            # subsample to keep this a report rather than a job
            cols = size(Θ, 2) <= 400 ? collect(axes(Θ, 2)) :
                   round.(Int, range(1, size(Θ, 2), length = 400))
            for k in cols
                Vp = abs.(SNMRForward.forward(Fi, Θ[:, k]))
                push!(chi, sum(((Vp .- V) ./ σ̂) .^ 2) / N)
            end
        end
        q(v, p) = isempty(v) ? NaN : quantile(v, p)
        @printf("    %-4d %8.3f %6.3f %6.3f %10.2f %8.2f %8.2f\n", e,
                q(sf, 0.50), q(sf, 0.05), q(sf, 0.95),
                q(chi, 0.50), q(chi, 0.05), q(chi, 0.95))
    end

    @printf("\n(4) WINDOW WATER CONTENT, mean theta over %.2f-%.2f m\n", WTOP, WBOT)
    θt = collect(Float64, T["θ_smrgrid"])
    wt = window_integral(θt, zsmr, WTOP, WBOT)[1] /
         window_integral(ones(length(zsmr)), zsmr, WTOP, WBOT)[2]
    @printf("    truth %.4f v/v\n", wt)
    @printf("    %-4s %8s %8s %8s %10s\n", "ens", "p05", "p50", "p95", "p50/truth")
    allm = Float64[]
    for e in eachindex(ens)
        m = window_mean(Float64.(ens[e]), zsmr, WTOP, WBOT)
        append!(allm, m)
        @printf("    %-4d %8.4f %8.4f %8.4f %9.3fx\n", e,
                quantile(m, 0.05), quantile(m, 0.50), quantile(m, 0.95),
                quantile(m, 0.50) / wt)
    end
    @printf("    pooled %7.4f %8.4f %8.4f %9.3fx   truth inside p05-p95: %s\n",
            quantile(allm, 0.05), quantile(allm, 0.50), quantile(allm, 0.95),
            quantile(allm, 0.50) / wt,
            (wt >= quantile(allm, 0.05) && wt <= quantile(allm, 0.95)) ? "YES" : "**NO**")

    # Convergence OF THE WINDOW AVERAGE, across the four ensembles' T=1 streams.
    # Rhat across the chains of one PT ensemble would be meaningless, they target
    # p^(1/T); across independent ensembles it is the right comparison.
    chains = [window_mean(Float64.(ens[e]), zsmr, WTOP, WBOT) for e in eachindex(ens)]
    r = rhat_and_ess(chains)
    @printf("\n    convergence of the window average, across %d ensembles:\n", length(chains))
    @printf("      R-hat %.3f   ESS bulk %.0f   ESS tail %.0f   (draws per ensemble: %s)\n",
            r.rhat, r.ess_bulk, r.ess_tail, join(length.(chains), ", "))
    @printf("      %s\n", (isfinite(r.rhat) && r.rhat <= 1.05 && r.ess_bulk >= 100) ?
            "PASS at R-hat <= 1.05 and ESS >= 100" :
            "**FAILS** the R-hat <= 1.05 / ESS >= 100 gate")

    # The spread across ensembles, which is what the gate is really about.
    meds = [median(c) for c in chains]
    @printf("      ensemble medians %s -> spread %.3fx\n",
            join([@sprintf("%.4f", m) for m in meds], ", "),
            maximum(meds) / minimum(meds))
end

main()
