#!/usr/bin/env julia
# naive_probability.jl -- salinity classes from AEM conductance ALONE.
#
#   julia --project=. tools/naive_probability.jl --cols 3P,3S,2P,2S
#   julia --project=. tools/naive_probability.jl --cols 2Pv2,2Sv2 --sampling both
#
# FOR THE DISCUSSION SECTION ONLY. This is NOT a workflow stage: nothing reads
# its output, no stage depends on it, and it must never be used to produce a
# result. It exists to show what the standard practice -- invert AEM, assume
# rock properties, divide out Archie -- would have concluded on these columns,
# so the paper can say what the joint method adds.
#
# It is READ-ONLY. It assembles existing Stage 1 draws and runs no inversion.
#
# WHAT IT DOES. For every post-burn-in T=1 draw it takes the mean conductivity
# over the column's stated window, then samples the rock properties:
#
#   sigma_bulk = sigma_w * theta^m + clay       (Equation 1)
#   sigma_w    = (sigma_eff - clay) / theta^m
#   TDS        = sigma_w * 6400                 ([priors].mg_per_sm)
#
# theta ~ U(0.25, 0.40) and m ~ U(1.50, 1.90) from [priors]. A NEGATIVE
# sigma_w means the clay term alone already exceeds the measured conductance;
# it is counted as fresh, which is the generous reading.
#
# THE CLAY PRIOR'S SHAPE DECIDES THE ANSWER, so both are offered and the
# default is the configured one. [priors.clay].product is [1e-4, 1.0] S/m,
# declared LOG-UNIFORM at config.toml:225-231.
#
#   --sampling loguniform  (default) the configured prior
#   --sampling uniform     uniform on the same range. DEGENERATE: the mean
#                          clay term is then 0.5 S/m, which exceeds the whole
#                          measured conductance of every column in this study
#                          -- 3P's is 0.0104 S/m -- so most draws go negative
#                          and the median salinity is 0 mg/L everywhere. It
#                          answers a question about the prior, not the water.
#   --sampling both        both, one table each.
#
# --scales multiplies the clay upper limit, to show how far the answer moves
# with that one assumption. Default "0.5,1,2".
using JLD2, Printf, Statistics, Random, TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
cd(ROOT)
include(joinpath(ROOT, "src", "aem_setup.jl"))

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) ?
            (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (i += 1)
    end
end
arg(k, d) = get(ARG, k, d)

# column -> (Stage 1 setting, window). The settings are tools/m4_run.sh's.
const COLDEF = Dict(
    "3P"   => ("m4c1", (7.5, 23.5)),  "2P"   => ("m4c2", (17.5, 33.5)),
    "3S"   => ("m4c3", (7.5, 23.5)),  "2S"   => ("m4c4", (17.5, 33.5)),
    "2Pv2" => ("m4c5", (17.5, 33.5)), "2Sv2" => ("m4c6", (17.5, 33.5)))

const COLS   = [strip(c) for c in split(arg("cols", "3P,3S,2P,2S"), ",")]
const NENS   = parse(Int, arg("nens", "4"))
const NREP   = parse(Int, arg("nrep", "200"))
const SEED   = parse(Int, arg("seed", "20261122"))
const SCALES = parse.(Float64, split(arg("scales", "0.5,1,2"), ","))
const WHICH  = arg("sampling", "loguniform")
WHICH in ("uniform", "loguniform", "both") ||
    error("--sampling must be uniform, loguniform or both; got $WHICH")

const CFG   = TOML.parsefile(joinpath(ROOT, "config.toml"))
const PHI   = Tuple(Float64.(CFG["priors"]["porosity_range"]))
const CLAY  = Tuple(Float64.(CFG["priors"]["clay"]["product"]))
const MGL   = Float64(CFG["priors"]["mg_per_sm"])
const MR    = (1.50, 1.90)        # [priors].m_lithology = "clayey_sand"

cell_edges(za) = begin
    e = zeros(length(za) + 1)
    for i in eachindex(za)
        e[i+1] = i < length(za) ? (za[i] + za[i+1]) / 2 : za[i] + (za[i] - za[i-1]) / 2
    end
    e
end

"Mean conductivity over [a,b], per draw, from one ensemble folder."
function sigma_eff(dir, a, b)
    S = load(joinpath(dir, "aem_setup.jld2"))
    zb = cell_edges(S["zall_aem"])
    opt = make_aem_opt(S["znall_aem"], "hypsynth_aem_";
                       nmax = haskey(S, "nmax_k") ? S["nmax_k"] : 40)
    opt.xall[:] .= S["znall_aem"]
    old = pwd(); cd(dir)
    F = try transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
        finally cd(old) end
    σ = 10 .^ (-hcat([f[:] for f in F]...))
    map(axes(σ, 2)) do j
        t = 0.0
        for i in axes(σ, 1)
            ov = min(zb[i+1], b) - max(zb[i], a); ov > 0 && (t += σ[i, j] * ov)
        end
        t / (b - a)
    end
end

"Class fractions and median TDS for a vector of sigma_eff."
function classify(seff, clay, logc; nrep = NREP)
    lo, hi = clay; llo, lhi = log10(lo), log10(hi)
    tds = Vector{Float64}(undef, length(seff) * nrep)
    k = 0
    for s in seff, _ in 1:nrep
        φ  = PHI[1] + rand() * (PHI[2] - PHI[1])
        m  = MR[1]  + rand() * (MR[2]  - MR[1])
        cl = logc ? 10^(llo + rand() * (lhi - llo)) : lo + rand() * (hi - lo)
        σw = (s - cl) / φ^m
        tds[k += 1] = σw <= 0 ? 0.0 : σw * MGL
    end
    (fresh  = count(<(1500.0), tds) / length(tds),
     brack  = count(x -> 1500.0 <= x < 3000.0, tds) / length(tds),
     saline = count(>=(3000.0), tds) / length(tds),
     med    = median(tds))
end

function main()
    Random.seed!(SEED)
    println("NAIVE SALINITY FROM AEM ALONE -- discussion only, not a workflow stage")
    @printf("theta U(%.2f,%.2f)  m U(%.2f,%.2f)  clay product %.0e-%.1f S/m  TDS = sigma_w x %.0f\n",
            PHI..., MR..., CLAY..., MGL)
    println("negative sigma_w counted as fresh; seed $SEED, $NREP property draws per model draw\n")

    for col in COLS
        haskey(COLDEF, col) || (@printf("%s: not a known column\n", col); continue)
        setting, (a, b) = COLDEF[col]
        runs = Vector{Vector{Float64}}()
        for e in 1:NENS
            d = joinpath(ROOT, "data", "aem$(col)_$(setting)e$(e)")
            isdir(d) || (@printf("  %s e%d: MISSING (%s)\n", col, e, basename(d)); continue)
            push!(runs, sigma_eff(d, a, b))
        end
        isempty(runs) && continue
        all = vcat(runs...)
        @printf("=== %-5s window %.1f-%.1f m   %d draws from %d run(s)   sigma_eff median %.5f S/m\n",
                col, a, b, length(all), length(runs), median(all))
        for (nm, logc) in (WHICH == "both" ? (("uniform", false), ("log-uniform", true)) :
                           WHICH == "uniform" ? (("uniform", false),) : (("log-uniform", true),))
            for sc in SCALES
                r = classify(all, (CLAY[1], CLAY[2] * sc), logc)
                @printf("  %-11s clay upper %5.2f   fresh %5.1f%%  brackish %5.1f%%  saline %5.1f%%  median %7.0f mg/L\n",
                        nm, CLAY[2] * sc, 100r.fresh, 100r.brack, 100r.saline, r.med)
            end
            @printf("  %-11s per run, clay upper %.2f: ", nm, CLAY[2])
            for (i, v) in enumerate(runs)
                r = classify(v, CLAY, logc)
                @printf("e%d %4.1f/%4.1f/%4.1f  ", i, 100r.fresh, 100r.brack, 100r.saline)
            end
            println("(fresh/brackish/saline %)")
        end
        println()
    end
end
main()
