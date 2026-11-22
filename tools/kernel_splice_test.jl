#!/usr/bin/env julia
# kernel_splice_test.jl -- which part of the conductivity error moves the SNMR
# response? A SAMPLER-FREE forward test.
#
#   julia --project=. tools/kernel_splice_test.jl --cols 3P,2P
#
# !!! NOT VALIDATED -- DO NOT USE ITS NUMBERS. 23 Nov 2026.
#
# It fails its own consistency check. Case (b), the pooled kernel everywhere,
# should reproduce the bias Stage 2 actually shows: on 2P the pooled run
# overstates window water content by 9.6% (0.2672 against 0.2438), so (b)
# ought to return a scale factor near 1.10. It returns 0.9970. A linear
# forward test with a least-squares projection onto the window's own response
# therefore cannot decompose this bias, and the per-splice numbers below mean
# nothing. It also failed to reproduce the earlier v1 summary figures, but
# that is the lesser evidence: the internal check is the decisive one.
#
# Superseded by running Stage 2 itself on spliced kernels, which needs no
# code change to the stage -- the default kernel path reads zall_aem and
# sigma50 from the file --label points at.
#
# ANALYSIS SCRIPT, NOT A WORKFLOW STAGE. Nothing reads its output and no stage
# depends on it. It runs no inversion: it forward-models the TRUE water content
# through conductivity models assembled from the truth and from the pooled
# Stage 1 median, and compares amplitudes.
#
# THE FIVE CASES. Water content is the truth's theta_SMR throughout; only the
# conductivity differs.
#
#   (a) true       the true sigma_bulk on the kernel grid (sigma_gen)
#   (b) pooled     the pooled Stage 1 median (sigma_kernel_used)
#   (c) above      pooled ABOVE the window, truth elsewhere
#   (d) window     pooled INSIDE the window, truth elsewhere
#   (e) below      pooled BELOW the window, truth elsewhere
#
# THE SCALE FACTOR. For each case we report the single multiplier on the
# WINDOW's water content that best reproduces that case's amplitudes using the
# TRUE conductivity -- i.e. how much apparent water the kernel error is worth.
# The SNMR response is LINEAR in water content at fixed kernel, so
#
#     V(theta + (s-1) * theta_win, sigma_true) = V_a + (s-1) * V_win
#
# and the least-squares s has a closed form; no search is needed. V_win is the
# response of the window's water alone, through the true conductivity.
using JLD2, Printf, Statistics, SNMRForward

const ROOT = normpath(joinpath(@__DIR__, ".."))
cd(ROOT)

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) ?
            (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (i += 1)
    end
end
arg(k, d) = get(ARG, k, d)

# column -> (SNMR data tag, pooled Stage 2 folder prefix, window)
const COLDEF = Dict(
  "3P"   => ("3PT",   "3PT_pooled",   (7.5, 23.5)),
  "3S"   => ("3ST",   "3ST_pooled",   (7.5, 23.5)),
  "2P"   => ("2PT",   "2PT_pooled",   (17.5, 33.5)),
  "2S"   => ("2ST",   "2ST_pooled",   (17.5, 33.5)),
  "2Pv2" => ("2Pv2T", "2Pv2T_pooled", (17.5, 33.5)),
  "2Sv2" => ("2Sv2T", "2Sv2T_pooled", (17.5, 33.5)))

"Newest Stage 2 folder matching `pre` that holds a .jld2."
function newest(pre)
    for d in sort(readdir(joinpath(ROOT, "data", "stage2")); rev = true)
        startswith(d, pre) || continue
        p = joinpath(ROOT, "data", "stage2", d)
        fs = filter(x -> endswith(x, ".jld2"), readdir(p))
        isempty(fs) || return joinpath(p, fs[1])
    end
    ""
end

"Layer thicknesses from layer centres/nodes, last layer the half-space."
thick_of(z) = diff(z)

function run_column(col)
    tag, pre, (a, b) = COLDEF[col]
    D = load(joinpath(ROOT, "data", "smr_data_murray$(tag).jld2"))
    S2 = newest(pre)
    isempty(S2) && (@printf("%-5s no pooled Stage 2 run on disk\n", col); return nothing)
    P = load(S2)

    zsmr = collect(Float64, D["zsmr"])
    w    = collect(Float64, D["w_true"])            # the TRUE water content
    qg   = collect(Float64, D["qgrid"])
    Be, inc, R = Float64(D["Be"]), Float64(D["inc_deg"]), Float64(D["R"])
    ϕ = deg2rad(inc)

    zk  = collect(Float64, D["z_gen"])              # kernel depth grid, 65 cells
    σt  = collect(Float64, D["sigma_gen"])          # TRUE sigma on that grid
    tk  = collect(Float64, D["t_gen"])
    σp  = collect(Float64, P["sigma_kernel_used"])  # pooled Stage 1 median
    length(σp) == length(σt) ||
        error("$col: pooled kernel has $(length(σp)) cells, truth has $(length(σt))")

    splice(mask) = [mask[i] ? σp[i] : σt[i] for i in eachindex(σt)]
    above  = zk .<  a
    inside = (zk .>= a) .& (zk .< b)
    below  = zk .>= b

    V(σ, water) = abs.(SNMRForward.forward(
        SNMRForward.MRSForward_square(R, zsmr, qg, ϕ, 0.0, Be,
            SNMRForward.ConductivityModel(σ, tk)), water))

    # water of the window alone, on the SNMR layer grid (bottoms at zsmr)
    top = vcat(0.0, zsmr[1:end-1])
    frac = [max(0.0, min(zsmr[i], b) - max(top[i], a)) / (zsmr[i] - top[i])
            for i in eachindex(zsmr)]
    wwin = w .* frac

    Va   = V(σt, w)
    Vwin = V(σt, wwin)          # linear response of the window's water alone
    denom = sum(Vwin .^ 2)

    cases = (("(b) pooled everywhere", σp),
             ("(c) pooled above only", splice(above)),
             ("(d) pooled window only", splice(inside)),
             ("(e) pooled below only", splice(below)))

    @printf("\n=== %-5s window %.1f-%.1f m, %d pulse moments, %d kernel cells\n",
            col, a, b, length(qg), length(σt))
    @printf("  %-24s %-9s %-9s %-9s  %s\n", "case", "min ratio", "max ratio", "mean", "window-water scale")
    out = Dict{String,Float64}()
    for (nm, σc) in cases
        Vc = V(σc, w)
        r  = Vc ./ Va
        s  = 1.0 + sum((Vc .- Va) .* Vwin) / denom
        @printf("  %-24s %-9.4f %-9.4f %-9.4f  %.4f\n", nm, minimum(r), maximum(r), mean(r), s)
        out[nm] = s
    end
    out
end

function main()
    cols = [strip(c) for c in split(arg("cols", "3P,3S,2P,2S,2Pv2,2Sv2"), ",")]
    println("SNMR forward splice test -- analysis only, no inversion, no stage changes")
    println("water content is the TRUTH in every case; only the conductivity differs")
    for c in cols
        haskey(COLDEF, c) ? run_column(c) : @printf("%s: unknown column\n", c)
    end
end
main()
