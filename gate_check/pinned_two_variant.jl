# pinned_two_variant.jl -- pinned cells by depth, both bound variants, with the
# cells below the AEM sensitivity depth excluded from the counts.
#
#   julia --project=. gate_check/pinned_two_variant.jl --tag 3P \
#         --wide <store.jld2> --sed <store.jld2> --sens <sensdepth.jld2> \
#         --top 7.5 --base 23.5
#
# WHY BY DEPTH. A count of "cells on a bound" over all 65 cells mixes two
# different things. A shallow cell on a bound is a statement about the data: the
# fit wanted a conductivity outside the stated limits where the data had
# something to say. A cell at 250 m on a bound is a statement about nothing --
# the data cannot see it, so the regularisation puts it wherever the gradient
# points, and it lands on a limit. Only the shallow count bears on whether a
# node's acceptance owes anything to where the clamps were set.
#
# The split depth is the noise-weighted sensitivity depth from
# gate_check/sensitivity_depth.jl, recomputed here from the stored normalised
# sensitivity so that `frac` can be varied without another forward solve.
using JLD2, Statistics, Printf
const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) ? (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (i += 1)
    end
end
arg(k, d) = get(ARG, k, d)

const TAG  = arg("tag", "col")
const WTOP = parse(Float64, arg("top", "7.5"))
const WBOT = parse(Float64, arg("base", "23.5"))

S  = load(arg("sens", ""))
sn = S["sens_norm"]; CT = S["celltop"]; CB = S["cellbot"]

"Deepest cell whose normalised sensitivity reaches `f`; cells below it are the
ones excluded. This is the definition in sensitivity_depth.jl -- last crossing."
function split_at(f)
    idx = findall(>=(f), sn)
    isempty(idx) && return (0, 0.0)
    j = maximum(idx)
    (j, CB[j])
end

function load_variant(path)
    isempty(path) && return nothing
    isfile(path) || return nothing
    P = load(path)
    chi = haskey(P, "min_chi2") ? P["min_chi2"] : P["chi2"]
    (C = haskey(P, "conductance") ? P["conductance"] : (10 .^ P["lt"]) .* P["thick"],
     chi = chi, M = P["models"], mb = P["mbounds"], ndata = P["ndata"],
     chiref = P["chiref"], path = basename(path))
end

solved(m) = !isempty(m) && !any(isnan, m) && !all(x -> x == -2.0, m)
onbound(m, mb) = findall(x -> x <= mb[1] + 1e-9 || x >= mb[2] - 1e-9, m)

function report(V, nsplit, zsplit, f)
    tabs = V.ndata + 3.84
    @printf("\n  %s   bounds %s -> %.1e to %.1e S/m   (%s)\n",
            uppercase(arg("tag","")), string(V.mb), 10.0^V.mb[1], 10.0^V.mb[2], V.path)
    @printf("  absolute threshold %.2f ; split at %.1f m (cell %d), frac %.3f\n",
            tabs, zsplit, nsplit, f)
    @printf("  %-10s %10s %8s %8s %8s  %s\n",
            "C (S)", "min chi2", "pinned", "shallow", "deep", "shallow cells (depth range, m)")
    best = (C = NaN, chi = NaN)
    for i in eachindex(V.C)
        solved(V.M[i]) || continue
        isfinite(V.chi[i]) && V.chi[i] <= tabs || continue
        ob = onbound(V.M[i], V.mb)
        sh = filter(j -> j <= nsplit, ob); dp = filter(j -> j > nsplit, ob)
        desc = isempty(sh) ? "none" :
            join([@sprintf("%d(%.0f-%.0f)", j, CT[j], CB[j]) for j in sh], " ")
        @printf("  %-10.4f %10.3f %8d %8d %8d  %s\n",
                V.C[i], V.chi[i], length(ob), length(sh), length(dp), desc)
        isempty(sh) && (isnan(best.C) || V.C[i] > best.C) && (best = (C = V.C[i], chi = V.chi[i]))
    end
    if isnan(best.C)
        println("  largest accepted C with no SHALLOW cell on a bound: none")
    else
        @printf("  largest accepted C with no SHALLOW cell on a bound: %.4f S (chi2 %.3f)\n",
                best.C, best.chi)
    end
    best
end

function main()
    W = load_variant(arg("wide", "")); D = load_variant(arg("sed", ""))
    @printf("pinned_two_variant  tag=%s  window %.2f-%.2f m\n", TAG, WTOP, WBOT)
    println("sensitivity split, from ", basename(arg("sens","")), ":")
    for f in (0.01, 0.02, 0.05)
        j, z = split_at(f)
        @printf("  frac %.3f -> %6.1f m, cell %2d, %2d of 65 cells excluded as deep\n",
                f, z, j, 65 - j)
    end
    out = Dict{String,Any}()
    for f in (0.01, 0.02)
        nsplit, zsplit = split_at(f)
        @printf("\n================ frac = %.3f, split %.1f m ================\n", f, zsplit)
        W === nothing || (out["wide_f$(f)"] = report(W, nsplit, zsplit, f))
        D === nothing || (out["sed_f$(f)"]  = report(D, nsplit, zsplit, f))
    end
    jldsave(joinpath("gate_check","results","pinned_two_variant_$(TAG).jld2");
            tag = TAG, window = [WTOP, WBOT],
            best = Dict(k => (v.C, v.chi) for (k, v) in out),
            sens_file = arg("sens",""), wide = arg("wide",""), sed = arg("sed",""))
    println("\nwrote gate_check/results/pinned_two_variant_$(TAG).jld2")
end

main()
