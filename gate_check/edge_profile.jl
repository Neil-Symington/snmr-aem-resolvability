# edge_profile.jl -- M1.1. chi2 against FORCED window conductance, continuing
# past both reported edges, for test one (3P) and its saline control (3S).
#
#   julia --project=. gate_check/edge_profile.jl
#
# Reads stored Stage 4 output only; runs no solves. Writes a table and the data
# for the figure into gate_check/results/.
#
# The question: is each reported edge a DATA limit (chi2 crosses chiref + 3.84
# inside the search range) or a SEARCH/SOLVER limit (the profile is still below
# the threshold where the search stops, or the node was dropped for a solver
# reason)? The two must be labelled differently and are not at present.
using JLD2, Printf

const P = "/home/symo/Documents/GA/AEM_SMR_Resolvability_Quickrun/outputs"
const HERE = normpath(@__DIR__)
mkpath(joinpath(HERE, "results"))

"Nodes, chi2, solver flags for one column, from whichever record the file has."
function profile(tag)
    D = load(joinpath(P, "synth_$tag", "stage4_synth_$tag.jld2"))
    thr = D["chiref"] + 3.84
    if haskey(D, "coarse")                      # two-pass: the coarse sweep
        C = D["coarse"]
        (lt = C["ltc"], chi2 = C["chi2"], good = C["good"], nbound = C["nbound"],
         thr = thr, chiref = D["chiref"], L = D["thick"], D = D, which = "coarse grid")
    else                                        # single-pass: the lambda-hat sweep
        lam = string(D["lambda_hat"]); S = D["sweeps"]
        k = haskey(S, lam) ? lam : first(keys(S)); sw = S[k]
        (lt = D["lt"], chi2 = sw.chi2, good = sw.good, nbound = sw.nbound,
         thr = thr, chiref = D["chiref"], L = D["thick"], D = D, which = "sweep lambda2=$k")
    end
end

io = open(joinpath(HERE, "results", "edge_profile.md"), "w")
emit(s="") = (println(s); println(io, s))
emit("# M1.1 — is each gate edge a data limit?")
emit()
emit("chi2 against forced window conductance, continued past both reported edges.")
emit("Threshold is chiref + 3.84. Read only from stored Stage 4 output; no solves.")
emit()

csv = open(joinpath(HERE, "results", "edge_profile.csv"), "w")
println(csv, "tag,log10_sigma,sigma_Spm,conductance_S,chi2,threshold,nbound,good")

for tag in ("3P", "3S")
    p = profile(tag); D = p.D
    lo_rep = haskey(D,"edge_lo_dex") ? 10^D["edge_lo_dex"]*p.L : NaN
    hi_rep = haskey(D,"edge_hi_dex") ? 10^D["edge_hi_dex"]*p.L : NaN
    emit("## $tag  ($(p.which))")
    emit()
    emit(@sprintf("chiref %.3f, threshold %.3f, window thickness %.1f m, sigma search bounds %s S/m.",
                  p.chiref, p.thr, p.L, string(D["sigma_bounds"])))
    isnan(lo_rep) || emit(@sprintf("Reported band %.4f - %.4f S (edge_solves = %s).",
                                   lo_rep, hi_rep, string(get(D,"edge_solves","-"))))
    emit()
    emit("| C (S) | sigma (S/m) | chi2 | vs threshold | nbound | solver good |")
    emit("|---|---|---|---|---|---|")
    acc = Int[]
    for i in eachindex(p.lt)
        C = 10^p.lt[i]*p.L
        below = p.chi2[i] <= p.thr
        below && p.good[i] && push!(acc, i)
        @printf(csv, "%s,%.4f,%.6f,%.6f,%.4f,%.4f,%d,%d\n", tag, p.lt[i], 10^p.lt[i], C,
                p.chi2[i], p.thr, p.nbound[i], p.good[i])
        # print a readable subset: everything for 3P, every other node for 3S
        (tag == "3P" || i <= 22) &&
          emit(@sprintf("| %.4f | %.5f | %.2f | %s | %d | %s |", C, 10^p.lt[i], p.chi2[i],
                        below ? "**below**" : "above", p.nbound[i], p.good[i] ? "yes" : "**NO**"))
    end
    emit()
    if isempty(acc)
        emit("No node is both below threshold and solver-good.")
    else
        emit(@sprintf("Nodes below threshold AND solver-good: C = %.4f to %.4f S.",
                      10^p.lt[first(acc)]*p.L, 10^p.lt[last(acc)]*p.L))
    end
    # is the lowest node below threshold? then the low edge is NOT a data limit
    emit(@sprintf("At the LOWEST search node (%.4f S) chi2 = %.2f, which is %s the threshold -> low edge is %s.",
                  10^p.lt[1]*p.L, p.chi2[1], p.chi2[1] <= p.thr ? "BELOW" : "above",
                  p.chi2[1] <= p.thr ? "**a search limit, not a data limit**" : "**a data limit**"))
    emit(@sprintf("At the HIGHEST search node (%.4f S) chi2 = %.2f, which is %s the threshold -> high edge is %s.",
                  10^p.lt[end]*p.L, p.chi2[end], p.chi2[end] <= p.thr ? "BELOW" : "above",
                  p.chi2[end] <= p.thr ? "**a search limit, not a data limit**" : "**a data limit**"))
    emit()
end
close(csv); close(io)
println("\nwrote ", joinpath(HERE,"results","edge_profile.md"), " and .csv")
