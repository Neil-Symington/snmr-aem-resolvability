# boundary_sweep.jl -- how much of the package-conductance mismatch is the
# lower boundary sitting on a contrast?
#
#   julia --project=. convergence_test/boundary_sweep.jl
#
# Reads B1 output only. No inversions. Also checks whether the TRUE window
# conductance lies inside the gate's accepted range, which is the quantity the
# M5 decision path uses.
using JLD2, Statistics, Printf

const HERE = normpath(@__DIR__)
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(ROOT, "data")
const GADI_SRC = joinpath(ROOT, "src")
include(joinpath(ROOT, "src", "aem_setup.jl"))
using HiQGA.transD_GP

filesha(p) = bytes2hex(reinterpret(UInt8, [hash(read(p))]))
cell_edges(ztop, zcen) = (n = length(ztop);
    [i <= n ? ztop[i] : ztop[n] + 2*(zcen[n]-ztop[n]) for i in 1:n+1])
cond(s, e, a, b) = sum(max(0.0, min(e[i+1],b) - max(e[i],a)) * s[i] for i in eachindex(s))

"All T=1 draws of sigma(z) for one B1 ensemble, plus the cell edges."
function ensemble_sigma(tag, e)
    dir = joinpath(OUT, "aem$(tag)_b1e$(e)")
    isdir(dir) || return (nothing, nothing)
    S = load(joinpath(dir, "aem_setup.jld2"))
    znall = S["znall_aem"]; edges = cell_edges(S["zboundaries_aem"], S["zall_aem"])
    opt = make_aem_opt(znall, "hypsynth_aem_"; nmax = get(S, "nmax_k", 40))
    opt.xall[:] .= znall
    old = pwd(); cd(dir)
    F = try transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
        finally cd(old) end
    (10 .^ (-hcat([f[:] for f in F]...)), edges)
end

# column, window, package top, package bases to sweep, gate band in S (conductance)
CASES = [
 ("3P", (7.5, 23.5), 6.0, [36.0,38.0,40.0,42.0,44.0],
      (lo = 10^-2.7333333333333334 * 16, hi = 10^-1.6666666666666667 * 16,
       openlo = true, openhi = true, src = "stage4_synth_3P.jld2")),
 ("2P", (17.5, 33.5), 4.0, [30.0,32.0,34.0,36.0,38.0],
      (lo = 0.03192, hi = 7.079, openlo = true, openhi = false,
       src = "verdict_synth_2P.csv (1-2 Nov)")),
]

for (tag, win, ptop, bases, g) in CASES
    tf = joinpath(OUT, "truth_murray$tag.jld2")
    T = load(tf); z = T["z"]; sb = T["σ_bulk"]; dz = T["dz"]
    Ctruth(a,b) = sum(sb[i]*dz for i in eachindex(z) if a <= z[i] < b)

    println("="^86)
    @printf("%s   config %s   truth %s  hash %s\n", tag,
            tag == "3P" ? "synthetics/config/paper_3P.toml" : "synthetics/config/confined_aquifer.toml",
            basename(tf), filesha(tf))
    println("="^86)

    E = [ensemble_sigma(tag, e) for e in 1:4]
    ok = [x[1] !== nothing for x in E]
    @printf("ensembles found: %d\n\n", count(ok))

    @printf("package conductance, top fixed at %.1f m, base swept (S)\n", ptop)
    @printf("%-8s %9s | %-39s | %8s\n", "base", "truth", "ensemble p50 (1..4)", "mean/truth")
    for b in bases
        ct = Ctruth(ptop, b)
        ms = Float64[]
        for (σ, ed) in E
            σ === nothing && continue
            v = [cond(view(σ,:,j), ed, ptop, b) for j in axes(σ,2)]
            push!(ms, median(v))
        end
        @printf("%6.1f m %9.4f | %-39s | %7.2fx\n", b, ct,
                join([@sprintf("%.4f", m) for m in ms], "  "), mean(ms)/ct)
    end

    println()
    @printf("window %.1f-%.1f m: truth %.4f S\n", win..., Ctruth(win...))
    @printf("  gate accepted range %.4f - %.4f S  (lower edge %s, upper edge %s; %s)\n",
            g.lo, g.hi, g.openlo ? "OPEN" : "closed", g.openhi ? "OPEN" : "closed", g.src)
    tw = Ctruth(win...)
    @printf("  TRUE window conductance is %s the gate range\n",
            (tw >= g.lo && tw <= g.hi) ? "INSIDE" : "**OUTSIDE**")
    for (i,(σ,ed)) in enumerate(E)
        σ === nothing && continue
        v = [cond(view(σ,:,j), ed, win...) for j in axes(σ,2)]
        @printf("  ensemble %d sampled p05-p95 %.4f - %.4f S : truth %s\n", i,
                quantile(v,0.05), quantile(v,0.95),
                (tw >= quantile(v,0.05) && tw <= quantile(v,0.95)) ? "inside" : "**OUTSIDE**")
    end
    println()
end
