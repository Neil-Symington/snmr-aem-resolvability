# make_report.jl -- build the one-page 3P report: AEM ensemble, SNMR ensemble
# and the gate's models on a shared depth axis, with the misfit profile below.
#
#   julia --project=. tools/make_report.jl \
#       --aem   data/aem3P_b1e1,data/aem3P_b1e2,data/aem3P_b1e3,data/aem3P_b1e4 \
#       --snmr  ../synthetics/outputs/HYP_SYNTH/cond_C2_typical3PT.jld2 \
#       --gate  gate_check/results/profile_3P_sed.jld2,gate_check/results/edges_3P_sed2.jld2
#
# Defaults for all three are the sources that exist on this machine today, so
# the bare command reproduces the page. Point --aem at the M4 folders and
# --snmr at the new Stage 2 output when they land and the page is rebuilt with
# one command.
#
# THIS RUNS NO INVERSION AND NO SOLVE. It re-assembles stored chains, bins
# them, reads the gate's stored per-node models, and writes one self-contained
# HTML file. Nothing is written outside --out.
#
# THE HISTOGRAMS ARE BINNED HERE, NOT IN THE BROWSER. Eight thousand AEM draws
# and twenty thousand SNMR draws are not going into a web page; the page gets
# per-cell counts and the percentile lines computed from the full set.
#
# WHY THE PERCENTILES ARE COMPUTED FROM THE POOLED DRAWS AND NOT READ FROM
# aem_posterior.jld2: that file keeps 30 subsampled draws per folder plus its
# own sigma10/50/90, which are per-folder. Pooling four folders changes the
# percentiles, so they are recomputed over the pooled set -- the same quantity
# the histogram shows. `tools/extract_full_draws.jl` is where the assembleTat1
# idiom comes from; stage2_snmr.jl's `pooled_sigma` is the same call.
#
# THE AIR LAYER. `sigma_aem_truth` carries 66 values and the first is the air
# layer the forward model owns (src/aem_setup.jl, and config.toml's note on
# z_max). The truth on the 65-cell grid is therefore `[2:end]`. Getting this
# wrong puts the whole truth profile one cell out of register with the
# ensemble, which looks like a small depth shift rather than an error.
#
# CELLS ON A STATED LIMIT are recomputed from each stored model with the same
# predicate `gate_check/pinned_two_variant.jl:55` uses. The `nbound` field in
# the edge store does not agree with the models it sits beside (it reads 0 at
# C = 0.0506 S where the model has three cells on the floor), so the model is
# taken as the record and the stored field is ignored. The reproduction of
# pinned_two_variant.jl's published counts is asserted below.
using HiQGA, JLD2, Printf, Statistics, Dates
using HiQGA.transD_GP

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "aem_setup.jl"))

# ---------------------------------------------------------------- arguments
const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        if startswith(ARGS[i], "--") && i < length(ARGS) && !startswith(ARGS[i+1], "--")
            ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2
        elseif startswith(ARGS[i], "--")
            ARG[ARGS[i][3:end]] = "true"; i += 1
        else
            i += 1
        end
    end
end
arg(k, d) = get(ARG, k, d)
splitlist(s) = [strip(x) for x in split(s, ",") if !isempty(strip(x))]
abspath_rel(p) = isabspath(p) ? p : normpath(joinpath(ROOT, p))

const AEMDIRS = abspath_rel.(splitlist(arg("aem",
    "data/aem3P_b1e1,data/aem3P_b1e2,data/aem3P_b1e3,data/aem3P_b1e4")))
const SNMRF   = abspath_rel(arg("snmr",
    "../synthetics/outputs/HYP_SYNTH/cond_C2_typical3PT.jld2"))
const GATEF   = abspath_rel.(splitlist(arg("gate",
    "gate_check/results/profile_3P_sed.jld2,gate_check/results/edges_3P_sed2.jld2")))
const TRUTHF  = abspath_rel(arg("truth", "data/truth_murray3P.jld2"))
const SENSF   = abspath_rel(arg("sens", "gate_check/results/sensdepth_3P.jld2"))
const CFGF    = abspath_rel(arg("config", "configs/synth_3P.toml"))
const SMOOTHF = abspath_rel(arg("smooth", "gate_check/results/smooth_3P_sed.jld2"))
const DLIM    = arg("depthlimit", "")
const SMRDATA = arg("smrdata", "data/smr_data_murray3PT.jld2")
const TIER    = arg("tier", "")
const OUTF    = abspath_rel(arg("out", "figures/report_3P.html"))
const ZMAX    = parse(Float64, arg("zmax", "100"))
const WTOP    = parse(Float64, arg("top", "7.5"))
const WBOT    = parse(Float64, arg("base", "23.5"))
const SNMRNOTE = arg("snmr-note", "superseded SNMR data, placeholder")
const FORCE   = haskey(ARG, "force")
# --fragment emits body content only, for a host that supplies its own
# <html>/<head> wrapper. The DEFAULT is a complete standalone document,
# because the point of this file is that it can be sent to someone.
const FRAGMENT = haskey(ARG, "fragment")

isfile(OUTF) && !FORCE && error("""
    $OUTF exists. Pass --force to overwrite it.
    Refusing by default: this script's only output is that one file, and a
    rebuild from different sources should not silently replace a page whose
    provenance line says something else.""")

# ------------------------------------------------------------------ helpers
pct(v, p) = quantile(v, p)
"Bin `vals` into `nb` equal bins on [lo,hi]; values outside are clamped in."
function binrow(vals, lo, hi, nb)
    c = zeros(Int, nb); w = (hi - lo) / nb
    @inbounds for x in vals
        isfinite(x) || continue
        j = clamp(floor(Int, (x - lo) / w) + 1, 1, nb)
        c[j] += 1
    end
    c
end
jnum(x) = isfinite(x) ? (@sprintf("%.6g", x)) : "null"
jarr(v) = "[" * join(jnum.(v), ",") * "]"
jiarr(v) = "[" * join(string.(v), ",") * "]"
jstr(s) = "\"" * replace(string(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n") * "\""

println("make_report.jl  ", Dates.format(now(), "yyyy-mm-dd HH:MM"))

# =========================================================== 1. AEM ensemble
println("\nAEM ensemble")

"""
Pool the post-burn-in T = 1 draws of several Stage 1 folders.

In a function, not at top level, and deliberately: `docs/M2_NOTES.md` records
four soft-scope failures in M1 from exactly this loop shape. At top level
`zall = ...` inside the loop binds a new local each iteration, which here
fails loudly and elsewhere has failed silently.
"""
function pool_aem(dirs)
    mats = Matrix{Float64}[]; per = Int[]
    zall = nothing; zbnd = nothing; σtruth = nothing
    for d in dirs
        isdir(d) || error("AEM folder not found: $d")
        S = load(joinpath(d, "aem_setup.jld2"))
        znall = S["znall_aem"]
        if zall === nothing
            zall = collect(Float64, S["zall_aem"])
            zbnd = collect(Float64, S["zboundaries_aem"])
            t = collect(Float64, S["σ_aem_truth"])
            # drop the air layer; see the header note
            σtruth = length(t) == length(zall) + 1 ? t[2:end] : t
        end
        opt = make_aem_opt(znall, "hypsynth_aem_"; nmax = get(S, "nmax_k", 40))
        opt.xall[:] .= znall
        old = pwd(); cd(d)
        F = try
            transD_GP.assembleTat1(opt, :fstar; temperaturenum = 1, burninfrac = 0.5)
        finally
            cd(old)
        end
        m = 10 .^ (-hcat([f[:] for f in F]...))      # fstar is log10(rho)
        push!(mats, m); push!(per, size(m, 2))
        @printf("  %-28s %5d T=1 draws\n", basename(d), size(m, 2))
    end
    (σ = hcat(mats...), zall = zall, zbnd = zbnd, truth = σtruth, per = per)
end

const AEMPOOL = pool_aem(AEMDIRS)
σ_aem = AEMPOOL.σ
zall = AEMPOOL.zall
zbnd = AEMPOOL.zbnd
σtruth_aem = AEMPOOL.truth
const NZ = size(σ_aem, 1)
@printf("  pooled %d draws over %d cells, sigma %.3g to %.3g S/m\n",
        size(σ_aem, 2), NZ, minimum(σ_aem), maximum(σ_aem))

# CELL GEOMETRY. The gate defines its cells as `celltop = [0; zall[1:end-1]]`,
# `cellbot = zall` (gate_check/profile_minchi2.jl:66) -- that is, `zall_aem` is
# read as cell BOTTOMS, not as centres, and NOT as `zboundaries_aem`, which is
# the forward model's layer grid and a different set of depths. The same
# convention is what `pinned_two_variant.jl` prints its depth ranges against:
# cells 1-5 come out as 0-1, 1-2, 2-4, 4-5, 5-7 m, which reproduces
# PINNED_TWO_VARIANT.md exactly. Use zboundaries_aem here instead and the
# window picks up the wrong cells and every profile is drawn half a cell out.
celltop = vcat(0.0, zall[1:end-1])
cellbot = collect(Float64, zall)
dz_aem = cellbot .- celltop
winw = [max(0.0, min(cellbot[i], WBOT) - max(celltop[i], WTOP)) for i in 1:NZ]
Ctruth_grid = sum(σtruth_aem .* winw)
@printf("  window %.1f-%.1f m: support %.4f m over %d of %d cells\n",
        WTOP, WBOT, sum(winw), count(>(1e-9), winw), NZ)
@printf("  truth on this grid: window conductance %.4f S\n", Ctruth_grid)

const AEM_LO, AEM_HI, AEM_NB = -3.3, 0.6, 78
aem_hist = [binrow(log10.(view(σ_aem, i, :)), AEM_LO, AEM_HI, AEM_NB) for i in 1:NZ]
aem_p05 = [pct(view(σ_aem, i, :), 0.05) for i in 1:NZ]
aem_p50 = [pct(view(σ_aem, i, :), 0.50) for i in 1:NZ]
aem_p95 = [pct(view(σ_aem, i, :), 0.95) for i in 1:NZ]

# ========================================================== 2. SNMR ensemble
println("\nSNMR ensemble")
isfile(SNMRF) || error("SNMR file not found: $SNMRF")
SN = load(SNMRF)
haskey(SN, "ens") || error("$SNMRF carries no `ens`; this is not a Stage 2 output")
snmr_ens = SN["ens"]
zsmr = collect(Float64, SN["zsmr"])
θ_smr = hcat([Float64.(M) for M in snmr_ens]...)
@printf("  %-28s %d ensembles, %d draws over %d depths\n",
        basename(SNMRF), length(snmr_ens), size(θ_smr, 2), size(θ_smr, 1))
if size(θ_smr, 2) < 50
    @warn "only $(size(θ_smr,2)) SNMR draws: this file is a point estimate, not a posterior, and its histogram will be a single line"
end
const NS = size(θ_smr, 1)
const SNMR_LO, SNMR_HI, SNMR_NB = 0.0, 0.45, 60
snmr_hist = [binrow(view(θ_smr, i, :), SNMR_LO, SNMR_HI, SNMR_NB) for i in 1:NS]
snmr_p05 = [pct(view(θ_smr, i, :), 0.05) for i in 1:NS]
snmr_p50 = [pct(view(θ_smr, i, :), 0.50) for i in 1:NS]
snmr_p95 = [pct(view(θ_smr, i, :), 0.95) for i in 1:NS]

# ---- the true SNMR water content, on the ensemble's own depth grid
TR = load(TRUTHF)
θtrue_smr = if haskey(TR, "θ_smrgrid") && length(TR["θ_smrgrid"]) == NS
    collect(Float64, TR["θ_smrgrid"])
else
    zt = collect(Float64, TR["zsmr"]); vt = collect(Float64, TR["θ_smrgrid"])
    [vt[argmin(abs.(zt .- z))] for z in zsmr]
end

# ============================================================== 3. the gate
println("\nGate")
solved(m) = !isempty(m) && all(isfinite, m) && !any(isnan, m) && !all(x -> x == -2.0, m)
onbound(m, mb) = findall(x -> x <= mb[1] + 1e-9 || x >= mb[2] - 1e-9, m)

"""
Merge the gate stores into one node set.

A node present in both stores keeps the LOWER misfit and the model that goes
with it: a sweep returns an upper bound on min chi2, not min chi2, so the
lowest solve ever obtained at a node is the one that stands
(`gate_check/results/EDGE_2P.md` records the same rule).
"""
function merge_gate(paths)
    nodes = Dict{Float64,NamedTuple}()
    mbounds = nothing; ndata = nothing; chiref = Inf; thick = nothing
    unsolved = Tuple{Float64,String}[]
    for p in paths
        isfile(p) || error("gate store not found: $p")
        P = load(p)
        chi = haskey(P, "min_chi2") ? P["min_chi2"] : P["chi2"]
        C   = haskey(P, "conductance") ? P["conductance"] : (10 .^ P["lt"]) .* P["thick"]
        mbounds === nothing && (mbounds = collect(Float64, P["mbounds"]))
        collect(Float64, P["mbounds"]) == mbounds ||
            error("gate stores disagree on parameter bounds: $(P["mbounds"]) vs $mbounds")
        ndata === nothing && (ndata = Int(P["ndata"]))
        thick === nothing && (thick = Float64(P["thick"]))
        chiref = min(chiref, Float64(P["chiref"]))
        for i in eachindex(C)
            m = collect(Float64, P["models"][i])
            if !solved(m) || !isfinite(chi[i])
                push!(unsolved, (C[i], basename(p))); continue
            end
            k = round(C[i], digits = 6)
            if !haskey(nodes, k) || chi[i] < nodes[k].chi2
                nodes[k] = (chi2 = Float64(chi[i]), model = m, src = basename(p))
            end
        end
        @printf("  %-28s %d nodes\n", basename(p), length(C))
    end
    (nodes = nodes, mbounds = mbounds, ndata = ndata, chiref = chiref,
     thick = thick, unsolved = unsolved)
end

const G = merge_gate(GATEF)
nodes = G.nodes; mbounds = G.mbounds; ndata = G.ndata
chiref = G.chiref; thick = G.thick; unsolved = G.unsolved
const TABS = ndata + 3.84
Cs = sort(collect(keys(nodes)))
@printf("  bounds %.1f to %.5f (log10 sigma) -> %.3g to %.3g S/m\n",
        mbounds[1], mbounds[2], 10^mbounds[1], 10^mbounds[2])
@printf("  N = %d data, absolute threshold chi2 <= %.2f\n", ndata, TABS)

if !isempty(unsolved)
    println("\n  NODES WITH NO STORED MODEL -- no solve met the conductance constraint there.")
    println("  They are not part of the admissible range; nothing is drawn for them.")
    for (c, s) in sort(unsolved, by = first)
        @printf("    C = %8.4f S   (%s)\n", c, s)
    end
end

adm = [c for c in Cs if nodes[c].chi2 <= TABS]
isempty(adm) && error("no admissible node has a stored model")
bestC = Cs[argmin([nodes[c].chi2 for c in Cs])]
const CHIBEST = nodes[bestC].chi2
# chiref = min(unconstrained multi-start, min over constrained nodes) -- the
# definition in DECISIONS.md. The stores stopped half way: profile_3P_sed took
# the min over ITS OWN twelve nodes, and edges_3P_sed2 was seeded with that
# number and recorded it unchanged as its own `chiref_unconstrained`, so the
# bisection nodes it then added were never allowed to lower it. Applying the
# definition to the merged set is what this line does.
const CHIREF_USED = min(chiref, CHIBEST)
const TREL = CHIREF_USED + 3.84
@printf("\n  %d nodes with models, %d of them admissible: %.4f to %.4f S\n",
        length(Cs), length(adm), minimum(adm), maximum(adm))
@printf("  lowest misfit over the merged nodes: %.4f at C = %.4f S\n", CHIBEST, bestC)

# chiref is what the stores recorded and what the PUBLISHED relative edge was
# measured against, so the threshold drawn on the page stays chiref + 3.84.
# Merging the bisection nodes in can turn up a node below it -- the same
# pathology DECISIONS.md records for the unconstrained solve, one level down --
# and when it does the relative threshold is too high by that much. Said here
# rather than quietly re-based, because re-basing would restate a published
# edge from a script whose job is to draw one.
@printf("\n  REFERENCE MISFIT\n")
@printf("  chiref recorded in the stores : %.4f\n", chiref)
@printf("  lowest over the merged nodes  : %.4f at C = %.4f S\n", CHIBEST, bestC)
@printf("  chiref USED (the lower)       : %.4f  -> relative threshold %.3f\n",
        CHIREF_USED, TREL)
if CHIBEST < chiref - 1e-9
    @printf("  The stores' chiref was NOT the minimum; the lower value is used.\n")
    @printf("  The absolute rule, which drives, is unaffected.\n")
end

# ---- reproduce the published pinned counts, as a check on the predicate
println("\n  cells on a stated limit, recomputed from the stored models:")
for c in Cs
    nb = length(onbound(nodes[c].model, mbounds))
    @printf("    C = %7.4f S  chi2 %9.3f  %2d on a limit  (%s)\n",
            c, nodes[c].chi2, nb, nodes[c].src)
end
let want = Dict(0.0506 => 3, 0.2014 => 2, 0.2260 => 5, 0.2845 => 4, 0.4319 => 13),
    bad = String[]
    for (c, n) in want
        haskey(nodes, c) || continue
        got = length(onbound(nodes[c].model, mbounds))
        got == n || push!(bad, @sprintf("C=%.4f got %d want %d", c, got, n))
    end
    if isempty(bad)
        println("  CHECK: counts reproduce gate_check/results/PINNED_TWO_VARIANT.md")
    else
        @warn "pinned counts do NOT reproduce PINNED_TWO_VARIANT.md: $(join(bad, "; "))"
    end
end

# ---- SNMR DATA FIT.
#
# THE FORWARD IS EXACT, NOT APPROXIMATED. `tools/smr_data.jl` calls
# SNMRForward.forward and stores `Matrix(F.kernel)` (nz x nq) and `F.dz`
# (nz-1). Forwarding a water content w off those stored arrays reproduces the
# stored `V0_clean` to a relative error of 0.000e+00 when, and only when, the
# depth integral is a TRAPEZOID over the kernel nodes:
#
#     V_q = | sum_i  dz_i * ( K[i,q] w[i] + K[i+1,q] w[i+1] ) / 2 |
#
# A lower- or upper-rectangle sum is 2.8% out and a plain diagonal dz weighting
# 1.5%, against a noise level that is 0.3-0.5% of signal -- so either would show
# every model as a 3-5 sigma misfit for a reason that is arithmetic, not data.
# The weights are folded into one nq x nz matrix so the whole ensemble is a
# single multiply.
#
# TWO INPUT SHAPES. The synthetic `smr_data_*.jld2` carries qgrid, <tier>_V,
# <tier>_sigma and the kernel, so the predicted envelope can be drawn. The
# field `smr_transd_*.jld2` carries q, phi and sigma_hat but NO kernel, so the
# measurements and the noise are drawn and the envelope is not.
function trapweights(K, dz)
    nz, nq = size(K); A = zeros(ComplexF64, nq, nz)
    for q in 1:nq
        A[q, 1] += dz[1] * K[1, q] / 2
        for i in 2:nz-1; A[q, i] += (dz[i-1] + dz[i]) * K[i, q] / 2; end
        A[q, nz] += dz[nz-1] * K[nz, q] / 2
    end
    A
end

const FIT = let path = abspath_rel(SMRDATA)
    if !isfile(path)
        @printf("\nSNMR data fit: %s not found; the tab says so.\n", path)
        nothing
    else
        D = load(path)
        if haskey(D, "qgrid")                      # synthetic dataset
            tier = isempty(TIER) ? String(get(SN, "tier", "typical")) : TIER
            vk, sk = tier * "_V", tier * "_σ̂"
            (haskey(D, vk) && haskey(D, sk)) ||
                error("$path has no tier `$tier`; keys are $(sort(collect(keys(D))))")
            q, V, σ = Float64.(D["qgrid"]), Float64.(D[vk]), Float64.(D[sk])
            size(D["kernel"], 1) == NS ||
                error("kernel has $(size(D["kernel"],1)) depths, the ensemble has $NS")
            A = trapweights(D["kernel"], Float64.(D["kernel_dz"]))
            P = abs.(A * θ_smr)                    # nq x ndraws
            χ² = vec(sum(((V .- P) ./ σ) .^ 2, dims = 1))
            @printf("\nSNMR data fit: %s, tier `%s`, %d gates\n", basename(path), tier, length(q))
            @printf("  signal %.0f-%.0f nV, noise %.2f-%.2f nV (%.2f-%.2f%% of signal)\n",
                    1e9minimum(V), 1e9maximum(V), 1e9minimum(σ), 1e9maximum(σ),
                    100minimum(σ ./ V), 100maximum(σ ./ V))
            @printf("  ensemble forwarded through the STORED kernel: chi2/N p05 %.2f p50 %.2f p95 %.2f over N=%d\n",
                    pct(χ²,0.05)/length(q), pct(χ²,0.50)/length(q), pct(χ²,0.95)/length(q), length(q))
            if haskey(D, "V0_clean")
                c0 = sum(((V .- Float64.(D["V0_clean"])) ./ σ) .^ 2) / length(q)
                @printf("  for reference, the generating model's own chi2/N = %.3f\n", c0)
            end
            (q = q, V = V, sig = σ,
             p05 = [pct(view(P,i,:),0.05) for i in 1:length(q)],
             p50 = [pct(view(P,i,:),0.50) for i in 1:length(q)],
             p95 = [pct(view(P,i,:),0.95) for i in 1:length(q)],
             chi2n = pct(χ², 0.50)/length(q), tier = tier, src = basename(path), haspred = true)
        elseif haskey(D, "q") && haskey(D, "phi")  # field sounding
            q, V, σ = Float64.(D["q"]), Float64.(D["phi"]), Float64.(D["sigma_hat"])
            @printf("\nSNMR data fit: %s (field shape), %d gates, no kernel stored\n",
                    basename(path), length(q))
            println("  measurements and noise drawn; the predicted envelope needs a kernel.")
            (q = q, V = V, sig = σ, p05 = Float64[], p50 = Float64[], p95 = Float64[],
             chi2n = NaN, tier = get(D, "sounding_id", ""), src = basename(path), haspred = false)
        else
            @printf("\nSNMR data fit: %s has neither (qgrid,<tier>_V) nor (q,phi); skipped\n",
                    basename(path))
            nothing
        end
    end
end

# ---- the smoothest admissible model per node, when the store exists.
# Only the MIN-CHI2 model is kept by the gate drivers, and it is systematically
# the least smooth: both stores ran lams = [0.01, 0.1] and the winners sit at
# the low end. DECISIONS.md's plan for this figure wanted the smoothest model
# that still passes at each node, which needs new solves. Nothing is solved
# here -- the file is read if it is there and the toggle stays greyed if it is
# not. Accepted shape: the same as the gate stores (`models` plus `lt` or
# `conductance`).
const SMOOTH = if isfile(SMOOTHF)
    S = load(SMOOTHF)
    Cx = haskey(S, "conductance") ? S["conductance"] : (10 .^ S["lt"]) .* S["thick"]
    d = Dict(round(Float64(Cx[i]), digits = 6) => collect(Float64, S["models"][i])
             for i in eachindex(Cx) if solved(collect(Float64, S["models"][i])))
    @printf("\nSmoothest models: %s, %d nodes\n", basename(SMOOTHF), length(d))
    d
else
    @printf("\nNo smoothest-model store at %s; the toggle is greyed.\n", SMOOTHF)
    nothing
end

# ---- the Stage 2 depth limit. Stage 2 does not write one; see the report.
const DEPTHLIM = if !isempty(DLIM)
    v = tryparse(Float64, DLIM)
    if v !== nothing
        @printf("Stage 2 depth limit: %.1f m, given on the command line\n", v)
        (z = v, src = "--depthlimit")
    else
        f = abspath_rel(DLIM); D = load(f)
        k = findfirst(x -> occursin(r"depth_limit|depthlim"i, x), sort(collect(keys(D))))
        k === nothing && error("no depth-limit key in $f; keys are $(sort(collect(keys(D))))")
        kk = sort(collect(keys(D)))[k]
        @printf("Stage 2 depth limit: %.1f m, from %s[%s]\n", D[kk], basename(f), kk)
        (z = Float64(D[kk]), src = basename(f))
    end
else
    local has = haskey(SN, "depth_limit") || haskey(SN, "zlimit")
    if has
        z = Float64(haskey(SN, "depth_limit") ? SN["depth_limit"] : SN["zlimit"])
        @printf("Stage 2 depth limit: %.1f m, from %s\n", z, basename(SNMRF))
        (z = z, src = basename(SNMRF))
    else
        println("\nStage 2 depth limit: NOT STORED.")
        println("  $(basename(SNMRF)) carries no depth-limit key, and stage2_snmr.jl")
        println("  writes none. `depth_limit = \"min_detect_responding\"` is declared at")
        println("  config.toml:490 but the string does not appear in any .jl in this")
        println("  repository, so no stage computes it. Pass --depthlimit <metres|file>")
        println("  to draw one. The SNMR panel says so in plain view.")
        nothing
    end
end

# ---- sensitivity depth, if its store is there.
# Below it the profile is where the regularisation puts it, not where the data
# do, and a cell landing on a limit down there says nothing. Drawn so the
# bottom of the panel is not read as information. frac is a documented
# parameter and the store records which one was used; it is NOT re-chosen here.
const SENS = if isfile(SENSF)
    S = load(SENSF)
    @printf("\nSensitivity depth %.1f m at frac %.3g (%s), %d cells below\n",
            S["zsens"], S["frac"], basename(SENSF), S["ncells_below"])
    (z = Float64(S["zsens"]), frac = Float64(S["frac"]))
else
    @printf("\nNo sensitivity store at %s; the depth marker is omitted.\n", SENSF)
    nothing
end

# ---- petrophysics, READ FROM THE SETTINGS.
# mg_per_sm, the porosity range, the m range behind the lithology keyword and
# the clay ceiling all live in the config. None of them is written into this
# script or into the page: a reader who changes config.toml and rebuilds must
# see the page change, or the page is quoting numbers that are no longer the
# study's.
include(joinpath(ROOT, "src", "config.jl"))
const PET = try
    cfg = SiteConfig.resolve(CFGF)
    pr, ce = cfg["priors"], cfg["ceilings"]
    ml, mh = haskey(pr, "m_range") ? Tuple(Float64.(pr["m_range"])) :
             SiteConfig.LITHOLOGY_M[pr["m_lithology"]]
    @printf("\nPetrophysics from %s: mg_per_sm %.0f, theta %.2f-%.2f, m %.2f-%.2f\n",
            basename(CFGF), pr["mg_per_sm"], pr["porosity_range"]..., ml, mh)
    @printf("  clay ceiling %.1f x %.1f x (1 - theta/%.2f) S/m, clay prior %.0e-%.1f S/m\n",
            ce["sigma_clay_max"], ce["clay_c0"], ce["clay_phi_ref"], pr["clay"]["product"]...)
    (mgl = Float64(pr["mg_per_sm"]), phi = Float64.(pr["porosity_range"]), m = [ml, mh],
     scmax = Float64(ce["sigma_clay_max"]), c0 = Float64(ce["clay_c0"]),
     phiref = Float64(ce["clay_phi_ref"]), cp = Float64.(pr["clay"]["product"]),
     src = basename(CFGF))
catch e
    error("could not read petrophysics from $CFGF: $e\n" *
          "Pass --config <site toml>. Nothing is defaulted here on purpose.")
end

# ---- the TRUE clay term in the window, computed the robust way.
# NOT `sigma_clay .* clay`: in these truth files `sigma_clay` IS the clay term
# already, and multiplying it by the clay fraction again is the error that put
# conditions_table out by 5-16x (DECISIONS.md). Taken as the residual instead,
# sigma_bulk - sigma_w * theta_total^m, which is the relation Stage 5 uses.
const CLAYTRUE = let z = collect(Float64, TR["z"]),
                     k = findall(q -> WTOP <= q <= WBOT, collect(Float64, TR["z"]))
    m   = Float64(TR["m_archie"])
    σb  = mean(collect(Float64, TR["σ_bulk"])[k])
    wat = mean((collect(Float64, TR["σw"]) .* collect(Float64, TR["θ_total"]) .^ m)[k])
    frac = (σb - wat) / σb
    tds  = mean(collect(Float64, TR["σw"])[k]) * PET.mgl
    @printf("  true clay term in the window %.5f S/m, %.1f%% of sigma_bulk %.5f\n",
            σb - wat, 100frac, σb)
    @printf("  true pore water %.0f mg/L\n", tds)
    (term = σb - wat, frac = frac, sigb = σb, tds = tds)
end

# ============================================================ 4. the payload
rowmax(h) = [max(1, maximum(r)) for r in h]

"""Collapse `aem3P_b1e1 … aem3P_b1e4` to `aem3P_b1e1-e4`; leave anything else
alone. The provenance line has to name the folders and also has to stay one
line, and four near-identical names is mostly repetition."""
function collapse_names(ns)
    length(ns) < 2 && return join(ns, ", ")
    k = 0
    while k < minimum(length.(ns)) && all(n -> n[k+1] == ns[1][k+1], ns)
        k += 1
    end
    k < 4 ? join(ns, ", ") : ns[1] * "-" * ns[end][k+1:end]
end

gate_models = [nodes[c].model for c in Cs]
gate_chi2   = [nodes[c].chi2 for c in Cs]
gate_ob     = [onbound(nodes[c].model, mbounds) for c in Cs]
gate_adm    = [nodes[c].chi2 <= TABS for c in Cs]

function js_payload()
    io = IOBuffer()
    print(io, "const D={")
    print(io, "gen:", jstr(Dates.format(now(), "yyyy-mm-dd HH:MM")), ",")
    print(io, "win:[", jnum(WTOP), ",", jnum(WBOT), "],zmax:", jnum(ZMAX), ",")
    print(io, "thick:", jnum(thick), ",ndata:", ndata, ",")
    print(io, "tabs:", jnum(TABS), ",trel:", jnum(TREL), ",chiref:", jnum(CHIREF_USED),
          ",chibest:", jnum(CHIBEST), ",chirefstore:", jnum(chiref), ",")
    print(io, "mb:", jarr(mbounds), ",")
    # AEM
    print(io, "aem:{z:", jarr(zall), ",top:", jarr(celltop), ",bot:", jarr(cellbot),
          ",lo:", jnum(AEM_LO), ",hi:", jnum(AEM_HI), ",nb:", AEM_NB,
          ",h:[", join([jiarr(r) for r in aem_hist], ","), "]",
          ",rmax:", jiarr(rowmax(aem_hist)),
          ",p05:", jarr(aem_p05), ",p50:", jarr(aem_p50), ",p95:", jarr(aem_p95),
          ",truth:", jarr(σtruth_aem),
          ",ndraw:", size(σ_aem, 2), ",nens:", length(AEMDIRS),
          ",src:", jstr(collapse_names(basename.(AEMDIRS))), "},")
    # SNMR
    print(io, "snmr:{z:", jarr(zsmr),
          ",lo:", jnum(SNMR_LO), ",hi:", jnum(SNMR_HI), ",nb:", SNMR_NB,
          ",h:[", join([jiarr(r) for r in snmr_hist], ","), "]",
          ",rmax:", jiarr(rowmax(snmr_hist)),
          ",p05:", jarr(snmr_p05), ",p50:", jarr(snmr_p50), ",p95:", jarr(snmr_p95),
          ",truth:", jarr(θtrue_smr),
          ",ndraw:", size(θ_smr, 2), ",nens:", length(snmr_ens),
          ",src:", jstr(basename(SNMRF)), ",note:", jstr(SNMRNOTE),
          ",dlim:", DEPTHLIM === nothing ? "null" : jnum(DEPTHLIM.z), "},")
    # gate
    print(io, "gate:{C:", jarr(Cs), ",chi2:", jarr(gate_chi2),
          ",adm:[", join([x ? "1" : "0" for x in gate_adm], ","), "]",
          ",best:", jnum(bestC),
          ",m:[", join([jarr(m) for m in gate_models], ","), "]",
          ",ob:[", join([jiarr(v .- 1) for v in gate_ob], ","), "]",
          ",ms:", SMOOTH === nothing ? "null" :
              (all(haskey(SMOOTH, c) for c in adm) ?
                  string("[", join([haskey(SMOOTH, c) ? jarr(SMOOTH[c]) : "null" for c in Cs], ","), "]") :
                  "null"),
          ",msn:", SMOOTH === nothing ? 0 : count(c -> haskey(SMOOTH, c), adm),
          ",msneed:", length(adm),
          ",src:", jstr(join(basename.(GATEF), " + ")), "},")
    print(io, "truthC:", jnum(Ctruth_grid), ",")
    print(io, "sens:", SENS === nothing ? "null" :
        string("{z:", jnum(SENS.z), ",frac:", jnum(SENS.frac), "}"), ",")
    print(io, "pet:{mgl:", jnum(PET.mgl), ",phi:", jarr(PET.phi), ",m:", jarr(PET.m),
          ",scmax:", jnum(PET.scmax), ",c0:", jnum(PET.c0), ",phiref:", jnum(PET.phiref),
          ",cp:", jarr(PET.cp), ",src:", jstr(PET.src), "},")
    print(io, "claytrue:{term:", jnum(CLAYTRUE.term), ",frac:", jnum(CLAYTRUE.frac),
          ",tds:", jnum(CLAYTRUE.tds), "},")
    print(io, "fit:", FIT === nothing ? "null" : string(
        "{q:", jarr(FIT.q), ",V:", jarr(FIT.V), ",sig:", jarr(FIT.sig),
        ",p05:", jarr(FIT.p05), ",p50:", jarr(FIT.p50), ",p95:", jarr(FIT.p95),
        ",chi2n:", jnum(FIT.chi2n), ",tier:", jstr(FIT.tier),
        ",src:", jstr(FIT.src), ",haspred:", FIT.haspred ? "true" : "false", "}"))
    print(io, "};")
    String(take!(io))
end

"""
ONE DATA SOURCE. The payload the page embeds is also written to
report_<column>.json, and the page embeds THAT STRING -- not a second
serialisation of the same numbers. The figures read the JSON, so a figure and
the page cannot disagree about what the data are: there is one serialiser and
two consumers.

js_payload emits a JavaScript object literal with bare keys, which is not JSON.
The keys are quoted here by a transformation on the generated string, which is
safe because every string VALUE in the payload is already JSON-quoted and none
contains the `,identifier:` or `{identifier:` shape the pattern matches. That is
not left to argument: the result is parsed back before anything is written, and
the run stops if it is not valid JSON.
"""
function json_payload()
    js = js_payload()
    body = js[findfirst("{", js)[1]:end]          # drop the `const D=` prefix
    body = rstrip(rstrip(body), ';')              # and the statement's semicolon
    j = replace(body, r"([{,])([A-Za-z_][A-Za-z0-9_]*):" => s"\1\"\2\":")
    j
end

# ============================================================== 5. the page
const PAGE = raw"""
<title>Test One, Depth by Depth</title>
<style>
:root{
  --bg:#f1f4f4;--panel:#fff;--ink:#16262b;--soft:#4a5f65;--mute:#7b9095;
  --rule:#ccd8d8;--line:#9fb3b3;--acc:#0d6a68;--acc2:#1d7fa8;--hot:#b4502f;
  --warn:#8a5a0c;--warnbg:#f7edd8;--grid:#e3eaea;color-scheme:light;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --bg:#0d1517;--panel:#162023;--ink:#e3ebea;--soft:#aec0c0;--mute:#8098a0;
  --rule:#2a393d;--line:#4a6066;--acc:#51b5ad;--acc2:#58a9d0;--hot:#dd7d60;
  --warn:#dcae57;--warnbg:#2a2314;--grid:#223034;color-scheme:dark;}}
:root[data-theme="dark"]{
  --bg:#0d1517;--panel:#162023;--ink:#e3ebea;--soft:#aec0c0;--mute:#8098a0;
  --rule:#2a393d;--line:#4a6066;--acc:#51b5ad;--acc2:#58a9d0;--hot:#dd7d60;
  --warn:#dcae57;--warnbg:#2a2314;--grid:#223034;color-scheme:dark;}
[hidden]{display:none!important}
img{max-width:100%}
html,body{height:100%}
body{background:var(--bg);color:var(--ink);margin:0;
  font:13px/1.45 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  -webkit-font-smoothing:antialiased}
.app{display:flex;flex-direction:column;min-height:100%;padding:12px 16px 18px;box-sizing:border-box;gap:10px}
h1{font-size:15px;font-weight:650;letter-spacing:-.01em;margin:0}
.head{display:flex;flex-wrap:wrap;align-items:baseline;gap:6px 16px}
.prov{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:10.5px;color:var(--mute);
  line-height:1.5;word-break:break-word;flex:1 1 320px;min-width:0}
.flag{color:var(--warn);background:var(--warnbg);padding:1px 5px;border-radius:3px;font-weight:600}
.panels{display:grid;grid-template-columns:repeat(3,1fr);gap:10px;flex:1 1 auto;min-height:330px}
.pan{background:var(--panel);border:1px solid var(--rule);border-radius:4px;
  display:flex;flex-direction:column;min-width:0;overflow:hidden}
.pt{font-size:11.5px;font-weight:650;padding:7px 10px 0;color:var(--ink);
  display:flex;justify-content:space-between;align-items:baseline;gap:8px}
.pt em{font-style:normal;font-weight:400;font-size:10px;color:var(--mute);
  font-family:ui-monospace,Menlo,Consolas,monospace}
.tabs{display:inline-flex;gap:2px}
.tab{font:inherit;font-size:11.5px;font-weight:600;color:var(--mute);background:none;
  border:none;border-bottom:2px solid transparent;padding:0 2px 2px;cursor:pointer}
.tab.on{color:var(--ink);border-bottom-color:var(--acc)}
.tab:focus-visible{outline:2px solid var(--acc);outline-offset:2px}
.cw{flex:1 1 auto;position:relative;min-height:0}
canvas{position:absolute;inset:0;width:100%;height:100%;display:block}
.strip{display:grid;grid-template-columns:minmax(260px,1.05fr) minmax(300px,1.35fr);
  gap:10px;align-items:stretch}
.box{background:var(--panel);border:1px solid var(--rule);border-radius:4px;padding:9px 11px;
  display:flex;flex-direction:column;gap:7px;min-width:0}
.chiw{position:relative;height:124px;min-height:124px}
.ctl{display:flex;align-items:center;gap:9px;font-size:11px;color:var(--soft)}
.ctl label{white-space:nowrap;font-size:10.5px;color:var(--mute);min-width:88px}
.ctl output{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:11px;color:var(--ink);
  min-width:66px;text-align:right;font-variant-numeric:tabular-nums}
input[type=range]{flex:1 1 auto;min-width:60px;accent-color:var(--acc);height:18px;margin:0}
/* FIXED layout, explicit column widths and tabular figures: the numbers
   change on every slider step and nothing may move when they do. */
#say{width:100%;table-layout:fixed;border-collapse:collapse;font-size:12px}
#say col.k{width:34%} #say col.v{width:36%} #say col.t{width:30%}
#say th{font-size:10px;font-weight:600;letter-spacing:.06em;text-transform:uppercase;
  color:var(--mute);text-align:right;padding:0 0 3px;border-bottom:1px solid var(--rule)}
#say th.k{text-align:left}
#say td{padding:3px 0;border-bottom:1px solid var(--rule);vertical-align:baseline}
#say td.k{color:var(--mute);font-size:11px}
#say td.v,#say td.t{text-align:right;font-variant-numeric:tabular-nums;
  font-family:ui-monospace,Menlo,Consolas,monospace;white-space:nowrap;overflow:hidden;
  text-overflow:ellipsis}
#say td.v{color:var(--ink);font-weight:600}
#say td.t{color:var(--soft)}
#say td.v em{font-style:normal;font-weight:400;color:var(--mute);font-size:10.5px}
#say tr:last-child td{border-bottom:none}
.btn{font:inherit;font-size:11px;padding:4px 10px;border-radius:3px;cursor:pointer;
  border:1px solid var(--line);background:var(--panel);color:var(--ink)}
.btn:hover{border-color:var(--acc);color:var(--acc)}
.btn:focus-visible{outline:2px solid var(--acc);outline-offset:1px}
.btn[disabled]{opacity:.45;cursor:not-allowed;border-style:dashed}
.btn[disabled]:hover{border-color:var(--line);color:var(--ink)}
.seg{min-width:120px;text-align:left}
.hint{font-size:10px;color:var(--mute);font-family:ui-monospace,Menlo,Consolas,monospace}
.row{display:flex;gap:8px;align-items:center;flex-wrap:wrap}
textarea{width:100%;box-sizing:border-box;height:72px;font-family:ui-monospace,Menlo,Consolas,monospace;
  font-size:10px;background:var(--bg);color:var(--ink);border:1px solid var(--rule);
  border-radius:3px;padding:5px;resize:vertical}
.key{display:flex;flex-wrap:wrap;gap:3px 12px;font-size:10px;color:var(--mute);
  font-family:ui-monospace,Menlo,Consolas,monospace}
.key i{display:inline-block;width:13px;height:2.5px;vertical-align:middle;margin-right:4px;border-radius:2px}
@media (max-width:860px){
  .panels{grid-template-columns:1fr;min-height:0}
  .cw{height:260px}
  .strip{grid-template-columns:1fr}
}
</style>

<div class="app">
  <div class="head">
    <h1>Test one &mdash; conductivity, water content, gate models</h1>
    <div class="prov" id="prov"></div>
  </div>

  <div class="panels">
    <div class="pan">
      <div class="pt"><span>Conductivity, AEM ensemble</span><em id="n1"></em></div>
      <div class="cw"><canvas id="c1"></canvas></div>
    </div>
    <div class="pan">
      <div class="pt">
        <span class="tabs"><button class="tab on" id="t2a">Water content</button><button class="tab" id="t2b">Data fit</button></span>
        <em id="n2"></em>
      </div>
      <div class="cw"><canvas id="c2"></canvas></div>
    </div>
    <div class="pan">
      <div class="pt"><span>Conductivity, gate models</span><em id="n3"></em></div>
      <div class="cw"><canvas id="c3"></canvas></div>
    </div>
  </div>

  <div class="strip">
    <div class="box">
      <div class="chiw"><canvas id="c4"></canvas></div>
      <div class="ctl">
        <label for="sC">window conductance</label>
        <input type="range" id="sC" min="0" max="1" step="1" value="0">
        <output id="oC"></output>
      </div>
      <div class="ctl">
        <label for="bM">model shown</label>
        <button class="btn seg" id="bM" aria-pressed="false"></button>
        <span class="hint" id="hM"></span>
      </div>
    </div>
    <div class="box">
      <table id="say"><tbody></tbody></table>
      <div class="ctl">
        <label for="sP">porosity &theta;</label>
        <input type="range" id="sP" min="0.25" max="0.40" step="0.005" value="0.325">
        <output id="oP"></output>
      </div>
      <div class="ctl">
        <label for="sM">exponent m</label>
        <input type="range" id="sM" min="1.5" max="1.9" step="0.01" value="1.7">
        <output id="oM"></output>
      </div>
      <div class="ctl">
        <label for="sK">clay carries</label>
        <input type="range" id="sK" min="0" max="95" step="1" value="0">
        <output id="oK"></output>
      </div>
      <div class="row">
        <button class="btn" id="bC">Export CSV</button>
        <button class="btn" id="bS" hidden>Select all</button>
        <span class="key" id="key"></span>
      </div>
      <textarea id="ta" hidden readonly aria-label="CSV"></textarea>
    </div>
  </div>
</div>
"""

const SCRIPT = raw"""
<script>
(function(){
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const F = (x,n) => x.toFixed(n===undefined?3:n);
const clamp=(a,b,c)=>Math.max(b,Math.min(c,a));

// ---- depth grids -> cell edges
const aemTop = D.aem.top.slice();
const aemBot = D.aem.bot.slice();
const smrEdge = (()=>{const z=D.snmr.z, e=[Math.max(0,z[0]-(z[1]-z[0])/2)];
  for(let i=0;i<z.length-1;i++) e.push((z[i]+z[i+1])/2);
  e.push(z[z.length-1]+(z[z.length-1]-z[z.length-2])/2); return e;})();

const PAD={l:42,r:10,t:10,b:30};

// ---- inferno, the matplotlib ramp, as 9 stops linearly interpolated.
// Density is drawn as an image on its own dark ground rather than as ink on
// the page, so the panel reads the same in either theme and the overlay lines
// have one background to be legible against instead of two.
const INFERNO=[[0,0,4],[31,12,72],[85,15,109],[136,34,106],[186,54,85],
               [227,89,51],[249,142,9],[249,212,45],[252,255,164]];
function ramp(t){
  t=Math.max(0,Math.min(1,t));
  const u=t*(INFERNO.length-1), i=Math.min(INFERNO.length-2,Math.floor(u)), f=u-i;
  const a=INFERNO[i], b=INFERNO[i+1];
  return 'rgb('+Math.round(a[0]+(b[0]-a[0])*f)+','+Math.round(a[1]+(b[1]-a[1])*f)+','
               +Math.round(a[2]+(b[2]-a[2])*f)+')';
}
const GROUND=ramp(0);
// Inferno runs black -> purple -> red -> orange -> near-white yellow, so the
// overlays are cold hues the ramp never reaches, and the truth carries a dark
// halo: the p50 sits exactly where density peaks, which is the palest part of
// the ramp, and a plain white line would disappear into it.
//   white  the true model, in every panel
//   cyan   what the data produced -- percentiles, and the selected gate model
//   amber  the lowest-misfit model, a reference and nothing else
//   green  a cell resting on a stated limit
const C_P50='#00e5ff', C_PTL='rgba(0,229,255,.72)', C_TRUE='#ffffff',
      C_FAINT='rgba(255,255,255,.22)', C_BEST='#ffb703',
      C_SEL='#00e5ff', C_LIM='#ffffff', C_GRID='rgba(255,255,255,.13)',
      C_HALO='rgba(0,0,0,.6)';

function setup(cv){
  const r=cv.getBoundingClientRect(), d=window.devicePixelRatio||1;
  cv.width=Math.max(1,Math.round(r.width*d)); cv.height=Math.max(1,Math.round(r.height*d));
  const g=cv.getContext('2d'); g.setTransform(d,0,0,d,0,0);
  g.clearRect(0,0,r.width,r.height);
  return {g,w:r.width,h:r.height};
}
// axis furniture sits on the PAGE and is themed; the plot interior is the
// plasma ground and everything drawn inside it is light
function frame(g,w,h,xlab,xticks,xof,ylab){
  const x0=PAD.l,x1=w-PAD.r,y0=PAD.t,y1=h-PAD.b;
  g.fillStyle=GROUND; g.fillRect(x0,y0,x1-x0,y1-y0);
  g.font='9.5px ui-monospace,monospace';
  g.strokeStyle=C_GRID; g.lineWidth=1;
  for(const t of xticks){const X=xof(t); if(X<x0+1||X>x1-1) continue;
    g.beginPath(); g.moveTo(X,y0); g.lineTo(X,y1); g.stroke();}
  const step = D.zmax<=30?5:(D.zmax<=60?10:20);
  for(let d=step; d<D.zmax; d+=step){
    const Y=y0+(d/D.zmax)*(y1-y0);
    g.beginPath(); g.moveTo(x0,Y); g.lineTo(x1,Y); g.stroke();}
  g.fillStyle=css('--mute'); g.textAlign='center'; g.textBaseline='top';
  for(const t of xticks){const X=xof(t); if(X<x0-1||X>x1+1) continue;
    g.fillText(String(t.lab!==undefined?t.lab:t.v),X,y1+4);}
  g.textAlign='right'; g.textBaseline='middle';
  for(let d=0; d<=D.zmax; d+=step){
    const Y=y0+(d/D.zmax)*(y1-y0); g.fillText(String(d),x0-5,Y);}
  g.save(); g.translate(11,(y0+y1)/2); g.rotate(-Math.PI/2);
  g.textAlign='center'; g.textBaseline='middle'; g.fillText(ylab,0,0); g.restore();
  g.strokeStyle=css('--line'); g.lineWidth=1; g.strokeRect(x0,y0,x1-x0,y1-y0);
  g.fillStyle=css('--mute'); g.textAlign='right'; g.textBaseline='top';
  g.fillText(xlab,x1,y1+15);
  return {x0,x1,y0,y1};
}
const yOf=(A,d)=>A.y0+(clamp(d,0,D.zmax)/D.zmax)*(A.y1-A.y0);
function band(g,A){
  const a=yOf(A,D.win[0]), b=yOf(A,D.win[1]);
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  g.fillStyle='rgba(255,255,255,.07)'; g.fillRect(A.x0,a,A.x1-A.x0,b-a);
  g.strokeStyle='rgba(255,255,255,.55)'; g.lineWidth=1; g.setLineDash([4,3]);
  g.beginPath(); g.moveTo(A.x0,a); g.lineTo(A.x1,a); g.moveTo(A.x0,b); g.lineTo(A.x1,b);
  g.stroke(); g.setLineDash([]); g.restore();
}
// Below the sensitivity depth the models are drawn, then dimmed: they are
// still what the solver returned, but they are not a statement about the
// ground, so they should not compete with the part that is.
function veilAt(g,A,z,label){
  if(z===null||z===undefined||z>=D.zmax) return;
  const Y=yOf(A,z);
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  g.fillStyle='rgba(0,0,0,.52)'; g.fillRect(A.x0,Y,A.x1-A.x0,A.y1-Y);
  g.strokeStyle='rgba(255,255,255,.65)'; g.lineWidth=1.2; g.setLineDash([2,3]);
  g.beginPath(); g.moveTo(A.x0,Y); g.lineTo(A.x1,Y); g.stroke(); g.setLineDash([]);
  g.fillStyle='rgba(255,255,255,.8)'; g.font='9px ui-monospace,monospace';
  g.textAlign='left'; g.textBaseline='bottom';
  g.fillText(label+' '+z.toFixed(0)+' m', A.x0+4, Y-2);
  g.restore();
}
function note(g,A,txt){
  g.save(); g.fillStyle='rgba(255,255,255,.62)'; g.font='9px ui-monospace,monospace';
  g.textAlign='left'; g.textBaseline='bottom'; g.fillText(txt,A.x0+4,A.y1-4); g.restore();
}
function veil(g,A){
  if(!D.sens || D.sens.z>=D.zmax) return;
  const Y=yOf(A,D.sens.z);
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  g.fillStyle='rgba(0,0,0,.52)'; g.fillRect(A.x0,Y,A.x1-A.x0,A.y1-Y);
  g.restore();
  sensline(g,A);
}
function sensline(g,A){
  if(!D.sens || D.sens.z>=D.zmax) return;
  const Y=yOf(A,D.sens.z);
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  g.strokeStyle='rgba(255,255,255,.65)'; g.lineWidth=1.2; g.setLineDash([2,3]);
  g.beginPath(); g.moveTo(A.x0,Y); g.lineTo(A.x1,Y); g.stroke(); g.setLineDash([]);
  g.fillStyle='rgba(255,255,255,.8)'; g.font='9px ui-monospace,monospace';
  g.textAlign='left'; g.textBaseline='bottom';
  g.fillText('sensitivity depth '+D.sens.z.toFixed(0)+' m', A.x0+4, Y-2);
  g.restore();
}
// piecewise-constant: for a layered model, which is what it is
function steps(g,A,xof,v,top,bot,col,wd,dash,halo,alpha){
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  if(alpha!==undefined) g.globalAlpha=alpha;
  const path=()=>{g.beginPath(); let started=false;
    for(let i=0;i<v.length;i++){
      if(top[i]>D.zmax) break;
      if(!isFinite(v[i])) continue;
      const X=xof(v[i]), ya=yOf(A,top[i]), yb=yOf(A,Math.min(bot[i],D.zmax));
      if(!started){g.moveTo(X,ya);started=true;} else {g.lineTo(X,ya);}
      g.lineTo(X,yb);}};
  if(halo){g.strokeStyle=C_HALO; g.lineWidth=wd+2.2; g.setLineDash([]); path(); g.stroke();}
  g.strokeStyle=col; g.lineWidth=wd; g.setLineDash(dash||[]); path(); g.stroke();
  g.setLineDash([]); g.restore();
}
// smooth: a Catmull-Rom through the per-cell values at their cell centres.
// The percentiles are a summary of a distribution that varies continuously
// with depth, not a layered model, so they are drawn as a curve. The truth
// and the gate's models stay as steps, because those ARE layered models.
function smooth(g,A,xof,v,top,bot,col,wd,dash){
  const pts=[];
  for(let i=0;i<v.length;i++){
    const zc=(top[i]+bot[i])/2;
    if(!isFinite(v[i])) continue;
    pts.push({x:xof(v[i]),y:yOf(A,zc)});
    if(zc>D.zmax) break;
  }
  if(pts.length<2) return;
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  g.strokeStyle=col; g.lineWidth=wd; g.setLineDash(dash||[]); g.lineJoin='round';
  g.beginPath(); g.moveTo(pts[0].x,pts[0].y);
  for(let i=0;i<pts.length-1;i++){
    const p0=pts[i>0?i-1:0], p1=pts[i], p2=pts[i+1],
          p3=pts[i+2<pts.length?i+2:pts.length-1];
    g.bezierCurveTo(p1.x+(p2.x-p0.x)/6, p1.y+(p2.y-p0.y)/6,
                    p2.x-(p3.x-p1.x)/6, p2.y-(p3.y-p1.y)/6, p2.x, p2.y);
  }
  g.stroke(); g.setLineDash([]); g.restore();
}
// per-depth normalisation: each row is the marginal at that depth scaled to
// its own maximum, so a well-determined depth and a poorly determined one are
// both readable. A single global scale would wash out everything below the
// depth of investigation.
function hist(g,A,H,rmax,lo,hi,nb,top,bot){
  const bw=(hi-lo)/nb;
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  for(let i=0;i<H.length;i++){
    if(top[i]>D.zmax) break;
    const ya=yOf(A,top[i]), yb=yOf(A,Math.min(bot[i],D.zmax)), hgt=Math.max(1,yb-ya);
    const row=H[i], mx=rmax[i];
    for(let j=0;j<nb;j++){
      const c=row[j]; if(!c) continue;
      const xa=A.x0+((j*bw))/(hi-lo)*(A.x1-A.x0);
      const xb=A.x0+(((j+1)*bw))/(hi-lo)*(A.x1-A.x0);
      g.fillStyle=ramp(Math.sqrt(c/mx));
      g.fillRect(xa,ya,Math.max(1,xb-xa)+0.4,hgt+0.4);
    }
  }
  g.restore();
}

const LOGT=[{v:-3,lab:'.001'},{v:-2,lab:'.01'},{v:-1,lab:'0.1'},{v:0,lab:'1'}];
const sxLog = A => (s => A.x0+((Math.log10(s)-D.aem.lo)/(D.aem.hi-D.aem.lo))*(A.x1-A.x0));

function drawAEM(){
  const {g,w,h}=setup(document.getElementById('c1')); if(w<8) return;
  const xof=t=>PAD.l+((t.v-D.aem.lo)/(D.aem.hi-D.aem.lo))*(w-PAD.l-PAD.r);
  const A=frame(g,w,h,'conductivity S/m',LOGT,xof,'depth m');
  const X=sxLog(A);
  hist(g,A,D.aem.h,D.aem.rmax,D.aem.lo,D.aem.hi,D.aem.nb,aemTop,aemBot);
  band(g,A);
  smooth(g,A,X,D.aem.p05,aemTop,aemBot,C_PTL,1.3,[5,3]);
  smooth(g,A,X,D.aem.p95,aemTop,aemBot,C_PTL,1.3,[5,3]);
  smooth(g,A,X,D.aem.p50,aemTop,aemBot,C_P50,2.1);
  steps(g,A,X,D.aem.truth,aemTop,aemBot,C_TRUE,2,null,true);
  sensline(g,A);
}
function drawSNMR(){
  const {g,w,h}=setup(document.getElementById('c2')); if(w<8) return;
  const T=[{v:0,lab:'0'},{v:0.1,lab:'0.1'},{v:0.2,lab:'0.2'},{v:0.3,lab:'0.3'},{v:0.4,lab:'0.4'}];
  const xof=t=>PAD.l+((t.v-D.snmr.lo)/(D.snmr.hi-D.snmr.lo))*(w-PAD.l-PAD.r);
  const A=frame(g,w,h,'water content v/v',T,xof,'depth m');
  const X=s=>A.x0+((s-D.snmr.lo)/(D.snmr.hi-D.snmr.lo))*(A.x1-A.x0);
  const top=smrEdge.slice(0,-1), bot=smrEdge.slice(1);
  hist(g,A,D.snmr.h,D.snmr.rmax,D.snmr.lo,D.snmr.hi,D.snmr.nb,top,bot);
  band(g,A);
  smooth(g,A,X,D.snmr.p05,top,bot,C_PTL,1.3,[5,3]);
  smooth(g,A,X,D.snmr.p95,top,bot,C_PTL,1.3,[5,3]);
  smooth(g,A,X,D.snmr.p50,top,bot,C_P50,2.1);
  steps(g,A,X,D.snmr.truth,top,bot,C_TRUE,2,null,true);
  if(D.snmr.dlim===null||D.snmr.dlim===undefined) note(g,A,'depth limit not stored');
  else veilAt(g,A,D.snmr.dlim,'depth limit');
}
function drawFit(){
  const {g,w,h}=setup(document.getElementById('c2')); if(w<8) return;
  const F=D.fit;
  if(!F){ const A=frame(g,w,h,'',[],()=>-99,'');
    note(g,A,'SNMR data file not found'); return; }
  const qs=F.q, all=F.V.concat(F.sig).concat(F.haspred?F.p05.concat(F.p95):[]).filter(v=>v>0);
  const vlo=Math.log10(Math.min(...all))-0.08, vhi=Math.log10(Math.max(...all))+0.08;
  const qlo=Math.log10(Math.min(...qs))-0.05, qhi=Math.log10(Math.max(...qs))+0.05;
  const T=[]; for(let d=Math.ceil(qlo);d<=Math.floor(qhi);d++)
    T.push({v:d,lab:(d<0?'1e'+d:(d<4?String(Math.pow(10,d)):'1e'+d))});
  const xof=t=>PAD.l+((t.v-qlo)/(qhi-qlo))*(w-PAD.l-PAD.r);
  // Amplitude spans three decades while the noise is 0.07-2.5% of signal, so
  // a single log axis shows the sounding curve and hides the fit completely.
  // The panel is therefore split: the curve above, and below it the residuals
  // in units of the noise, which is where a 5-sigma misfit is actually visible.
  const x0=PAD.l,x1=w-PAD.r,y0=PAD.t,yB=h-PAD.b, SPLIT=0.60;
  const y1=y0+(yB-y0)*SPLIT-8, y2=y0+(yB-y0)*SPLIT+8;
  const X=q=>x0+((Math.log10(q)-qlo)/(qhi-qlo))*(x1-x0);
  const Y=v=>y1-((Math.log10(Math.max(v,1e-30))-vlo)/(vhi-vlo))*(y1-y0);
  g.fillStyle=GROUND; g.fillRect(x0,y0,x1-x0,y1-y0); g.fillRect(x0,y2,x1-x0,yB-y2);
  g.font='9.5px ui-monospace,monospace';
  g.strokeStyle=C_GRID; g.lineWidth=1;
  for(const t of T){const px=xof(t); if(px<x0+1||px>x1-1)continue;
    g.beginPath();g.moveTo(px,y0);g.lineTo(px,y1);g.moveTo(px,y2);g.lineTo(px,yB);g.stroke();}
  for(let d=Math.ceil(vlo);d<=Math.floor(vhi);d++){
    const py=Y(Math.pow(10,d)); if(py<y0||py>y1)continue;
    g.beginPath();g.moveTo(x0,py);g.lineTo(x1,py);g.stroke();}
  g.fillStyle=css('--mute'); g.textAlign='center'; g.textBaseline='top';
  for(const t of T){const px=xof(t); if(px<x0-1||px>x1+1)continue; g.fillText(t.lab,px,yB+4);}
  g.textAlign='right'; g.textBaseline='middle';
  for(let d=Math.ceil(vlo);d<=Math.floor(vhi);d++){
    const py=Y(Math.pow(10,d)); if(py<y0||py>y1)continue;
    g.fillText(d<=0?'1e'+d:(d<4?String(Math.pow(10,d)):'1e'+d),x0-5,py);}
  g.save(); g.translate(11,(y0+y1)/2); g.rotate(-Math.PI/2);
  g.textAlign='center'; g.textBaseline='middle'; g.fillText('amplitude V',0,0); g.restore();

  // ---- upper: the sounding curve
  g.save(); g.beginPath(); g.rect(x0,y0,x1-x0,y1-y0); g.clip();
  if(F.haspred){
    g.strokeStyle=C_P50; g.lineWidth=1.8; g.beginPath();
    qs.forEach((q,i)=>{const px=X(q),py=Y(F.p50[i]); i?g.lineTo(px,py):g.moveTo(px,py);});
    g.stroke();
  }
  g.strokeStyle=C_BEST; g.lineWidth=1.4; g.setLineDash([5,3]); g.beginPath();
  qs.forEach((q,i)=>{const px=X(q),py=Y(F.sig[i]); i?g.lineTo(px,py):g.moveTo(px,py);});
  g.stroke(); g.setLineDash([]);
  g.fillStyle=C_TRUE;
  qs.forEach((q,i)=>{g.beginPath(); g.arc(X(q),Y(F.V[i]),2.6,0,6.284); g.fill();});
  g.restore();
  g.strokeStyle=css('--line'); g.lineWidth=1; g.strokeRect(x0,y0,x1-x0,y1-y0);
  g.fillStyle=css('--mute'); g.font='9px ui-monospace,monospace';
  g.textAlign='left'; g.textBaseline='top'; g.fillText('noise',x0+4,Y(F.sig[0])+3);

  // ---- lower: residuals in units of the noise
  const res = F.haspred ? qs.map((q,i)=>(F.V[i]-F.p50[i])/F.sig[i]) : qs.map(()=>0);
  const rhi = F.haspred
    ? Math.max(3, Math.ceil(Math.max(...res.map(Math.abs),
        ...qs.map((q,i)=>Math.abs((F.V[i]-F.p05[i])/F.sig[i])),
        ...qs.map((q,i)=>Math.abs((F.V[i]-F.p95[i])/F.sig[i])))))
    : 3;
  const R=r=>y2+(yB-y2)/2-(r/rhi)*((yB-y2)/2-3);
  g.save(); g.beginPath(); g.rect(x0,y2,x1-x0,yB-y2); g.clip();
  g.fillStyle='rgba(255,255,255,.13)'; g.fillRect(x0,R(1),x1-x0,R(-1)-R(1));
  if(F.haspred){
    g.fillStyle='rgba(0,229,255,.22)'; g.beginPath();
    qs.forEach((q,i)=>{const px=X(q),py=R((F.V[i]-F.p95[i])/F.sig[i]); i?g.lineTo(px,py):g.moveTo(px,py);});
    for(let i=qs.length-1;i>=0;i--) g.lineTo(X(qs[i]),R((F.V[i]-F.p05[i])/F.sig[i]));
    g.closePath(); g.fill();
  }
  g.strokeStyle='rgba(255,255,255,.55)'; g.lineWidth=1;
  g.beginPath(); g.moveTo(x0,R(0)); g.lineTo(x1,R(0)); g.stroke();
  if(F.haspred){
    g.strokeStyle=C_TRUE; g.lineWidth=1.2; g.fillStyle=C_TRUE;
    qs.forEach((q,i)=>{const px=X(q);
      g.beginPath(); g.moveTo(px,R(0)); g.lineTo(px,R(res[i])); g.stroke();
      g.beginPath(); g.arc(px,R(res[i]),2.4,0,6.284); g.fill();});
  }
  g.restore();
  g.strokeStyle=css('--line'); g.lineWidth=1; g.strokeRect(x0,y2,x1-x0,yB-y2);
  g.fillStyle=css('--mute'); g.font='9.5px ui-monospace,monospace';
  g.textAlign='right'; g.textBaseline='middle';
  g.fillText('+'+rhi,x0-5,R(rhi*0.92)); g.fillText('0',x0-5,R(0)); g.fillText('-'+rhi,x0-5,R(-rhi*0.92));
  g.save(); g.translate(11,(y2+yB)/2); g.rotate(-Math.PI/2);
  g.textAlign='center'; g.textBaseline='middle';
  g.fillText('residual / '+String.fromCharCode(963),0,0); g.restore();
  g.textAlign='right'; g.textBaseline='top';
  g.fillText('pulse moment q',x1,yB+15);
  // top of the upper panel: the bottom-left corner already carries 'noise'
  g.fillStyle='rgba(255,255,255,.72)'; g.font='9.5px ui-monospace,monospace';
  g.textAlign='right'; g.textBaseline='top';
  g.fillText(F.haspred
    ? ('median '+String.fromCharCode(967)+String.fromCharCode(178)+'/N = '+F.chi2n.toFixed(2))
    : 'no kernel stored: measurements and noise only', x1-5, y0+4);
}
function drawGate(){
  const {g,w,h}=setup(document.getElementById('c3')); if(w<8) return;
  const xof=t=>PAD.l+((t.v-D.aem.lo)/(D.aem.hi-D.aem.lo))*(w-PAD.l-PAD.r);
  const A=frame(g,w,h,'conductivity S/m',LOGT,xof,'depth m');
  const X=sxLog(A);
  band(g,A);
  const MS = (useSmooth && D.gate.ms) ? D.gate.ms : D.gate.m;
  for(const k of admIdx){
    if(k===sel) continue;
    steps(g,A,X,MS[k].map(v=>Math.pow(10,v)),aemTop,aemBot,C_FAINT,1);
  }
  const bi=D.gate.C.indexOf(D.gate.best);
  if(bi>=0 && bi!==sel)
    steps(g,A,X,MS[bi].map(v=>Math.pow(10,v)),aemTop,aemBot,C_BEST,1.4,[5,3]);
  steps(g,A,X,D.aem.truth,aemTop,aemBot,C_TRUE,2,null,true);
  const m=MS[sel].map(v=>Math.pow(10,v));
  steps(g,A,X,m,aemTop,aemBot,C_SEL,3.2,null,true);
  for(const i of obOf(sel)){
    if(aemTop[i]>D.zmax) continue;
    const Y=(yOf(A,aemTop[i])+yOf(A,Math.min(aemBot[i],D.zmax)))/2;
    g.beginPath(); g.arc(X(m[i]),Y,4.4,0,6.284);
    g.fillStyle='rgba(0,0,0,.8)'; g.fill();
    g.lineWidth=1.8; g.strokeStyle=C_LIM; g.stroke();
  }
  veil(g,A);
  g.strokeStyle=css('--line'); g.lineWidth=1;
  g.strokeRect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0);
}
function drawChi(){
  const {g,w,h}=setup(document.getElementById('c4')); if(w<8) return;
  const C=D.gate.C, lC=C.map(Math.log10);
  const lo=Math.min(...lC)-0.08, hi=Math.max(...lC)+0.08;
  const y0=8,y1=h-26,x0=40,x1=w-8;
  const chis=D.gate.chi2.filter(isFinite);
  // The three rejected nodes reach chi2 49, 234 and 932. Scaling to them
  // flattens everything that matters; CLAMPING them to the top -- which this
  // did -- draws three different misfits as one horizontal line at a value
  // none of them has. Instead the axis covers the decision region, the curve
  // is clipped so it genuinely leaves the top, and each node above the axis
  // is marked with its own value.
  const cmin=Math.min(...chis)-1.5, cmax=D.tabs+8;
  const X=c=>x0+((Math.log10(c)-lo)/(hi-lo))*(x1-x0);
  const Y=v=>y1-((v-cmin)/(cmax-cmin))*(y1-y0);
  g.font='9.5px ui-monospace,monospace'; g.fillStyle=css('--mute');
  g.strokeStyle=css('--grid'); g.lineWidth=1;
  const TK=[0.01,0.02,0.05,0.1,0.2,0.5,1,2,5].filter(v=>Math.log10(v)>=lo&&Math.log10(v)<=hi);
  for(const v of TK){const x=X(v);
    g.beginPath();g.moveTo(x,y0);g.lineTo(x,y1);g.stroke();
    g.textAlign='center';g.textBaseline='top';
    g.fillText(v<1?String(v).replace(/^0/,''):String(v),x,y1+4);}
  g.strokeStyle=css('--hot'); g.setLineDash([5,3]); g.beginPath();
  g.moveTo(x0,Y(D.tabs)); g.lineTo(x1,Y(D.tabs)); g.stroke();
  g.strokeStyle=css('--acc2'); g.beginPath();
  g.moveTo(x0,Y(D.trel)); g.lineTo(x1,Y(D.trel)); g.stroke(); g.setLineDash([]);
  g.textAlign='left'; g.textBaseline='bottom'; g.fillStyle=css('--hot');
  g.fillText('N + 3.84 = '+F(D.tabs,2), x0+3, Y(D.tabs)-1);
  g.fillStyle=css('--acc2');
  g.fillText('best fit + 3.84 = '+F(D.trel,2), x0+3, Y(D.trel)-1);
  g.strokeStyle=css('--grid'); g.textAlign='right'; g.textBaseline='middle';
  g.fillStyle=css('--mute');
  for(const v of [30,40,50]){ if(v<cmin||v>cmax) continue;
    g.beginPath(); g.moveTo(x0,Y(v)); g.lineTo(x1,Y(v)); g.stroke();
    g.fillText(String(v),x0-4,Y(v)); }
  g.strokeStyle=css('--line'); g.beginPath();
  g.moveTo(x0,y0); g.lineTo(x0,y1); g.lineTo(x1,y1); g.stroke();
  g.fillStyle=css('--mute'); g.textAlign='right'; g.textBaseline='top';
  g.fillText('window conductance S',x1,y1+14);
  g.save(); g.translate(11,(y0+y1)/2); g.rotate(-Math.PI/2);
  g.textAlign='center'; g.textBaseline='middle'; g.fillText('min '+String.fromCharCode(967)+String.fromCharCode(178),0,0); g.restore();
  g.save(); g.beginPath(); g.rect(x0,y0-2,x1-x0,y1-y0+2); g.clip();
  g.strokeStyle=css('--acc'); g.lineWidth=1.8; g.beginPath();
  C.forEach((c,i)=>{const px=X(c),py=Y(D.gate.chi2[i]); i?g.lineTo(px,py):g.moveTo(px,py);});
  g.stroke();
  C.forEach((c,i)=>{ if(D.gate.chi2[i]>cmax) return;
    g.beginPath();g.arc(X(c),Y(D.gate.chi2[i]),2.6,0,6.284);
    g.fillStyle=D.gate.adm[i]?css('--acc'):css('--mute');g.fill();});
  g.restore();
  // nodes off the top, each with its own value
  g.fillStyle=css('--mute'); g.font='9px ui-monospace,monospace';
  g.textAlign='center'; g.textBaseline='top';
  C.forEach((c,i)=>{ if(D.gate.chi2[i]<=cmax) return;
    const px=X(c);
    g.beginPath(); g.moveTo(px,y0+1); g.lineTo(px-3.4,y0+6); g.lineTo(px+3.4,y0+6);
    g.closePath(); g.fill();
    g.fillText(String(Math.round(D.gate.chi2[i])),px,y0+7);});
  g.strokeStyle=css('--hot'); g.lineWidth=1.3; g.setLineDash([2,2]); g.beginPath();
  g.moveTo(X(D.truthC),y0); g.lineTo(X(D.truthC),y1); g.stroke(); g.setLineDash([]);
  // the slider's position, as a marker that moves along the curve
  const sc=D.gate.C[sel], sy=Y(D.gate.chi2[sel]), sx=X(sc);
  g.strokeStyle=css('--ink'); g.globalAlpha=.45; g.lineWidth=1.2;
  g.beginPath(); g.moveTo(sx,y0); g.lineTo(sx,y1); g.stroke();
  g.beginPath(); g.moveTo(x0,sy); g.lineTo(sx,sy); g.stroke(); g.globalAlpha=1;
  g.beginPath(); g.arc(sx,sy,7.5,0,6.284);
  g.fillStyle=C_SEL; g.fill();
  g.lineWidth=2.4; g.strokeStyle=css('--panel'); g.stroke();
  g.lineWidth=1.2; g.strokeStyle=css('--ink'); g.beginPath(); g.arc(sx,sy,8.8,0,6.284); g.stroke();
  g.font='10px ui-monospace,monospace'; g.fillStyle=css('--ink');
  g.textAlign = sx>(x0+x1)/2 ? 'right' : 'left';
  g.textBaseline='bottom';
  g.fillText(String.fromCharCode(967)+String.fromCharCode(178)+' '+F(D.gate.chi2[sel],2),
             sx+(sx>(x0+x1)/2?-13:13), sy+4);
}

// ---- conductance over a depth range, for one model (S/m per cell)
function cond(m,a,b){let s=0;
  for(let i=0;i<m.length;i++){
    const ov=Math.max(0,Math.min(aemBot[i],b)-Math.max(aemTop[i],a));
    if(ov>0) s+=m[i]*ov;}
  return s;}

const admIdx = D.gate.C.map((c,i)=>i).filter(i=>D.gate.adm[i]);
let useSmooth=false;
// the stored `ob` lists belong to the min-chi2 models; for the smoothest ones
// recompute with the same predicate the gate uses
const obOf = k => (useSmooth && D.gate.ms)
  ? D.gate.ms[k].map((v,i)=>[v,i]).filter(([v])=>v<=D.mb[0]+1e-9||v>=D.mb[1]-1e-9).map(([,i])=>i)
  : D.gate.ob[k];
// open on the node nearest the geometric centre of the admissible range: a
// well-fitting model in the middle of what the data allow, picked by the
// range alone and not by where the answer happens to be
const admMid = Math.sqrt(D.gate.C[admIdx[0]]*D.gate.C[admIdx[admIdx.length-1]]);
let sel = admIdx.reduce((a,b)=>
  Math.abs(Math.log(D.gate.C[b]/admMid))<Math.abs(Math.log(D.gate.C[a]/admMid))?b:a, admIdx[0]);

// Conductance above and below the window, as a factor of the lowest-misfit
// model: the compensation the forced conductance is paid for with. Exported
// per node rather than shown, so the table stays still.
const ZB=Math.min(D.zmax,80.7);
function ratios(k){
  const SRC=(useSmooth && D.gate.ms)?D.gate.ms:D.gate.m;
  const bi=D.gate.C.indexOf(D.gate.best);
  const mk=SRC[k].map(v=>Math.pow(10,v)), bm=SRC[bi].map(v=>Math.pow(10,v));
  return [cond(mk,0,D.win[0])/cond(bm,0,D.win[0]),
          cond(mk,D.win[1],ZB)/cond(bm,D.win[1],ZB)];
}
function say(){
  const C=D.gate.C[sel], chi=D.gate.chi2[sel];
  const SRC=(useSmooth && D.gate.ms)?D.gate.ms:D.gate.m;
  const mm=SRC[sel].map(v=>Math.pow(10,v));
  const th=+document.getElementById('sP').value, ex=+document.getElementById('sM').value;
  const kf=+document.getElementById('sK').value/100;
  const sigb=C/D.thick, clay=kf*sigb;
  // sigma_w = (sigma_bulk - clay) / theta^m, the relation Stage 5 uses. Clay
  // can only ADD conductance, so clay = 0 gives the MOST saline reading the
  // data allow; every step of this slider lowers it.
  const tds=(sigb-clay)/Math.pow(th,ex)*D.pet.mgl;
  const ceilK=D.pet.scmax*D.pet.c0*Math.max(0,1-th/D.pet.phiref);
  const overK=clay>Math.min(ceilK,D.pet.cp[1]);
  const ob=obOf(sel);
  const row=(k,v,t)=>'<tr><td class="k">'+k+'</td><td class="v">'+v+
    '</td><td class="t">'+(t===undefined?'&middot;':t)+'</td></tr>';
  document.getElementById('say').innerHTML =
    '<colgroup><col class="k"><col class="v"><col class="t"></colgroup>'+
    '<thead><tr><th class="k"></th><th>selected</th><th>true</th></tr></thead><tbody>'+
    row('window conductance', F(C,4)+' S', F(D.truthC,4)+' S')+
    row('misfit', String.fromCharCode(967)+String.fromCharCode(178)+' '+F(chi,2))+
    row('on a limit', ob.length===0 ? '0 <em>of '+mm.length+'</em>'
        : ob.length+' <em>of '+mm.length+', from '+F(aemTop[ob[0]],1)+' m</em>')+
    row('clay term', F(clay,4)+' S/m <em>'+Math.round(kf*100)+'%'+
        (overK?', past the ceiling':'')+'</em>',
        F(D.claytrue.term,4)+' S/m  '+Math.round(D.claytrue.frac*100)+'%')+
    row('pore water', Math.round(tds).toLocaleString()+' mg/L'+
        (kf===0?' <em>most saline allowed</em>':''),
        Math.round(D.claytrue.tds).toLocaleString()+' mg/L')+
    '</tbody>';
  document.getElementById('oC').textContent=F(C,4)+' S';
  document.getElementById('oP').textContent=F(th,3);
  document.getElementById('oM').textContent=F(ex,2);
  document.getElementById('oK').textContent=Math.round(kf*100)+'% / '+F(clay,4);
}

function csv(){
  const L=['block,key1,key2,value'];
  const bw1=(D.aem.hi-D.aem.lo)/D.aem.nb;
  for(let i=0;i<D.aem.h.length;i++)for(let j=0;j<D.aem.nb;j++){const c=D.aem.h[i][j];
    if(c) L.push('aem_hist,'+F(D.aem.z[i],3)+','+F(D.aem.lo+(j+0.5)*bw1,4)+','+c);}
  for(let i=0;i<D.aem.z.length;i++) L.push('aem_lines,'+F(D.aem.z[i],3)+',p05_p50_p95_truth,'+
    F(D.aem.p05[i],6)+' '+F(D.aem.p50[i],6)+' '+F(D.aem.p95[i],6)+' '+F(D.aem.truth[i],6));
  const bw2=(D.snmr.hi-D.snmr.lo)/D.snmr.nb;
  for(let i=0;i<D.snmr.h.length;i++)for(let j=0;j<D.snmr.nb;j++){const c=D.snmr.h[i][j];
    if(c) L.push('snmr_hist,'+F(D.snmr.z[i],3)+','+F(D.snmr.lo+(j+0.5)*bw2,4)+','+c);}
  for(let i=0;i<D.snmr.z.length;i++) L.push('snmr_lines,'+F(D.snmr.z[i],3)+',p05_p50_p95_truth,'+
    F(D.snmr.p05[i],6)+' '+F(D.snmr.p50[i],6)+' '+F(D.snmr.p95[i],6)+' '+F(D.snmr.truth[i],6));
  for(let k=0;k<D.gate.C.length;k++){
    L.push('gate_node,'+F(D.gate.C[k],6)+',chi2_admissible,'+F(D.gate.chi2[k],6)+' '+D.gate.adm[k]);
    if(D.gate.adm[k]){ const r=ratios(k);
      L.push('gate_ratio,'+F(D.gate.C[k],6)+',above_0_'+F(D.win[0],1)+'_below_'+
             F(D.win[1],1)+'_'+F(ZB,1)+',' +F(r[0],6)+' '+F(r[1],6)); }
    const ob=new Set(D.gate.ob[k]);
    for(let i=0;i<D.gate.m[k].length;i++)
      L.push('gate_model,'+F(D.gate.C[k],6)+','+F(D.aem.z[i],3)+','+
        F(Math.pow(10,D.gate.m[k][i]),6)+' '+(ob.has(i)?1:0));}
  return L.join('\n');
}

let tab2='wc';
function draw2(){ tab2==='wc' ? drawSNMR() : drawFit(); }
function redraw(){drawAEM();draw2();drawGate();drawChi();say();}

const sC=document.getElementById('sC');
sC.max=String(admIdx.length-1);
sC.value=String(admIdx.indexOf(sel));
sC.addEventListener('input',()=>{sel=admIdx[+sC.value];drawGate();drawChi();say();});
(function(){
  const P=document.getElementById('sP'), M=document.getElementById('sM');
  P.min=D.pet.phi[0]; P.max=D.pet.phi[1]; P.value=(D.pet.phi[0]+D.pet.phi[1])/2;
  M.min=D.pet.m[0];  M.max=D.pet.m[1];  M.value=(D.pet.m[0]+D.pet.m[1])/2;
})();
for(const id of ['sP','sM','sK']) document.getElementById(id).addEventListener('input',say);
(function(){
  const b=document.getElementById('bM'), h=document.getElementById('hM');
  const have=!!D.gate.ms;
  b.textContent='lowest misfit';
  b.disabled=!have;
  if(have){ h.textContent='';
    b.addEventListener('click',()=>{
      useSmooth=!useSmooth; b.setAttribute('aria-pressed',String(useSmooth));
      b.textContent=useSmooth?'smoothest passing':'lowest misfit';
      drawGate(); say();});
  } else {
    b.title='the smoothest admissible model per node is not stored yet';
    h.textContent=D.gate.msn>0 ? D.gate.msn+'/'+D.gate.msneed+' nodes' : 'not stored yet';
  }
})();

document.getElementById('bC').addEventListener('click',()=>{
  const text=csv(), ta=document.getElementById('ta');
  ta.value=text; ta.hidden=false; document.getElementById('bS').hidden=false;
  try{
    const u=URL.createObjectURL(new Blob([text],{type:'text/csv'}));
    const a=document.createElement('a'); a.href=u; a.download='report_3P_data.csv';
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(()=>URL.revokeObjectURL(u),2000);
  }catch(e){}
});
document.getElementById('bS').addEventListener('click',()=>{
  const ta=document.getElementById('ta'); ta.focus(); ta.select();});

document.getElementById('prov').innerHTML =
  'AEM '+D.aem.src+' &times;'+D.aem.nens+', '+D.aem.ndraw.toLocaleString()+' draws &middot; '+
  'SNMR '+D.snmr.src+', '+D.snmr.ndraw.toLocaleString()+' draws <span class="flag">'+D.snmr.note+'</span> &middot; '+
  'gate '+D.gate.src+', '+Math.pow(10,D.mb[0]).toPrecision(2)+'&ndash;'+
  Math.pow(10,D.mb[1]).toPrecision(2)+' S/m &middot; '+D.gen;
document.getElementById('n1').textContent=D.aem.ndraw.toLocaleString()+' draws';
function setTab(which){
  tab2=which;
  const a=document.getElementById('t2a'), b=document.getElementById('t2b');
  a.classList.toggle('on',which==='wc'); b.classList.toggle('on',which==='fit');
  a.setAttribute('aria-pressed',String(which==='wc'));
  b.setAttribute('aria-pressed',String(which==='fit'));
  document.getElementById('n2').textContent = which==='wc'
    ? D.snmr.ndraw.toLocaleString()+' draws'
    : (D.fit ? D.fit.q.length+' gates, '+D.fit.src : 'no data file');
  draw2();
}
document.getElementById('t2a').addEventListener('click',()=>setTab('wc'));
document.getElementById('t2b').addEventListener('click',()=>setTab('fit'));
if(!D.fit) document.getElementById('t2b').disabled=true;
document.getElementById('n2').textContent=D.snmr.ndraw.toLocaleString()+' draws';
document.getElementById('n3').textContent=admIdx.length+' of '+D.gate.C.length+' nodes';
document.getElementById('key').innerHTML=
  '<span><i style="background:'+C_P50+';outline:1px solid '+css('--rule')+'"></i>p50</span>'+
  '<span><i style="background:'+C_PTL+';outline:1px solid '+css('--rule')+'"></i>p05 p95</span>'+
  '<span><i style="background:'+C_TRUE+';outline:1px solid '+css('--rule')+'"></i>true</span>'+'<span><i style="background:'+C_BEST+'"></i>best fit</span>'+'<span><i style="background:rgba(255,255,255,.65)"></i>sensitivity depth</span>'+
  '<span><i style="background:transparent;border:1.6px solid '+C_LIM+';outline:1px solid '+css('--rule')+';height:8px;width:8px;border-radius:50%"></i>on a stated limit</span>'+
  '<span>window '+D.win[0]+'&ndash;'+D.win[1]+' m</span>';

redraw();
let t=null;
addEventListener('resize',()=>{clearTimeout(t);t=setTimeout(redraw,120);});
try{
  const mq=matchMedia('(prefers-color-scheme:dark)');
  const onTheme=()=>setTimeout(redraw,30);
  mq.addEventListener ? mq.addEventListener('change',onTheme) : mq.addListener(onTheme);
}catch(e){}
})();
</script>
"""

mkpath(dirname(OUTF))
# PAGE carries <title> + <style> and then the markup. For a standalone file
# the first part belongs in <head> and the rest in <body>, so it is split at
# the end of the style block rather than emitted as one run of tag soup.
let i = findfirst("</style>", PAGE)
    i === nothing && error("PAGE has no </style>; cannot split head from body")
    headpart = PAGE[1:last(i)]
    bodypart = PAGE[last(i)+1:end]
    # Write and VALIDATE the JSON before the page, so the page can only ever
    # embed a payload that parsed.
    JSONSTR = json_payload()
    JSONF = replace(OUTF, r"\.html$" => ".json")
    open(JSONF, "w") do f; print(f, JSONSTR); end
    let chk = read(`python3 -c "import json,sys; json.load(open(sys.argv[1])); print('ok')" $JSONF`, String)
        strip(chk) == "ok" || error("report JSON did not parse: $JSONF")
    end
    @printf("wrote %s  (%.1f kB, valid JSON)\n", relpath(JSONF, ROOT), filesize(JSONF)/1024)

    open(OUTF, "w") do f
        if !FRAGMENT
            print(f, "<!doctype html>\n<html lang=\"en\">\n<head>\n",
                     "<meta charset=\"utf-8\">\n",
                     "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n")
            print(f, headpart, "\n</head>\n<body>\n", bodypart)
        else
            print(f, PAGE)
        end
        print(f, "\n<script>\nconst D=", JSONSTR, ";\n</script>\n")
        print(f, SCRIPT)
        FRAGMENT || print(f, "\n</body>\n</html>\n")
    end
end
@printf("\nwrote %s  (%.1f KB, %s)\n", OUTF, filesize(OUTF) / 1024,
        FRAGMENT ? "body fragment" : "standalone document")
