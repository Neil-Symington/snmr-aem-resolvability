# detectability_cellsize.jl -- does the constraint depth measure the data or the grid?
#
#   julia --project=. fieldrun/test/detectability_cellsize.jl --config <toml> \
#         [--truth 3P=data/truth_murray3P.jld2,3S=...] [--sites Mulan31,Mulan32] \
#         [--levels 25,50,100] [--kcache <dir>] [--save <file.jld2>]
#
# Rebuilds each SMR kernel on the SNMR depth grid as it stands and with every
# cell split into 2 and into 4, and evaluates on each:
#
#   cell     detectability()          one cell perturbed alone (the current rule)
#   blk10    detectability_block()    a 10 m block, THICK_MIN
#   blk5     detectability_block()    a 5 m block, [detect] min_thick_m in the BDU files
#   sdoi     sens_doi(w_ref = 0.05)   the cumulative DOI the stage reports (quadrature)
#   sdoiC    sens_doi_coherent        the same, section summed coherently over depth
#   blkP     detectability_projected  blk10 with the kernel projected onto the response
#                                     phase at w_ref = 0.05, as sdoiC does: blk10 minus
#                                     its no-cancellation |sum K̃|
#
# A rule that measures the data gives the same depth at every refinement. The
# functions are lifted from detect_aquifer.jl by name, not copied, so this tests
# the code the stage runs.
#
# SYNTHETIC COLUMNS (--truth) use the truth's sigma_bulk for the kernel and the
# noise tools/smr_data.jl generates on m1-production-rerun: the EK-measured
# shape, median 4.63 nV (typical) and 11.73 nV (poor), on its 24-moment q grid.
# That is s = 1, the synthetics' convention. FIELD SITES (--sites) read the
# SNMR trans-D output under [paths] outputs and use both sigma-hat and the
# inferred sqrt(s50) * sigma-hat, on the kernel detect_aquifer.jl builds.
#
# Needs truth files and field outputs that are not in the repository, so it is
# a script, not part of test/runtests.jl.
using JLD2, Statistics, Printf, Random, LinearAlgebra

include(joinpath(@__DIR__, "..", "common", "stage_prelude.jl"))
include(joinpath(COMMON, "prerun_data.jl"))
include(joinpath(COMMON, "udf_smr.jl"))

# ── the functions under test, taken from the stage by name ──────────────────
const DETECT_SRC = joinpath(@__DIR__, "..", "stages", "detect_aquifer.jl")
let want = Set([:cellsof, :detectability, :detectability_block, :sens_doi, :sens_doi_coherent,
                :THICK_MIN, :THETA_DETECT])
    ex = Meta.parseall(read(DETECT_SRC, String))
    name(e) = e isa Expr && e.head === :macrocall ? name(e.args[end]) :
              e isa Expr && e.head === :const ? name(e.args[1]) :
              e isa Expr && e.head in (:function, :(=)) ?
                  (a = e.args[1]; a isa Expr && a.head === :call ? a.args[1] : a) : nothing
    got = Set{Symbol}()
    for e in ex.args
        n = name(e)
        n in want || continue
        Core.eval(Main, e); push!(got, n)
    end
    got == want || error("not found in detect_aquifer.jl: $(setdiff(want, got))")
end

# ── grids and kernels ───────────────────────────────────────────────────────
"Split every interval of the node grid into k equal parts."
refine(z, k) = k == 1 ? collect(Float64, z) :
    vcat([collect(range(z[i], z[i+1], length = k + 1))[1:end-1] for i in 1:length(z)-1]..., z[end])

"K̃ exactly as SMRPInversion.SMRDetOperator forms it: trapezoidal weights times the kernel."
function ktilde(F)
    dz = F.dz; nz = size(F.kernel, 1)
    wt = zeros(nz); wt[1] = 0.5dz[1]; wt[end] = 0.5dz[end]
    for i in 2:nz-1; wt[i] = 0.5 * (dz[i-1] + dz[i]); end
    permutedims(wt .* F.kernel)
end

# smr_data.jl (m1-production-rerun), verbatim apart from layout
const QGRID = collect(10.0 .^ range(log10(0.1), log10(15.0), length = 24))
# --levels 25,50,100 (nV) adds levels with the same shape. At the two tiers the
# synthetic columns resolve the whole 147 m grid under every rule, so they
# cannot show a cell-size effect; field-like levels put the limit in the column.
const TIERS = vcat([("typical", 4.63e-9), ("poor", 11.73e-9)],
    [("level", parse(Float64, x) * 1e-9) for x in split(get(ARGN, "levels", ""), ",") if !isempty(x)])
function sigma_shape(nq, level; rise = 2.5, cv = 0.60, seed = 20261012)
    ramp = collect(range(1.0, rise, length = nq))
    f(w) = (Random.seed!(seed); s = ramp .* exp.(w .* randn(nq)); std(s) / mean(s))
    lo, hi = 0.0, 2.0
    for _ in 1:60; mid = (lo + hi) / 2; f(mid) < cv ? (lo = mid) : (hi = mid); end
    Random.seed!(seed)
    s = ramp .* exp.(((lo + hi) / 2) .* randn(nq))
    s .* (level / median(s))
end

"""
detectability_block with sens_doi_coherent's phase convention: the block's data
change is the component along the current response, real(conj(u) .* K̃) summed
over the block, so depths whose contributions are out of phase cancel. Only the
projection differs from blk10; the block, threshold and base are the same.
"""
function detectability_projected(K̃, σ, z; w_ref = 0.05, θdetect = THETA_DETECT, Lmin = THICK_MIN)
    zc = K̃ * fill(w_ref, size(K̃, 2))
    u = conj.(zc) ./ max.(abs.(zc), eps())
    detectability_block(real.(u .* K̃), σ, z; θdetect = θdetect, Lmin = Lmin)
end

"fig04_snmr_posteriors.jl's sensitivity DOI (m1): the cell rule, on the same K̃."
fig04_doi(K̃, σ, z) = detectability(K̃, σ, z).zlim

const REFINE = (1, 2, 4)
row(label, f, K̃, σ, z) = (label = label, z = z,
    cell = detectability(K̃, σ, z).zlim,
    blk10 = detectability_block(K̃, σ, z; Lmin = 10.0).zlim,
    blk5 = detectability_block(K̃, σ, z; Lmin = 5.0).zlim,
    sdoi = sens_doi(K̃, σ, z, 0.05).doi,
    sdoiC = sens_doi_coherent(K̃, σ, z, 0.05).doi,
    blkP = detectability_projected(K̃, σ, z; Lmin = 10.0).zlim)

function report(title, rows)
    println("\n", title)
    @printf("  %-22s %6s %5s | %8s %8s %8s %8s %8s %8s\n", "noise", "refine", "dz0", "cell", "blk10", "blk5", "sdoi", "sdoiC", "blkP")
    for (lab, rs) in rows, (k, r) in zip(REFINE, rs)
        @printf("  %-22s %6s %5.2f | %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f\n", k == 1 ? lab : "",
                "x$k", r.z[2] - r.z[1], r.cell, r.blk10, r.blk5, r.sdoi, r.sdoiC, r.blkP)
    end
    for (lab, rs) in rows
        a, b = rs[1], rs[2]
        ratio(x, y) = x == 0 ? (y == 0 ? "0->0" : "0->>0") : @sprintf("%.3f", y / x)
        @printf("  %-22s half/current: cell %s  blk10 %s  blk5 %s  sdoi %s  sdoiC %s  blkP %s\n", lab,
                ratio(a.cell, b.cell), ratio(a.blk10, b.blk10), ratio(a.blk5, b.blk5),
                ratio(a.sdoi, b.sdoi), ratio(a.sdoiC, b.sdoiC), ratio(a.blkP, b.blkP))
    end
end

ALL = Dict{String,Any}()

# ── synthetic columns ───────────────────────────────────────────────────────
if haskey(ARGN, "truth")
    for spec in split(ARGN["truth"], ",")
        col, path = split(spec, "=")
        T = load(path)
        zt = collect(Float64, T["z"]); σt = collect(Float64, T["σ_bulk"])
        # One conductivity model serves every grid: the truth's 0.25 m sigma_bulk,
        # geometric-meaned onto 2 m layers. NOT the 0.25 m layers themselves:
        # SNMRForward overflows silently on a fine stack -- 800 x 0.25 m gave
        # |kernel| ~ 1e55 and V0 ~ 1e62 nV, 200 x 1 m ~ 1e11 nV -- while 2 m
        # and 5 m layers agree (3P: V0 356-2920 nV). Every grid passed every
        # rule on the broken kernel, so check V0 below before trusting a row.
        bnd = collect(0.0:2.0:zt[end] + (zt[2] - zt[1]) / 2)
        σl = [exp(mean(log.(σt[(zt .>= bnd[i]) .& (zt .< bnd[i+1])]))) for i in 1:length(bnd)-1]
        cm = SNMRForward.ConductivityModel(σl, diff(bnd)[1:end-1])
        R = T["LOOP_SIDE"] / 2; ϕ = deg2rad(T["INC_DEG"]); Be = T["BE_T"]
        rows = Pair{String,Vector}[]
        # kernels are minutes each on the 0.25 m truth layers; --kcache keeps them
        kc = haskey(ARGN, "kcache") ? joinpath(ARGN["kcache"], "K_$(col).jld2") : ""
        Ks = if !isempty(kc) && isfile(kc)
            load(kc, "Ks")
        else
            k_ = [ktilde(SMRPInversion.MRSForward_square(R, refine(T["zsmr"], k), QGRID, ϕ, 0.0, Be, cm))
                  for k in REFINE]
            isempty(kc) || (mkpath(dirname(kc)); jldsave(kc; Ks = k_))
            k_
        end
        V0 = abs.(Ks[1] * T["θ_smrgrid"])
        1e-8 < maximum(V0) < 1e-4 ||
            error("$col: forward V0 up to $(maximum(V0)) V is not an SNMR amplitude -- the kernel is broken")
        @printf("\n%s forward V0 on the truth: %.0f-%.0f nV\n", col, 1e9 * minimum(V0), 1e9 * maximum(V0))
        for (tier, level) in TIERS
            σ = sigma_shape(length(QGRID), level)
            push!(rows, "$tier $(round(level*1e9, digits=2)) nV" =>
                  [row(tier, nothing, Ks[i], σ, refine(T["zsmr"], k)) for (i, k) in enumerate(REFINE)])
        end
        report("SYNTHETIC $col  (truth sigma_bulk kernel, loop $(T["LOOP_SIDE"]) m, inc $(T["INC_DEG"]))", rows)
        ALL[String(col)] = rows
    end
end

# ── field sites ─────────────────────────────────────────────────────────────
if haskey(ARGN, "sites")
    for fid in split(ARGN["sites"], ",")
        fid = String(fid)
        D = load(joinpath(OUT, "smr_transd_$(fid).jld2"))
        r = udf_meta(fid)
        b = build_udf_sounding(merge(r, (matpath = udf_matfile(fid)[1],)))
        cm = haskey(D, "aem_sigma") ?
            SNMRForward.ConductivityModel(D["aem_sigma"], diff(D["aem_ztop"])) :
            (A = load(joinpath(OUT, "aem_occam_$(fid).jld2"));
             SNMRForward.ConductivityModel(A["sigma"], diff(A["zbnd"])))
        Ks = [ktilde(SMRPInversion.MRSForward_square(r.loop_side / 2, refine(b.zbnd, k), b.q,
                     deg2rad(r.inc_deg), 0.0, b.Be, cm)) for k in REFINE]
        # the unrefined kernel must be the stage's own, or this tests something else
        Kstage = SMRPInversion.SMRDetOperator(sounding_fixed(merge(b, (F = SMRPInversion.MRSForward_square(
                     r.loop_side / 2, b.zbnd, b.q, deg2rad(r.inc_deg), 0.0, b.Be, cm),)))).K̃
        maximum(abs.(Ks[1] .- Kstage)) <= 1e-12 * maximum(abs.(Kstage)) ||
            error("$fid: rebuilt K̃ differs from the stage's")
        σ̂ = D["sigma_hat"]; s50 = quantile(exp.(D["log_s"]), 0.5)
        rows = Pair{String,Vector}[]
        for (lab, σ) in (("sigma-hat (s = 1)", σ̂), (@sprintf("inferred, s50 %.0f", s50), sqrt(s50) .* σ̂))
            push!(rows, lab => [row(lab, nothing, Ks[i], σ, refine(b.zbnd, k)) for (i, k) in enumerate(REFINE)])
        end
        report("FIELD $fid  ($(D["cond_source"]), sigma-hat median $(round(median(σ̂)*1e9, digits=2)) nV)", rows)
        ALL[fid] = rows
    end
end

if haskey(ARGN, "save")
    jldsave(ARGN["save"]; results = ALL)
    println("\nwrote ", ARGN["save"])
end
