# conditions.jl -- Step 5 replacement: a conditions table per salinity class.
#
#   julia --project=. conditions_table/conditions.jl <TAG> [<TAG> ...]
#
# No inversions. Reads existing gate output, existing SNMR ensembles and the
# config. Writes only into conditions_table/.
#
# THE RELATION, exactly as stages/stage5_quality.jl:17-18 uses it:
#
#     sigma_bulk = sigma_w * theta_total^m + clay_term        (clay_term = sigma_clay * C)
#     sigma_w    = (sigma_bulk - clay_term) / theta_total^m
#     TDS        = sigma_w * mg_per_sm,     mg_per_sm = 6400
#
# so for a target TDS the clay term a combination REQUIRES is
#
#     clay_term = sigma_bulk - (TDS / 6400) * theta_total^m
#
# which is decreasing in TDS: a saltier water needs LESS clay to explain the
# same conductance. A class is RULED OUT when no combination inside the limits
# leaves a clay term that is both >= 0 and <= the ceiling.
using JLD2, TOML, Printf, Statistics

const HERE = normpath(@__DIR__)
const ROOT = normpath(joinpath(HERE, ".."))
const PARENT = normpath(joinpath(ROOT, ".."))
include(joinpath(ROOT, "src", "config.jl")); using .SiteConfig

const TAGS = isempty(ARGS) ? ["2P","3P","2S","3S"] : ARGS
const MGL  = 6400.0
# Classes as the BRIEF defines them. The code's own hypothesis threshold is a
# single 1000 mg/L ([hypothesis.salinity].value); see the report.
const CLASSES = [("fresh", 0.0, 1500.0), ("brackish", 1500.0, 3000.0),
                 ("saline", 3000.0, Inf)]

cfg0 = SiteConfig.resolve(joinpath(ROOT, "configs", "synth_2P.toml"))
const PR = cfg0["priors"]; const CE = cfg0["ceilings"]
const PHI_LO, PHI_HI = Float64.(PR["porosity_range"])
const M_LO, M_HI = 1.5, 1.9                 # clayey_sand, from [priors].m_lithology
const CP_LO, CP_HI = Float64.(PR["clay"]["product"])
const SCMAX = Float64(CE["sigma_clay_max"]); const C0 = Float64(CE["clay_c0"])
const PHIREF = Float64(CE["clay_phi_ref"])

clay_ceiling(θ) = SCMAX * C0 * max(1 - θ/PHIREF, 0.0)

# ── per-column inputs, all from existing output ────────────────────────────
WINDOW = Dict("2P"=>(17.5,33.5), "2S"=>(17.5,33.5), "3P"=>(7.5,23.5), "3S"=>(7.5,23.5))
# gate band in CONDUCTIVITY (S/m) and whether each edge is open, from stage4.
# 2P's 1 Nov stage4 predates the edge_*_dex keys, so its band comes from the
# verdict table (0.03192-7.079 S over a 16 m window).
GATE = Dict(
 "2P"=>(lo=0.03192/16, hi=7.079/16, openlo=true,  openhi=false, src="production 1-2 Nov, verdict table"),
 "3P"=>(lo=10^-2.7333333333333334, hi=10^-1.6666666666666667, openlo=true, openhi=true,
        src="production 2 Nov, stage4_synth_3P.jld2"),
 # 2S and 3S stage4 files predate the edge_*_dex keys too; bands from their
 # verdict tables. 3S's band is EXACTLY 1.000 dex wide (0.1271-1.271 S, a
 # factor of 10) with both edges open, which is the shape of a band set by the
 # search bounds rather than by the data -- treat its numbers with suspicion.
 "2S"=>(lo=0.03192/16, hi=7.185/16, openlo=true, openhi=false, src="production 2 Nov, verdict table"),
 "3S"=>(lo=0.1271/16,  hi=1.271/16, openlo=true, openhi=true,  src="production 2 Nov, verdict table (SUSPECT: exactly a factor of 10 wide)"),
)

"Window-mean water content per ensemble, from the trans-D SNMR ensemble."
function thetabyens(tag)
    # 2P has two vintages. The 13 Nov production ensemble (200,001 samples)
    # recovers the truth; the 22 Oct one (100,001 samples) sits 25% BELOW it on
    # all four ensembles, which is systematic, not sampling noise. Prefer the
    # newer where it exists and say so. The other three columns have only the
    # 22 Oct vintage, so their water content carries that same suspicion.
    cands = [joinpath(ROOT,"..","..","..","tmp_none"),
             "/tmp/claude-1000/-home-symo-Documents-GA-AEM-SMR-Resolvability-Quickrun/b6c90248-48ae-47d6-b7b7-ce3626a90b2e/scratchpad/prod/data/cond_C2_typical$(tag)T.jld2",
             joinpath(PARENT,"synthetics","outputs","HYP_SYNTH","cond_C2_typical$(tag)T.jld2")]
    f = findfirst(isfile, cands)
    f === nothing && return (ens=Vector{Vector{Float64}}(), src="NOT FOUND", path="")
    D = load(cands[f]); E = D["ens"]; z = D["zsmr"]
    w = WINDOW[tag]; idx = findall(q -> w[1] <= q <= w[2], z)
    (ens=[vec(mean(view(M,idx,:),dims=1)) for M in E], src=basename(cands[f]), path=cands[f])
end

function truthof(tag)
    T = load(joinpath(ROOT,"data","truth_murray$tag.jld2"))
    z = T["z"]; w = WINDOW[tag]; k = findall(q -> w[1] <= q <= w[2], z)
    ks = findall(q -> w[1] <= q <= w[2], T["zsmr"])
    (σw=mean(T["σw"][k]), tds=mean(T["σw"][k])*MGL, clay=mean((T["σ_clay"].*T["clay"])[k]),
     θtot=mean(T["θ_total"][k]), θsmr=mean(T["θ_smrgrid"][ks]), σb=mean(T["σ_bulk"][k]))
end

# ── the search ─────────────────────────────────────────────────────────────
"""
For one class, sweep the grid and return the range of clay term it requires,
plus which limit binds at each end.
`σrange`, `θrange`, `mrange` are (lo,hi). `capclay` optionally restricts the
clay term further (the gamma-log case in Part 2).
"""
function classrange(σlo, σhi, θlo, θhi, mlo, mhi, T1, T2; capclay = nothing, n = 41)
    σs = exp.(range(log(σlo), log(σhi), length = n))
    θs = range(θlo, θhi, length = n)
    ms = range(mlo, mhi, length = 17)
    best_lo = Inf; best_hi = -Inf; whylo = ""; whyhi = ""
    for σ in σs, θ in θs, m in ms
        θm = θ^m
        ceil_ = clay_ceiling(θ)
        capclay !== nothing && (ceil_ = min(ceil_, capclay[2]))
        floor_ = capclay === nothing ? 0.0 : max(0.0, capclay[1])
        ceil_ = min(ceil_, CP_HI)
        lo_raw = isfinite(T2) ? σ - (T2/MGL)*θm : -Inf     # saltier -> less clay
        hi_raw = σ - (T1/MGL)*θm
        lo = max(lo_raw, floor_); hi = min(hi_raw, ceil_)
        lo > hi && continue                                 # this combination fails
        if lo < best_lo
            best_lo = lo
            whylo = lo_raw >= floor_ ? (σ ≈ σlo ? "gate lower edge" : "class upper TDS") :
                    (capclay === nothing ? "clay term >= 0" : "clay cap (low)")
        end
        if hi > best_hi
            best_hi = hi
            whyhi = hi_raw <= ceil_ ? (σ ≈ σhi ? "gate upper edge" : "class lower TDS") :
                    (ceil_ ≈ CP_HI ? "clay prior ceiling 1.0" : "clay ceiling (theta)")
        end
    end
    best_hi < best_lo ? (supported=false, lo=NaN, hi=NaN, whylo="", whyhi="") :
                        (supported=true, lo=best_lo, hi=best_hi, whylo=whylo, whyhi=whyhi)
end

# ── report ─────────────────────────────────────────────────────────────────
mkpath(joinpath(HERE,"results"))
io = open(joinpath(HERE,"results","conditions.md"), "w")
emit(s="") = (println(s); println(io, s))

emit("# Conditions table and water thresholds")
emit()
emit("Relation, exactly as `stages/stage5_quality.jl:17-18` uses it:")
emit()
emit("```")
emit("sigma_bulk = sigma_w * theta_total^m + clay_term     clay_term = sigma_clay * C")
emit("sigma_w    = (sigma_bulk - clay_term) / theta_total^m")
emit("TDS        = sigma_w * 6400            [priors].mg_per_sm")
emit("```")
emit()
emit(@sprintf("Limits: theta_total ~ U(%.2f, %.2f); m ~ U(%.2f, %.2f); clay product prior %.0e-%.1f S/m.",
              PHI_LO, PHI_HI, M_LO, M_HI, CP_LO, CP_HI))
emit(@sprintf("Clay ceiling = %.1f x %.1f x (1 - theta/%.2f) S/m, so %.3f S/m at theta=0.25 and %.3f at theta=0.40.",
              SCMAX, C0, PHIREF, clay_ceiling(0.25), clay_ceiling(0.40)))
emit()
emit("Classes as the brief sets them: fresh < 1,500; brackish 1,500-3,000; saline > 3,000 mg/L.")
emit("**The code uses a single threshold of 1,000 mg/L** (`[hypothesis.salinity].value`), not these three.")
emit()

for tag in TAGS
    g = GATE[tag]; tr = truthof(tag); w = WINDOW[tag]; L = w[2]-w[1]
    Rw = thetabyens(tag); θe, θsrc = Rw.ens, Rw.src
    emit("## Column $tag")
    emit()
    emit(@sprintf("Window %.1f-%.1f m (thickness %.1f m). Gate band %.5f-%.5f S/m (conductance %.4f-%.4f S); lower edge %s, upper edge %s. Source: %s.",
                  w[1], w[2], L, g.lo, g.hi, g.lo*L, g.hi*L,
                  g.openlo ? "OPEN" : "closed", g.openhi ? "OPEN" : "closed", g.src))
    if !isempty(θe)
        p05 = minimum(quantile(v,0.05) for v in θe); p95 = maximum(quantile(v,0.95) for v in θe)
        emit(@sprintf("SNMR water content across %d ensembles: %.3f to %.3f v/v (lowest p05 to highest p95), from %s.",
                      length(θe), p05, p95, θsrc))
    else
        emit("SNMR ensemble NOT FOUND for this column.")
    end
    emit(@sprintf("Truth: sigma_bulk %.4f S/m, sigma_w %.4f S/m (TDS %.0f mg/L), clay term %.4f S/m, theta_total %.3f, theta_SMR %.3f.",
                  tr.σb, tr.σw, tr.tds, tr.clay, tr.θtot, tr.θsmr))
    emit()
    emit("### Part 1 — what each class requires")
    emit()
    emit("| class | supported | clay term required (S/m) | at central theta, m | binds low | binds high |")
    emit("|---|---|---|---|---|---|")
    for (nm, T1, T2) in CLASSES
        r  = classrange(g.lo, g.hi, PHI_LO, PHI_HI, M_LO, M_HI, T1, T2)
        θc, mc = (PHI_LO+PHI_HI)/2, (M_LO+M_HI)/2
        rc = classrange(g.lo, g.hi, θc, θc, mc, mc, T1, T2)
        emit(@sprintf("| %s | %s | %s | %s | %s | %s |", nm,
            r.supported ? "**supported**" : "**RULED OUT**",
            r.supported ? @sprintf("%.4f - %.4f", r.lo, r.hi) : "-",
            rc.supported ? @sprintf("%.4f - %.4f", rc.lo, rc.hi) : "-",
            r.whylo, r.whyhi))
    end
    emit()
    tcls = tr.tds < 1500 ? "fresh" : (tr.tds < 3000 ? "brackish" : "saline")
    rt = classrange(g.lo, g.hi, PHI_LO, PHI_HI, M_LO, M_HI,
                    tcls=="fresh" ? 0.0 : (tcls=="brackish" ? 1500.0 : 3000.0),
                    tcls=="fresh" ? 1500.0 : (tcls=="brackish" ? 3000.0 : Inf))
    inside = rt.supported && tr.clay >= rt.lo - 1e-12 && tr.clay <= rt.hi + 1e-12
    emit(@sprintf("Truth is **%s** (%.0f mg/L) with clay term %.4f S/m; that class is %s and the true clay term is %s its required range.",
                  tcls, tr.tds, tr.clay, rt.supported ? "supported" : "RULED OUT",
                  inside ? "**inside**" : "**OUTSIDE**"))
    emit()
end

# ── Part 2: adding a measurement ───────────────────────────────────────────
emit("## Part 2 — what one extra measurement would rule out (2P and 3P)")
emit()
emit("An induction log pins the window conductivity to the truth +/- 10%. A gamma log")
emit("pins the clay term to within a factor of 2 of the truth. Neither is used as an")
emit("input to anything; they are hypothetical constraints applied to the same grid.")
emit()
for tag in intersect(TAGS, ["2P","3P"])
    g = GATE[tag]; tr = truthof(tag)
    emit("### $tag")
    emit()
    emit("| scenario | fresh | brackish | saline |")
    emit("|---|---|---|---|")
    scen = [("survey alone",            g.lo, g.hi, nothing),
            ("+ induction log",         tr.σb*0.9, tr.σb*1.1, nothing),
            ("+ gamma log",             g.lo, g.hi, (tr.clay/2, tr.clay*2)),
            ("+ both",                  tr.σb*0.9, tr.σb*1.1, (tr.clay/2, tr.clay*2))]
    for (nm, slo, shi, cap) in scen
        cells = String[]
        for (cn, T1, T2) in CLASSES
            r = classrange(slo, shi, PHI_LO, PHI_HI, M_LO, M_HI, T1, T2; capclay = cap)
            push!(cells, r.supported ? @sprintf("%.4f-%.4f", r.lo, r.hi) : "**RULED OUT**")
        end
        emit("| $nm | $(cells[1]) | $(cells[2]) | $(cells[3]) |")
    end
    emit()
end

# ── Part 3: water thresholds ───────────────────────────────────────────────
emit("## Part 3 — water-content thresholds in the representative window")
emit()
emit("P(theta > x) computed separately in each ensemble, then summarised across them:")
emit("all > 0.95 -> \"greater than 0.95\"; all < 0.05 -> \"less than 0.05\"; otherwise the range.")
emit()
emit("| column | ensembles | P(theta > 0.20) | P(theta > 0.25) | P(theta > 0.30) | true theta_SMR | vintage |")
emit("|---|---|---|---|---|---|---|")
for tag in TAGS
    R = thetabyens(tag); θe, θsrc = R.ens, R.src; tr = truthof(tag)
    if isempty(θe); emit("| $tag | none found | - | - | - | - | - |"); continue; end
    cells = String[]
    for x in (0.20, 0.25, 0.30)
        ps = [count(>(x), v)/length(v) for v in θe]
        push!(cells, all(p -> p > 0.95, ps) ? "greater than 0.95" :
                     all(p -> p < 0.05, ps) ? "less than 0.05" :
                     @sprintf("%.1f - %.1f", minimum(ps), maximum(ps)))
    end
    vint = occursin("scratchpad", R.path) ? "13 Nov (200k)" : "22 Oct (100k)"
    emit(@sprintf("| %s | %d | %s | %s | %s | %.3f | %s |", tag, length(θe),
                  cells[1], cells[2], cells[3], tr.θsmr, vint))
end
emit()

# ── data for the figure ────────────────────────────────────────────────────
open(joinpath(HERE,"results","region.csv"), "w") do f
    println(f, "tag,tds,clay_lo,clay_hi,clay_ceiling_at_central_theta,truth_tds,truth_clay")
    for tag in TAGS
        g = GATE[tag]; tr = truthof(tag)
        for tds in exp.(range(log(50.0), log(60000.0), length = 160))
            r = classrange(g.lo, g.hi, PHI_LO, PHI_HI, M_LO, M_HI, tds, tds)
            r.supported || continue
            @printf(f, "%s,%.4f,%.6f,%.6f,%.6f,%.2f,%.6f\n", tag, tds, r.lo, r.hi,
                    clay_ceiling((PHI_LO+PHI_HI)/2), tr.tds, tr.clay)
        end
    end
end
println("wrote ", joinpath(HERE,"results","region.csv"))

close(io)
println("\nwrote ", joinpath(HERE,"results","conditions.md"))
