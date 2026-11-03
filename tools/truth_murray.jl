# truth_murray.jl — the one synthetic truth for the hypothesis paper (§1).
#
# MURRAY BASIN, RIVERINE PLAIN (e.g. Deniliquin–Shepparton). Supersedes the
# East Kimberley column in synthetics/superseded/, at the user's direction
# (12 Oct 2026). One truth, one site.
#
#   Shepparton Fm  clay-rich floodplain with palaeochannel sands — CONDUCTIVE cap
#   Calivil Fm     fluvial sand/gravel, the main irrigation aquifer — RESISTIVE
#   Renmark Group  sand/clay/lignite, more evolved water
#   weathered then unweathered Palaeozoic basement
#
# This is the classic Riverine Plain AEM signature: a resistive fresh aquifer
# sandwiched between two conductors. It is a harder and more interesting target
# than a simple layer cake, because the conductive cap is conductive for TWO
# reasons — clay and salt — and σ_bulk alone cannot separate them. That
# ambiguity is the subject of S-F and of Figure 3 panel f.
#
# Transform, unchanged from petrophysical_scenario.jl so the synthetic stays
# comparable with the existing scenarios:
#   θ_total = φ·S
#   θ_SMR   = θ_total·(1 − clay)^1.5          P_CLAY   = 1.5
#   σ_bulk  = σ_w·θ_total^1.7 + 0.5·clay      M_ARCHIE = 1.7, σ_SURF = 0.5 S/m
#
# SALINITY ARCHITECTURE: fresh shallow over saline deep, as chosen. Note this is
# the Mallee / western-Murray sense rather than the eastern Riverine Plain's
# shallow salinisation. Basement pore water is the MOST saline in the column,
# which is the normal evolution sequence (the EK column had this backwards).
#
# CALIBRATION: the aquifer is a deliberate BEST CASE, not a typical section —
# see `report`. Stated so no reader takes the resulting probabilities as typical
# Murray Basin performance.
#
#   julia --project=. synthetics/truth_murray.jl

using Random, Statistics, Printf, LinearAlgebra, JLD2

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(HERE, "..", "data")

include(joinpath(ROOT, "ek_prerun", "prerun_config.jl")); const CFG = PrerunConfig

# ── transform constants ──────────────────────────────────────────────────────
const P_CLAY     = 1.5
const M_ARCHIE   = 1.7
const σ_SURF     = 0.5      # S/m
const LOG_σ_φ    = 0.10
const LOG_σ_CLAY = 0.15

# ── AGRF, Riverine Plain (~35.5°S, 145°E) ───────────────────────────────────
# NOT a metadata edit: the SMR kernel is strongly inclination-dependent, and
# this is a 22.5° change from the East Kimberley value. Larmor = 42.5774 Hz/μT.
const BE_T       = 5.80e-5        # T  (58,000 nT)
const INC_DEG    = -66.5
const LARMOR_HZ  = 0.0425774 * 5.80e4
const LOOP_SIDE  = 88.6           # m, SMRPI convention

# ── grid ─────────────────────────────────────────────────────────────────────
const DZ_FINE = 0.25
const ZMAX    = 200.0
const ZF      = collect(DZ_FINE/2 : DZ_FINE : ZMAX)

struct Unit
    top::Float64; bot::Float64; name::String
    φ::Float64; S::Float64; clay::Float64; σw::Float64
end

const UNITS = [
    Unit(  0.0,   6.0, "Shepparton Fm, unsaturated", 0.33, 1.00, 0.40, 0.05),
    Unit(  6.0,  25.0, "Shepparton Fm, clay/sand",   0.35, 1.00, 0.45, 0.05),
    Unit( 25.0,  45.0, "Calivil Fm, sand/gravel",    0.30, 1.00, 0.08, 0.08),
    Unit( 45.0,  62.0, "Renmark Gp, upper sand",     0.28, 1.00, 0.25, 0.50),
    Unit( 62.0,  75.0, "Renmark Gp, lignitic",       0.30, 1.00, 0.35, 0.80),
    Unit( 75.0,  92.0, "weathered Palaeozoic",       0.12, 1.00, 0.30, 1.20),
    Unit( 92.0,  98.0, "saprock",                    0.06, 1.00, 0.12, 1.40),
    Unit( 98.0, ZMAX,  "unweathered Palaeozoic",     0.02, 1.00, 0.03, 1.50),
]

# Riverine Plain water tables are shallow. Clay-rich material also has a genuine
# multi-metre capillary rise, so unlike the EK column this label is honest.
const CAP_FRINGE_TOP = 0.15
const CAP_FRINGE_Z   = 4.0

function unit_field(z, trans, get)
    v = similar(z)
    @inbounds for (i, zi) in enumerate(z)
        k = findlast(u -> zi >= u.top, UNITS); k = k === nothing ? 1 : k
        val = get(UNITS[k])
        if k < length(UNITS)
            b = UNITS[k+1].top; h = trans[k] / 2
            if zi > b - h
                f = clamp((zi - (b - h)) / trans[k], 0.0, 1.0)
                val = (1 - f) * val + f * get(UNITS[k+1])
            end
        end
        if k > 1
            b = UNITS[k].top; h = trans[k-1] / 2
            if zi < b + h
                f = clamp((zi - (b - h)) / trans[k-1], 0.0, 1.0)
                val = (1 - f) * get(UNITS[k-1]) + f * get(UNITS[k])
            end
        end
        v[i] = val
    end
    v
end

sat_field(z) = [zi >= CAP_FRINGE_Z ? 1.0 :
                CAP_FRINGE_TOP + (1.0 - CAP_FRINGE_TOP) * zi / CAP_FRINGE_Z for zi in z]

function gp_draw(z; corr = 6.0, seed = 1)
    n = length(z); Random.seed!(seed)
    K = [exp(-0.5 * ((z[i] - z[j]) / corr)^2) for i in 1:n, j in 1:n] + 1e-8I
    cholesky(Symmetric(K)).L * randn(n)
end

function build(; seed = 20261012)
    Random.seed!(seed)
    trans = 3.0 .+ 3.0 * rand(length(UNITS) - 1)
    φ0   = unit_field(ZF, trans, u -> u.φ)
    clay0= unit_field(ZF, trans, u -> u.clay)
    σw   = unit_field(ZF, trans, u -> u.σw)
    S    = sat_field(ZF)
    φ    = clamp.(φ0   .* exp.(gp_draw(ZF; corr = 6.0, seed = seed + 1) * LOG_σ_φ),    0.005, 0.5)
    clay = clamp.(clay0.* exp.(gp_draw(ZF; corr = 6.0, seed = seed + 2) * LOG_σ_CLAY), 0.0,   1.0)
    θ_total = φ .* S
    θ_SMR   = θ_total .* clamp.((1 .- clay) .^ P_CLAY, 0.0, 1.0)
    σ_bulk  = σw .* (θ_total .^ M_ARCHIE) .+ σ_SURF .* clay
    (z = ZF, dz = DZ_FINE, trans = trans, φ = φ, S = S, clay = clay, σw = σw,
     θ_total = θ_total, θ_SMR = θ_SMR, σ_bulk = σ_bulk, seed = seed)
end

# ── hypothesis truths (§5) ───────────────────────────────────────────────────
colint(x, z, dz, z0, z1) = sum(x[i] * dz for i in eachindex(z) if z0 <= z[i] <= z1)

function crossing(x, z, cut)
    for i in 2:length(z)
        if (x[i-1] - cut) * (x[i] - cut) <= 0 && x[i] != x[i-1]
            return z[i-1] + (cut - x[i-1]) / (x[i] - x[i-1]) * (z[i] - z[i-1])
        end
    end
    NaN
end

function column_quantile(x, z, dz, zt, frac)
    k = findall(zi -> zi <= zt, z); c = cumsum(x[k] .* dz)
    j = findfirst(>=(frac * c[end]), c); j === nothing ? NaN : z[k[j]]
end

function maxrun(z, dz, zt, tests)
    run = 0.0; best = 0.0
    for i in eachindex(z)
        if z[i] <= zt && all(t -> t(i), tests)
            run += dz; best = max(best, run)
        else
            run = 0.0
        end
    end
    best
end

const σW_FRESH  = 0.15   # S/m = 1,500 μS/cm
const σW_SALINE = 0.30   # S/m = 3,000 μS/cm

# S-F's window moved from 30–38 m (EK) to 10–20 m. In the Murray column the
# clay-bound, lower-θ, conductive zone is the Shepparton Fm cap, not a mid-depth
# aquitard — and distinguishing clay from salt THERE is the real Riverine Plain
# AEM problem. Recorded as a dated change to the pre-registered window rather
# than an edit, per the §3 rule.
const SF_WINDOW = (10.0, 20.0)
const SF_WINDOW_NOTE = "moved from 30–38 m on 2026-10-12 with the basin change"

function truths(T)
    z, dz = T.z, T.dz
    fresh = [T.σw[i] < σW_FRESH for i in eachindex(z)]
    k = findall(zi -> SF_WINDOW[1] <= zi <= SF_WINDOW[2], z)
    (SA = colint(T.θ_SMR, z, dz, 0.0, 60.0),
     SB = column_quantile(T.θ_SMR, z, dz, 60.0, 0.5),
     SC = colint(T.θ_SMR .* fresh, z, dz, 0.0, 60.0),
     SD = crossing(T.σw, z, σW_FRESH),
     SE = maxrun(z, dz, 60.0, [i -> T.θ_SMR[i] > 0.15,
                               i -> T.σw[i] < σW_FRESH,
                               i -> T.clay[i] < 0.2]),
     SF_clay = mean(T.clay[k]), SF_σw = mean(T.σw[k]),
     SG = UNITS[7].top)
end

function report(T, tr)
    @printf("truth_murray — Murray Basin, Riverine Plain\n")
    @printf("seed %d, fine grid %.2f m to %.0f m (%d cells)\n", T.seed, T.dz, ZMAX, length(T.z))
    @printf("AGRF: Be = %.0f nT, inclination = %.1f°, Larmor = %.0f Hz; loop %.1f m square\n\n",
            BE_T*1e9, INC_DEG, LARMOR_HZ, LOOP_SIDE)
    println("transition thicknesses (m), drawn U(3,6) per boundary:")
    for (i, t) in enumerate(T.trans)
        @printf("  %5.1f m  %-28s | %-28s  %.2f m\n",
                UNITS[i+1].top, UNITS[i].name, UNITS[i+1].name, t)
    end
    # Trim 3 m off each end to measure clear of the transitions, but fall back
    # to the whole unit where it is too thin for that — otherwise thin units
    # return NaN and silently fail the monotonicity check below.
    um(x, u) = begin
        k = findall(zi -> u.top + 3 <= zi <= u.bot - 3, T.z)
        isempty(k) && (k = findall(zi -> u.top <= zi <= u.bot, T.z))
        isempty(k) ? NaN : mean(x[k])
    end
    println("\nunit means, measured away from the transitions:")
    @printf("  %-28s %7s %7s %7s %8s %9s\n", "unit", "phi", "clay", "th_SMR", "sig_w", "sig_bulk")
    for u in UNITS
        u.bot - u.top < 7 && continue
        @printf("  %-28s %7.3f %7.3f %7.3f %8.3f %9.4f\n", u.name,
                um(T.φ, u), um(T.clay, u), um(T.θ_SMR, u), um(T.σw, u), um(T.σ_bulk, u))
    end

    cap = um(T.σ_bulk, UNITS[2]); aq = um(T.σ_bulk, UNITS[3]); rn = um(T.σ_bulk, UNITS[4])
    println("\n§1 checks — the Riverine Plain signature is a RESISTOR BETWEEN CONDUCTORS:")
    @printf("  Shepparton cap σ_bulk        : %.4f S/m\n", cap)
    @printf("  Calivil aquifer σ_bulk       : %.4f S/m  (%.1f× more resistive than cap)\n", aq, cap/aq)
    @printf("  Renmark σ_bulk               : %.4f S/m  (%.1f× more conductive than aquifer)\n", rn, rn/aq)
    @printf("  aquifer visible to AEM?      : %s\n", (cap/aq > 2 && rn/aq > 2) ? "YES, bracketed both sides" : "marginal")
    @printf("  σ_w increases monotonically with depth? %s\n",
            all(i -> um(T.σw, UNITS[i]) <= um(T.σw, UNITS[i+1]) + 1e-9, 1:length(UNITS)-1) ?
            "yes — normal evolution sequence" : "NO — check, basement should be most evolved")

    println("\nCALIBRATION — this is a labelled BEST CASE, not a typical section:")
    @printf("  Calivil θ_SMR = %.3f. Realistic Riverine Plain SMR free water in the\n", um(T.θ_SMR, UNITS[3]))
    println("  Calivil is nearer 0.15–0.25; this sits at or above the top of that range.")
    println("  Chosen so C1–C4 separate visibly. Do not read the probabilities as")
    println("  typical Murray Basin performance.")

    println("\nhypothesis truths (§5):")
    @printf("  S-A  free-water column 0–60 m            = %7.3f m   > 6 m ? %s\n", tr.SA, tr.SA > 6 ? "TRUE" : "FALSE")
    @printf("  S-B  depth to 50%% of the 0–60 m column   = %7.2f m\n", tr.SB)
    @printf("  S-C  fresh (σ_w<%.2f) column 0–60 m      = %7.3f m   > 4 m ? %s\n", σW_FRESH, tr.SC, tr.SC > 4 ? "TRUE" : "FALSE")
    @printf("  S-D  depth σ_w crosses %.2f S/m          = %7.2f m\n", σW_FRESH, tr.SD)
    @printf("  S-E  producible fresh interval above 60 m= %7.2f m   ≥ 5 m ? %s\n", tr.SE, tr.SE >= 5 ? "TRUE" : "FALSE")
    @printf("  S-F  %.0f–%.0f m: clay = %.3f, σ_w = %.3f S/m → %s\n", SF_WINDOW...,
            tr.SF_clay, tr.SF_σw,
            tr.SF_clay > 0.4 && tr.SF_σw < σW_SALINE ? "CLAY-BOUND (true)" : "ambiguous")
    @printf("       (%s)\n", SF_WINDOW_NOTE)
    @printf("  S-G  weathered/fresh basement contact    = %7.1f m   < 100 m ? %s\n", tr.SG, tr.SG < 100 ? "TRUE" : "FALSE")
end

function main()
    mkpath(OUT)
    T = build(); tr = truths(T); report(T, tr)
    zsmr = CFG.grid()[3]
    smr_at(x) = [mean(x[findall(zi -> abs(zi - zs) <= 1.5, T.z)]) for zs in zsmr]
    jldsave(joinpath(OUT, "truth_murray.jld2");
        z = T.z, dz = T.dz, trans = T.trans, seed = T.seed,
        φ = T.φ, S = T.S, clay = T.clay, σw = T.σw,
        θ_total = T.θ_total, θ_SMR = T.θ_SMR, σ_bulk = T.σ_bulk,
        zsmr = zsmr, θ_smrgrid = smr_at(T.θ_SMR), σ_smrgrid = smr_at(T.σ_bulk),
        truths = tr, units = [(u.top, u.bot, u.name) for u in UNITS],
        P_CLAY = P_CLAY, M_ARCHIE = M_ARCHIE, σ_SURF = σ_SURF,
        σW_FRESH = σW_FRESH, σW_SALINE = σW_SALINE,
        BE_T = BE_T, INC_DEG = INC_DEG, LARMOR_HZ = LARMOR_HZ, LOOP_SIDE = LOOP_SIDE,
        SF_WINDOW = SF_WINDOW, basin = "Murray Basin, Riverine Plain")
    println("\nwrote ", joinpath(OUT, "truth_murray.jld2"))
end

main()
