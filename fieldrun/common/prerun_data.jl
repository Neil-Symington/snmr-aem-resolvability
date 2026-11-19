# Ported from udf_realdata/ into the repository. The include paths and the
# project are the only things changed from the working version: COMMON is
# fieldrun/common/ and the project is the repository root, not ../ek_prerun.
# prerun_data.jl — sounding assembly for the EFTF survey pre-run.
#
# Differs from ek_inversions/ek_common.jl in three ways that matter:
#   1. ONE common depth grid for every sounding (spec §6 wants posteriors on a
#      common grid), not the per-site z_max grid the Occam production run used.
#   2. σ̂ is per pulse moment, recovered from the FID residual, rather than a
#      single median σ (spec §2 shape × level; RUNBOOK §4 open item 2).
#   3. Bounds are [0.001, 0.45] per spec §1, not [1e-4, 0.5].
#
# RUNBOOK §1.1: every non-trivial loop lives in a function. Top-level soft scope
# silently discards assignments in a non-interactive session.

using MAT, Statistics, Printf
using HiQGA.transD_GP, SMRPInversion, SNMRForward
import SMRPInversion.ProcessingTools: analytic_signal, est_E0_T2

include(joinpath(@__DIR__, "..", "..", "src", "prerun_config.jl"))
using .PrerunConfig
const CFG = PrerunConfig

# THE EAST KIMBERLEY PATHS, which this file is wired to and the UDF workflow is
# not. They were four `const`s holding one developer's absolute paths, evaluated
# at LOAD TIME -- so a clone carried four dead paths that no error mentioned,
# because they are only used lazily. That was named as a porting blocker three
# rounds before this file moved.
#
# They are now resolved from the control file's optional [paths] ek_root, and
# absent that they are empty strings that say what is wrong when something tries
# to use them. Nothing in the UDF stages reads them: this file is included for
# its GENERIC machinery -- pick_coil, per_q_noise, the sigma-hat convention, the
# common depth grid, sounding_hier -- which is shared with East Kimberley so the
# two surveys' posteriors sit on one grid.
const EK_ROOT = let c = FieldConfig.CFG[]
    c === nothing ? "" : String(get(get(c, "paths", Dict()), "ek_root", ""))
end
ek_path(parts...) = isempty(EK_ROOT) ? error("""
    This is an East Kimberley path and [paths] ek_root is not set.

      Wanted: $(join(parts, "/"))

      The UDF stages do not need it -- they use this file's generic machinery
      only. If you are running East Kimberley, set ek_root in the control file.
    """) : joinpath(EK_ROOT, parts...)

SHEET()       = ek_path("EK_conductive_earth_inversion_spreadsheet.csv")
DATASHEET()   = ek_path("EK_SNMR_data_spreadsheet.csv")
PROCESSED()   = ek_path("Processed_data")
RESPROF_DIR() = ek_path("resistivity_profiles-20260925T031913Z-1-001",
                        "resistivity_profiles")

# ─────────────────────────────────────────────────────────────────────────────
# Spreadsheets
# ─────────────────────────────────────────────────────────────────────────────

"FID rows of the inversion spreadsheet."
function read_sheet(; sheet = SHEET)
    raw = readlines(sheet)
    col = Dict(strip(h) => i for (i, h) in enumerate(split(raw[1], ',')))
    rows = NamedTuple[]
    for line in raw[2:end]
        isempty(strip(line)) && continue
        f = split(line, ',')
        f[col["pulse_sequence"]] == "FID" || continue
        push!(rows, (
            field_id  = String(f[col["Field_ID"]]),
            filename  = String(f[col["filename"]]),
            pulse_ms  = parse(Float64, f[col["pulse_length"]]),
            coil      = parse(Int,     f[col["detect_coil"]]),
            loop_side = parse(Float64, f[col["coil_diameter"]]),
            inc_deg   = parse(Float64, f[col["bo_inc_angle_deg"]]),
            z_max     = parse(Float64, f[col["z_max"]]),
            condname  = String(f[col["loadcondname"]]),
        ))
    end
    rows
end

"""
Reference loop per (site, pulse length) from EK_SNMR_data_spreadsheet.csv.

This is the one §6 condition covariate that EK actually carries; TMI gradient,
distance to infrastructure, month and stack count do not exist in the archive.
"""
function read_refloops(; sheet = DATASHEET)
    isfile(sheet) || return Dict{String,String}()
    raw = readlines(sheet)
    col = Dict(strip(h) => i for (i, h) in enumerate(split(raw[1], ',')))
    haskey(col, "Ref_loop") || return Dict{String,String}()
    d = Dict{String,String}()
    for line in raw[2:end]
        isempty(strip(line)) && continue
        f = split(line, ',')
        length(f) < maximum(values(col)) && continue
        key = "$(strip(f[col["Field_ID"]]))__$(strip(f[col["Pulse_length_(ms)"]]))"
        d[key] = String(strip(f[col["Ref_loop"]]))
    end
    d
end

"""
One job per sounding: the EMI-removed file where the operator made one, else raw.

The Occam production run inverted both members of each pair; the pre-run wants
one row per sounding (spec §6), so the preferred file is taken. `emi` records
which, and `haspair` whether a counterpart exists.
"""
function build_jobs(rows = read_sheet())
    jobs = NamedTuple[]
    for r in rows
        base = replace(r.filename, "_noemi.mat" => ".mat")
        noemi = replace(base, ".mat" => "_noemi.mat")
        praw  = joinpath(PROCESSED(), r.field_id, base)
        pnoe  = joinpath(PROCESSED(), r.field_id, noemi)
        hasraw, hasnoe = isfile(praw), isfile(pnoe)
        (hasraw || hasnoe) || continue
        fn, tag, path = hasnoe ? (noemi, :noemi, pnoe) : (base, :raw, praw)
        push!(jobs, merge(r, (
            filename = fn, emi = tag, haspair = (hasraw && hasnoe),
            matpath  = path,
            sounding_id = "$(r.field_id)__$(replace(base, ".mat" => ""))")))
    end
    unique(j -> j.sounding_id, jobs)
end

# ─────────────────────────────────────────────────────────────────────────────
# Conductivity
# ─────────────────────────────────────────────────────────────────────────────

"Load a site resistivity profile as (ConductivityModel, ρ, thickness)."
function load_conductivity(condname)
    stem  = replace(condname, "_res_profile.txt" => "")
    cands = [condname, replace(stem, "_" => "-") * "_res_profile.txt"]
    idx   = findfirst(c -> isfile(joinpath(RESPROF_DIR(), c)), cands)
    idx === nothing && return nothing
    ρ = Float64[]; t = Float64[]
    for line in Iterators.drop(eachline(joinpath(RESPROF_DIR(), cands[idx])), 1)
        f = split(strip(line)); isempty(f) && continue
        push!(ρ, parse(Float64, f[1]))
        length(f) >= 2 && push!(t, parse(Float64, f[2]))
    end
    isempty(ρ) && return nothing
    # ρ = 0 is a missing-value placeholder: σ = Inf NaNs the kernel. Carry the
    # last valid value down. Affects 3 of 266 EK profiles.
    repaired = false
    for i in eachindex(ρ)
        if !(ρ[i] > 0)
            prev = findlast(j -> ρ[j] > 0, 1:i-1)
            prev === nothing && return nothing
            ρ[i] = ρ[prev]; repaired = true
        end
    end
    length(t) >= length(ρ) && (t = t[1:length(ρ)-1])
    length(t) == length(ρ) - 1 || return nothing
    (model = SNMRForward.ConductivityModel(1.0 ./ ρ, t),
     ρ = ρ, thick = t, repaired = repaired)
end

"Cumulative conductance ∫σ dz to `zmax` (siemens) — spec §0 stratifier, §6 covariate."
function conductance(ρ, thick, zmax)
    tops = vcat(0.0, cumsum(thick)); S = 0.0
    for i in eachindex(ρ)
        tops[i] >= zmax && break
        bot = i <= length(thick) ? min(tops[i+1], zmax) : zmax
        S += (bot - tops[i]) / ρ[i]
    end
    S
end

"Spec §0 conductivity class from conductance to 50 m."
cond_class(S50) = S50 < 1 ? "<1S" : (S50 <= 10 ? "1-10S" : ">10S")

# ─────────────────────────────────────────────────────────────────────────────
# Per-pulse-moment noise — the σ̂ⱼ shape
# ─────────────────────────────────────────────────────────────────────────────

"""
Per-pulse-moment σ̂ on V₀, recovered from the FID residual.

`est_E0_T2` fits log|analytic signal| against time and propagates the envelope
residual through to σ(E₀). That per-q σ(E₀) is the noise on the very quantity
being inverted, so it is σ̂ⱼ. Its *level* is deliberately not corrected here —
the hierarchical `s` is what estimates the level, and on real data `s` carries
noise plus representation error (RUNBOOK §4, item 3).

`analytic_signal` is a true Hilbert transform (zero the negative frequencies,
double the positive ones), so |z| = A with no factor-of-2 loss — verified
against a synthetic: E₀ = 1499.1 nV recovered from a 1500 nV truth. The ×2 trap
in RUNBOOK §1.7 applies to multiply-and-lowpass demodulation of raw .lvm, not
to this path.

Also returns the record tail noise (mean envelope over the last quarter of each
record, where the NMR signal has decayed) — spec §6's "record noise level (tail)".
"""
function per_q_noise(t, fid)
    nq = size(fid, 2)
    i0 = max(1, round(Int, 0.75 * size(fid, 1)))
    σ̂ = fill(NaN, nq); tailnV = fill(NaN, nq)
    for i in 1:nq
        col = @view fid[:, i]
        try
            σ̂[i] = est_E0_T2(t, collect(col))[3]
            tailnV[i] = mean(abs.(analytic_signal(collect(col)))[i0:end])
        catch
            # leave NaN; caller decides
        end
    end
    (σ̂ = σ̂, tail = tailnV)
end

# ─────────────────────────────────────────────────────────────────────────────
# Sounding assembly
# ─────────────────────────────────────────────────────────────────────────────

"Pick the data channel. detect_coil is authoritative when it holds float data."
function pick_coil(d, spec)
    hasdata(i) = (k = "coil_$(i)_fid";
                  haskey(d, k) && eltype(d[k]) <: AbstractFloat && !all(iszero, d[k]))
    hasdata(spec) && return spec, false
    avail = [i for i in 1:8 if hasdata(i)]
    length(avail) == 1 ? (avail[1], true) : (0, false)
end

"""
Read and prepare one sounding on the common grid.

Returns `(ok = false, why = ...)` when it cannot be inverted. On success every
field needed by both the trans-D and the Occam companion is present, plus the
§6 condition covariates that EK can supply.
"""
function build_sounding(job; refloops = Dict{String,String}())
    d = matread(job.matpath)
    coil, swapped = pick_coil(d, job.coil)
    coil == 0 && return (ok = false, why = "no channel holds data")
    fid = Float64.(d["coil_$(coil)_fid"])
    all(isfinite, fid) || return (ok = false, why = "non-finite samples")

    t  = vec(Float64.(d["time_fid"]))
    q  = vec(Float64.(d["pulse_moment"]))
    f  = Float64(d["detect_frequency"][1])
    Be = 2π * f / SMRPInversion.γh

    V0e, ϕe = SMRPInversion.get_sounding_curve(t, fid)
    V0 = vec(V0e); ϕ = vec(ϕe)
    pq = per_q_noise(t, fid)

    keep = isfinite.(V0) .& isfinite.(ϕ) .& (V0 .> 0) .&
           isfinite.(pq.σ̂) .& (pq.σ̂ .> 0)
    sum(keep) >= 6 || return (ok = false, why = "only $(sum(keep)) usable pulse moments")
    q, V0, ϕ = q[keep], V0[keep], ϕ[keep]
    σ̂, tail  = pq.σ̂[keep], pq.tail[keep]

    cond = load_conductivity(job.condname)
    cond === nothing && return (ok = false, why = "no usable resistivity profile")

    _, znall, zbnd = CFG.grid()
    F = SMRPInversion.MRSForward_square(job.loop_side / 2, zbnd, q,
            deg2rad(job.inc_deg), 0.0, Be, cond.model)

    rl = get(refloops, "$(job.field_id)__$(Int(job.pulse_ms))", "")

    (ok = true, zbnd = zbnd, znall = znall, q = q, V0 = V0, ϕ = ϕ, F = F,
     σ̂ = σ̂, σ̂_branch = String(CFG.SIGMA_HAT_BRANCH), why = "",
     Be = Be, freq = f, cond = cond.model,
     S50 = conductance(cond.ρ, cond.thick, 50.0),
     S100 = conductance(cond.ρ, cond.thick, 100.0),
     res_repaired = cond.repaired, coil_used = coil, coil_swapped = swapped,
     tail_nV = median(tail) * 1e9, nq = length(V0), ref_loop = rl)
end

"Hierarchical sounding: σᵢ² = s·σ̂ᵢ², s sampled with a log-normal prior (spec §2)."
function sounding_hier(b)
    SMRPInversion.newSMRSoundingHierarchical(b.V0, b.ϕ, b.F;
        σ̂_V0 = b.σ̂,
        s_prior_center   = CFG.S_PRIOR_CENTER,
        s_prior_logwidth = CFG.S_PRIOR_TAU,
        s_init = CFG.S_INIT,
        linearsat = true, amponly = true)
end

"Fixed-σ sounding at σⱼ² = s·σ̂ⱼ² — the Occam companion's noise (spec §5)."
function sounding_fixed(b; s = CFG.S_PRIOR_CENTER)
    σ = sqrt(s) .* b.σ̂
    SMRPInversion.newSMRSounding(b.V0, b.ϕ, b.F;
        σ_V0 = σ, σ_ϕ = σ ./ b.V0,
        linearsat = true, amponly = true, mult = false, showplot = false)
end

"Deterministic per-sounding RNG stream (spec §3: seeds logged, fixed per site ID)."
seed_for(sounding_id, ensemble) =
    CFG.SEED_BASE + (hash(sounding_id) % 100_000) * 10 + ensemble
