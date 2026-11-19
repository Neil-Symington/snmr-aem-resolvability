# udf_smr.jl — SNMR sounding assembly for the Upper Darling Floodplain.
#
# ek_prerun/prerun_data.jl's `build_sounding` is wired to the East Kimberley
# archive: it looks a site up in an EK spreadsheet and takes the conductive
# earth from an EK resistivity-profile directory. UDF has neither. What it does
# have is a metadata spreadsheet with the loop geometry and inclination, and
# `.mat` files in the same Vista Clara GMR format.
#
# So this reuses prerun_data.jl's *generic* machinery — pick_coil, per_q_noise
# (the sigma-hat convention), the common depth grid, sounding_hier — and
# supplies only the UDF-specific job assembly. Keeping those shared matters:
# the EK and UDF posteriors are then on the same grid, with the same noise
# model, and the two surveys can be compared without an apples-to-oranges
# caveat.
#
# The conductive earth comes from the AEM, always: there is no operator-
# supplied resistivity profile to fall back on, which is why `--supplied` is
# not offered here.

using MAT, Statistics, Printf, DelimitedFiles
using HiQGA.transD_GP, SMRPInversion, SNMRForward

# The processed-SNMR archive, from the control file rather than from one
# developer's home directory. [paths] snmr_dir is where the shipped two-site
# example lives; point it at the full archive to run more sites.
const UDF_SNMR = let c = FieldConfig.CFG[]
    c === nothing && error("udf_smr.jl needs the control file: include " *
                           "common/fieldconfig.jl and call load_config first")
    c["paths"]["snmr_dir"]
end
# The metadata CSV, likewise. The original pointed at an .xlsx beside the archive
# that a separate script converted; the shipped example carries the converted
# two-row CSV instead, which is what [paths] snmr_meta names.
const UDF_META = let c = FieldConfig.CFG[]
    c["paths"]["snmr_meta"]
end

"""
The UDF metadata row for one site.

The spreadsheet is .xlsx; `udf_meta.csv` is written beside it by
`select_udf_sites.py` so Julia does not need an Excel reader.
"""
function udf_meta(field_id)
    raw, hdr = readdlm(UDF_META, ',', String; header = true)
    ci = Dict(strip(h) => i for (i, h) in enumerate(vec(hdr)))
    i = findfirst(==(field_id), raw[:, ci["Field_ID"]])
    i === nothing && error("no metadata row for $field_id")
    (field_id = field_id,
     pulse_ms  = parse(Float64, raw[i, ci["pulse_width"]]),
     coil      = parse(Int,     raw[i, ci["detect_coil"]]),
     loop_side = parse(Float64, raw[i, ci["coil_diameter"]]),
     inc_deg   = parse(Float64, raw[i, ci["bo_inc_angle_deg"]]),
     easting   = parse(Float64, raw[i, ci["easting"]]),
     northing  = parse(Float64, raw[i, ci["northing"]]))
end

"""
Choose the processed `.mat` for a site.

Vista Clara writes several: the raw sounding, `_noemi` where EMI removal ran,
and `_harmonic_removed` where harmonics were stripped. The operator's preferred
product is the most-processed one, so the preference order is
harmonic+noemi > noemi > harmonic > raw, and which one was taken is recorded.
"""
function udf_matfile(field_id)
    # TWO LAYOUTS. The GA delivery nests each site's products under
    # <site>/Processed/, and that is what this function expected. The two-site
    # example shipped in fieldrun/example/snmr/ is FLAT -- SMR_UD21_noemi.mat
    # beside SMR_UD24_harmonic_removed_noemi.mat -- because only the preferred
    # product per site ships. Expecting only the archive layout meant the SNMR
    # stage could not read the data the repository ships with it: it failed with
    # "no Processed directory for SMR_UD21" on both sites.
    nested = joinpath(UDF_SNMR, field_id, "Processed")
    dir, cands = if isdir(nested)
        nested, filter(x -> endswith(x, ".mat"), readdir(nested))
    elseif isdir(UDF_SNMR)
        # flat: the site id is the filename prefix
        UDF_SNMR, filter(x -> endswith(x, ".mat") && startswith(x, field_id),
                         readdir(UDF_SNMR))
    else
        error("""
    no SNMR directory at $UDF_SNMR

      [paths] snmr_dir must name either the GA archive, which nests each site
      under <site>/Processed/, or a flat directory of per-site .mat files as
      fieldrun/example/snmr/ does.
    """)
    end
    isempty(cands) && error("""
    no .mat file for $field_id in $dir

      Looked for $(dir == nested ? "any .mat" : "a .mat beginning \"$field_id\"").
    """)
    # Vista Clara writes several products; the operator's preferred one is the
    # most processed, so this is a preference order and not a filename.
    for pat in ("harmonic_removed_noemi", "noemi", "harmonic_removed", "")
        i = findfirst(x -> occursin(pat, x), cands)
        i !== nothing && return joinpath(dir, cands[i]), cands[i]
    end
    error("no .mat for $field_id")
end

"""
Read and prepare one UDF sounding, everything except the kernel.

Mirrors prerun_data.jl `build_sounding` field for field, minus the supplied
conductivity. `F` and `cond` are filled in afterwards by the AEM draw.
"""
function build_udf_sounding(job)
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

    _, znall, zbnd = CFG.grid()
    (ok = true, zbnd = zbnd, znall = znall, q = q, V0 = V0, ϕ = ϕ,
     σ̂ = σ̂, σ̂_branch = String(CFG.SIGMA_HAT_BRANCH), why = "",
     Be = Be, freq = f, coil_used = coil, coil_swapped = swapped,
     tail_nV = median(tail) * 1e9, nq = length(V0), ref_loop = "")
end

"""
Build the kernel on one Stage-1 AEM conductivity draw.

`ens` selects which draw, so the pooled posterior marginalises over NENS
kernel realisations rather than conditioning on one — the same coarse version
of the paper's Stage 2 used at East Kimberley, and labelled the same way
through `n_kernel_draws`.

This matters far more here than it did there. The UDF trans-D AEM determines
conductance to about +/-15% but POINTWISE conductivity at 20-80 m only to a
factor of 4-64, because the response saturates over a thick conductor. The
kernel is a function of sigma(z), so kernel uncertainty is the dominant term
rather than a refinement.
"""
function udf_kernel(b, job, field_id, OUT; ens::Int = 0, occam::Bool = false)
    td = joinpath(OUT, "aem_td_$(field_id).jld2")
    oc = joinpath(OUT, "aem_occam_$(field_id).jld2")
    local σ, ztop, src, nk
    if isfile(td) && !occam
        A = load(td); D = A["sigma_draws"]; nk = size(D, 2)
        j = mod1(max(ens, 1), nk)
        σ = D[:, j]; ztop = A["zbnd"]; src = "trans-D AEM draw $j of $nk"
    else
        isfile(oc) || error("no AEM inversion for $field_id")
        A = load(oc); σ = A["sigma"]; ztop = A["zbnd"]; nk = 1
        src = "AEM Occam smooth model"
    end
    thick = diff(ztop)
    @assert length(thick) == length(σ) - 1
    cm = SNMRForward.ConductivityModel(σ, thick)
    F = SMRPInversion.MRSForward_square(job.loop_side / 2, b.zbnd, b.q,
            deg2rad(job.inc_deg), 0.0, b.Be, cm)
    merge(b, (F = F, cond = cm, cond_source = src,
              aem_sigma = σ, aem_ztop = ztop, n_kernel_draws = nk))
end
