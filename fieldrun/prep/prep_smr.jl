# prep_smr.jl -- turn a GMR field delivery into the fieldrun SNMR inputs.
#
#   julia --project=. fieldrun/prep/prep_smr.jl --src <delivery dir> --out <prepared dir>
#         [--sites Mulan31,Mulan32] [--prefix ""] [--year 2026]
#
# A SITE is found from its *_params.txt (antenna corners and side length, from
# the site-layout tool). Two layouts are accepted:
#
#   DELIVERY (the full dataset):   <site>/SiteInfo/<site>_params.txt
#                                  <site>/SiteInfo/<site>_agrf25.txt
#                                  <site>/Reprocessed_<mm_yy>/QC_FID <date>/QC_FID <date>[_noemi].mat
#                                  (+ photos, shapefiles, FIG/PNG we do not read)
#     The site root is the parent of SiteInfo; the agrf25 must sit beside the
#     params file and the .mat anywhere under the site root OUTSIDE SiteInfo.
#
#   FLAT (the first staging copy): <site>/{_params.txt, _agrf25.txt, .mat}
#
# Several .mat per site is normal (raw and _noemi, repeat reprocessings). The
# choice is made the way udf_smr.jl's udf_matfile makes it: under the newest
# Reprocessed_* folder if there are several, then by the processing preference
# harmonic_removed_noemi > noemi > harmonic_removed > raw. What was chosen and
# what was passed over are both printed.
#
# Writes, under --out:
#   snmr_metadata.csv          one row per site, the columns udf_smr.jl reads
#                              (Field_ID ... pulse_width) plus provenance
#   snmr/<Field_ID>_<mat>.mat  the .mat, copied (not converted) so udf_matfile
#                              finds it in the flat layout
#
# Point [paths] snmr_meta and snmr_dir at these. sites.csv (the SMR <-> AEM
# pairing) is NOT written here: it needs the AEM delivery.

using MAT, Printf, Dates

function parse_args(a)
    d = Dict{String,String}()
    i = 1
    while i <= length(a)
        startswith(a[i], "--") || error("unexpected argument $(a[i])")
        i < length(a) || error("$(a[i]) needs a value")
        d[a[i][3:end]] = a[i+1]; i += 2
    end
    haskey(d, "src") && haskey(d, "out") || error("usage: --src <dir> --out <dir> [--sites A,B] [--prefix P]")
    d
end

allfiles(dir) = [joinpath(r, f) for (r, _, fs) in walkdir(dir) for f in fs if !startswith(f, ".")]

function one_match(files, re, what, site)
    m = filter(f -> occursin(re, basename(f)), files)
    isempty(m) && error("$site: no $what")
    length(m) > 1 && error("$site: $(length(m)) candidates for $what, expected one:\n  " * join(m, "\n  "))
    only(m)
end

"""
The processed .mat for a site: newest Reprocessed_* folder first, then the
processing preference. Returns (chosen, passed_over).
"""
function pick_mat(mats, root, site)
    isempty(mats) && error("$site: no .mat under $root")
    # newest reprocessing: Reprocessed_<mm>_<yy>, compared as (yy, mm)
    stamp(f) = (m = match(r"Reprocessed_(\d{1,2})_(\d{2})", relpath(f, root));
                m === nothing ? (-1, -1) : (parse(Int, m[2]), parse(Int, m[1])))
    newest = maximum(stamp.(mats))
    pool = filter(f -> stamp(f) == newest, mats)
    for pat in ("harmonic_removed_noemi", "noemi", "harmonic_removed", "")
        c = filter(f -> occursin(pat, basename(f)), pool)
        isempty(c) && continue
        length(c) > 1 && pat != "" && error("$site: $(length(c)) candidates for \"$pat\":\n  " * join(c, "\n  "))
        # the raw ("") tier can hold several only if no processed one exists
        length(c) > 1 && error("$site: $(length(c)) .mat files and none marked as processed:\n  " * join(c, "\n  "))
        return only(c), setdiff(mats, c)
    end
end

"Antenna geometry from *_params.txt."
function read_params(f)
    L = readlines(f)
    val(key) = (i = findfirst(l -> startswith(strip(l), key), L);
                i === nothing && error("$f: no '$key' line"); strip(split(L[i], ':'; limit = 2)[2]))
    proj  = val("Projection")
    zm    = match(r"Zone_(\d+)([NS]?)", proj)
    zm === nothing && error("$f: cannot read a UTM zone from projection '$proj'")
    side  = parse(Float64, val("Side length"))
    cl    = filter(l -> occursin(r"_c\d+\s+[-\d.]+\s+[-\d.]+", l), L)
    sort!(cl, by = l -> parse(Int, match(r"_c(\d+)", l).captures[1]))
    xy    = [parse.(Float64, match(r"_c\d+\s+([-\d.]+)\s+([-\d.]+)", l).captures) for l in cl]
    length(xy) == 4 || error("$f: expected 4 antenna corners, found $(length(xy))")
    E = first.(xy); N = last.(xy)
    # every side should be close to the stated side length -- catches a
    # mistyped corner before it becomes a wrong kernel
    sides = [hypot(E[mod1(i+1,4)] - E[i], N[mod1(i+1,4)] - N[i]) for i in 1:4]
    maximum(abs.(sides .- side)) > 0.05side &&
        @warn "$f: corner spacing $(round.(sides, digits=1)) disagrees with side length $side"
    bearings = [mod(atand(E[mod1(i+1,4)] - E[i], N[mod1(i+1,4)] - N[i]), 360) for i in 1:4]
    (date = val("Date and time"), projection = proj, zone = parse(Int, zm[1]),
     hemisphere = isempty(zm[2]) ? "" : zm[2], side = side,
     easting = sum(E) / 4, northing = sum(N) / 4, bearings = bearings)
end

"Declination and inclination (degrees) from the AGRF25 output."
function read_agrf(f)
    L = readlines(f)
    i = findfirst(l -> occursin("Field", l) && occursin("=", l), L)
    i === nothing && error("$f: no 'Field =' line")
    _, _, _, _, F, D, I = parse.(Float64, split(L[i])[3:9])   # X Y Z H F D I
    lat = match(r"Latitude.*\(=\s*([-\d.]+)", join(L, "\n"))
    lon = match(r"Longitude.*\(=\s*([-\d.]+)", join(L, "\n"))
    (D = D, I = I, F = F, lat = lat === nothing ? NaN : parse(Float64, lat[1]),
     lon = lon === nothing ? NaN : parse(Float64, lon[1]))
end

"What the .mat holds, in either GMR layout."
function read_mat(f)
    d = matread(f)
    if haskey(d, "fid")
        fid = d["fid"]; ps = d["pulseSequence"]
        chans = sort([parse(Int, m[1]) for m in (match(r"^Channel_(\d+)$", k) for k in keys(fid)) if m !== nothing])
        live  = [c for c in chans if !all(iszero, fid["Channel_$c"])]
        (layout = "nested", pulse_s = Float64(ps["ExPulseLength_s"]),
         f0 = Float64(ps["DetectFrequency_Hz"]), seq = String(ps["SequenceType"]),
         nq = length(fid["pulseMoment"]["Magnitude"]), live = live)
    else
        live = [c for c in 1:8 if haskey(d, "coil_$(c)_fid") && !all(iszero, d["coil_$(c)_fid"])]
        (layout = "flat", pulse_s = Float64(d["T_pulse"]), f0 = Float64(d["detect_frequency"][1]),
         seq = "FID", nq = length(d["pulse_moment"]), live = live)
    end
end

function main(args)
    A   = parse_args(args)
    src = abspath(A["src"]); out = abspath(A["out"]); pre = get(A, "prefix", "")
    want = haskey(A, "sites") ? split(A["sites"], ',') : nothing
    files = allfiles(src)
    pfiles = filter(f -> occursin(r"_params\.txt$"i, basename(f)) &&
                         !occursin(r"process_params"i, basename(f)), files)
    isempty(pfiles) && error("no *_params.txt under $src")

    mkpath(joinpath(out, "snmr"))
    hdr = ["Field_ID", "easting", "northing", "zone", "azimuth", "bo_inc_angle_deg",
           "declination", "bo_dec_angle_deg", "coil_diameter", "detect_coil", "noise_coil",
           "pulse_sequence", "pulse_width", "year",
           # provenance, not read by the stages
           "site_name", "projection", "latitude", "longitude", "F_nT", "detect_freq_hz",
           "n_pulse_moments", "live_channels", "mat_layout", "src_mat", "src_params", "src_agrf"]
    rows = Vector{Vector{String}}()
    for pf in sort(pfiles)
        pdir = dirname(pf)
        site = replace(basename(pf), r"_params\.txt$"i => "")
        want !== nothing && site ∉ want && continue
        nested = lowercase(basename(pdir)) == "siteinfo"
        root = nested ? dirname(pdir) : pdir
        af = one_match(filter(f -> dirname(f) == pdir, files), r"_agrf25\.txt$"i, "*_agrf25.txt beside $(basename(pf))", site)
        mats = filter(f -> startswith(f, root * "/") && occursin(r"\.mat$"i, f) &&
                           !(nested && startswith(f, pdir * "/")), files)
        mf, passed = pick_mat(mats, root, site)
        isempty(passed) || println("  $site: chose $(relpath(mf, root)); passed over ",
                                   join(relpath.(passed, root), ", "))
        p, g, m = read_params(pf), read_agrf(af), read_mat(mf)
        isempty(m.live) && error("$site: no channel in $(basename(mf)) holds data")
        length(m.live) > 1 && @warn "$site: channels $(m.live) all hold data; using $(m.live[1]) as detect"

        # Loop orientation. azimuth = true bearing of the side nearest magnetic
        # north; bo_dec_angle_deg = declination - azimuth, which is how the UDF
        # sheet's columns relate (UD21: 9.42 - 45 = -35.58). smrp_invert.jl's
        # theta is the same angle with the opposite sign.
        dd  = [mod(b - g.D + 180, 360) - 180 for b in p.bearings]
        az  = p.bearings[argmin(abs.(dd))]
        rel = mod(g.D - az + 180, 360) - 180

        fid = pre * site
        dst = joinpath(out, "snmr", fid * "_" * replace(basename(mf), r"\s+" => "_"))
        # a site's earlier copies go first: udf_matfile picks by preference
        # among everything prefixed with the site, so a stale raw .mat from a
        # previous prep could be picked over nothing, or a stale _noemi over
        # a newer reprocessing
        for f in readdir(joinpath(out, "snmr"))
            startswith(f, fid * "_") && f != basename(dst) && rm(joinpath(out, "snmr", f))
        end
        cp(mf, dst; force = true)
        yr  = get(A, "year", string(year(DateTime(p.date[1:19], dateformat"yyyy-mm-dd HH:MM:SS"))))
        push!(rows, string.([fid, @sprintf("%.2f", p.easting), @sprintf("%.2f", p.northing), p.zone, round(az, digits = 3), g.I,
            g.D, round(rel, digits = 3), p.side, m.live[1], 0, m.seq == "FID" || occursin("FID", m.seq) ? "FID" : m.seq,
            round(1e3m.pulse_s, digits = 3), yr,
            site, p.projection, g.lat, g.lon, g.F, m.f0, m.nq, join(m.live, ";"), m.layout,
            relpath(mf, src), relpath(pf, src), relpath(af, src)]))
        @printf("%-10s zone %d%s  E %.0f N %.0f  side %.0f m  I %.2f D %.2f  az %.1f (rel %.2f)  %s %d q, %.0f ms, %.0f Hz, ch %s\n",
                fid, p.zone, p.hemisphere, p.easting, p.northing, p.side, g.I, g.D, az, rel,
                m.layout, m.nq, 1e3m.pulse_s, m.f0, join(m.live, ","))
    end
    isempty(rows) && error("no sites selected")
    q(s) = occursin(r"[,\"]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s
    open(joinpath(out, "snmr_metadata.csv"), "w") do io
        println(io, join(hdr, ","))
        foreach(r -> println(io, join(q.(r), ",")), rows)
    end
    zs = unique(r[4] for r in rows)
    length(zs) > 1 && @warn "sites span UTM zones $(join(zs, ", ")); AEM pairing must reproject"
    println("wrote $(length(rows)) site(s) -> $(joinpath(out, "snmr_metadata.csv"))")
end

main(ARGS)
