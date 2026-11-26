# prep_smr.jl -- turn a GMR field delivery into the fieldrun SNMR inputs.
#
#   julia --project=. fieldrun/prep/prep_smr.jl --src <delivery dir>[,<dir>...] --out <prepared dir>
#         [--sites Mulan31,Mulan32] [--prefix ""] [--year 2026]
#         [--emi cancelled|uncancelled] [--strict true] [--folder Reprocessed_11_26]
#
# --emi picks WHICH processed file per site (default "cancelled", as before).
# "uncancelled" takes the _noemi file instead -- the comparison set for the
# reference-loop signal-loss check. With --strict true a site that lacks the
# asked-for variant is SKIPPED (and listed) rather than falling back to the
# other, so the two prepared sets hold like-for-like files only.
# --folder restricts the .mat to one reprocessing folder (a path component),
# instead of "the newest Reprocessed_* folder".
#
# SITE NAMES. In the delivery layout a site is named after its FOLDER (the
# parent of SiteInfo), not its params file: at Birrindudu, Billiluna09 and
# Billiluna09_2DCDC both carry Billiluna09_params.txt, S_Mulan02's is
# S_Mulan0_params.txt and Jalywon1's is Jalywon-v2_params.txt.
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
# Several .mat per site is normal (with and without EMI cancellation, repeat
# reprocessings). The choice: the newest Reprocessed_* folder if there are
# several, then EMI-CANCELLED over not, then harmonics removed over not. What
# was chosen and what was passed over are both printed.
#
# "_noemi" MEANS NO EMI CANCELLATION (Carl, 2026-11-13) -- it is the file to
# AVOID, not to prefer. An earlier version of this script copied udf_smr.jl's
# order, harmonic_removed_noemi > noemi > harmonic_removed > raw, which reads
# "noemi" as "EMI removed" and so picked the UNcancelled file at every site.
# The Birrindudu sites are near a town and noisy, so the cancellation matters.
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
    haskey(d, "src") && haskey(d, "out") || error("usage: --src <dir>[,<dir>] --out <dir> [--sites A,B] [--prefix P] [--emi cancelled|uncancelled] [--strict true] [--folder F]")
    get(d, "emi", "cancelled") in ("cancelled", "uncancelled") || error("--emi must be cancelled or uncancelled")
    d
end

# Not walkdir: it trusts the directory entry's type, and an NTFS drive read
# through a third-party macOS driver (Birrindudu delivery, Nov 2026) reports
# some folders -- SiteInfo among them -- with no type, which walkdir lists as
# files and never enters. isdir() stats each entry instead.
function allfiles(dir)
    out = String[]
    for f in readdir(dir)
        startswith(f, ".") && continue
        p = joinpath(dir, f)
        isdir(p) ? append!(out, allfiles(p)) : push!(out, p)
    end
    out
end

function one_match(files, re, what, site)
    m = filter(f -> occursin(re, basename(f)), files)
    isempty(m) && error("$site: no $what")
    length(m) > 1 && error("$site: $(length(m)) candidates for $what, expected one:\n  " * join(m, "\n  "))
    only(m)
end

"""
The processed .mat for a site: newest Reprocessed_* folder first, then
EMI-cancelled (no "_noemi"), then harmonics removed. Returns (chosen, passed_over).
"""
function pick_mat(mats, root, site; emi = "cancelled", strict = false)
    isempty(mats) && error("$site: no .mat under $root")
    # newest reprocessing: Reprocessed_<mm>_<yy>, compared as (yy, mm)
    stamp(f) = (m = match(r"Reprocessed_(\d{1,2})_(\d{2})", relpath(f, root));
                m === nothing ? (-1, -1) : (parse(Int, m[2]), parse(Int, m[1])))
    newest = maximum(stamp.(mats))
    pool = filter(f -> stamp(f) == newest, mats)
    isnoemi(f) = occursin("noemi", basename(f))
    want_noemi = emi == "uncancelled"
    rank(f) = (isnoemi(f) == want_noemi ? 0 : 1,                    # the asked-for variant first
               occursin("harmonic_removed", basename(f)) ? 0 : 1)  # then harmonics removed
    best = minimum(rank.(pool))
    c = filter(f -> rank(f) == best, pool)
    length(c) > 1 && error("$site: $(length(c)) equally preferred .mat files:\n  " * join(c, "\n  "))
    if best[1] == 1
        strict && return nothing, mats
        @warn "$site: no $(emi == "cancelled" ? "EMI-cancelled" : "_noemi (uncancelled)") file; using $(basename(only(c)))"
    end
    return only(c), setdiff(mats, c)
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
    # Corners from the POINT (E N) geometry column where it is there: it is
    # never truncated and carries the decimals. The easting/northing columns
    # are a pandas print, which elides middle columns as "..." when the table
    # is wide (BDU Mulan37: no northing column at all) and rounds to 1 m.
    cl    = filter(l -> occursin(r"_c\d+\s", l) &&
                        (occursin(r"POINT \(", l) || occursin(r"_c\d+\s+[-\d.]+\s+[-\d.]+", l)), L)
    sort!(cl, by = l -> parse(Int, match(r"_c(\d+)", l).captures[1]))
    corner(l) = (m = match(r"POINT \(([-\d.]+) ([-\d.]+)\)", l);
                 m === nothing && (m = match(r"_c\d+\s+([-\d.]+)\s+([-\d.]+)", l));
                 parse.(Float64, m.captures))
    xy    = [corner(l) for l in cl]
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
    srcs = abspath.(split(A["src"], ',')); out = abspath(A["out"]); pre = get(A, "prefix", "")
    emi = get(A, "emi", "cancelled"); strict = get(A, "strict", "false") == "true"
    folder = get(A, "folder", "")
    want = haskey(A, "sites") ? split(A["sites"], ',') : nothing
    files = reduce(vcat, allfiles.(srcs))
    srcof(f) = srcs[findfirst(s -> startswith(f, s * "/"), srcs)]
    relsrc(f) = (s = srcof(f); length(srcs) > 1 ? joinpath(basename(s), relpath(f, s)) : relpath(f, s))
    skipped = String[]
    pfiles = filter(f -> occursin(r"_params\.txt$"i, basename(f)) &&
                         !occursin(r"process_params"i, basename(f)), files)
    # Delivery layout: only SiteInfo/*_params.txt names a site. Acquisition
    # folders can carry their own params copies (BDU Jalywon1/Partial
    # acquisition/...), which would otherwise be read as flat-layout sites.
    # And only <src>/<site>/SiteInfo: a nested SiteInfo inside an acquisition
    # folder is a partial acquisition, not a site.
    insi = filter(f -> lowercase(basename(dirname(f))) == "siteinfo" &&
                       dirname(dirname(dirname(f))) in srcs, pfiles)
    if !isempty(insi) && length(insi) < length(pfiles)
        println("ignoring params files outside SiteInfo: ",
                join(relsrc.(setdiff(pfiles, insi)), ", "))
        pfiles = insi
    end
    isempty(pfiles) && error("no *_params.txt under $(join(srcs, ", "))")

    mkpath(joinpath(out, "snmr"))
    hdr = ["Field_ID", "easting", "northing", "zone", "azimuth", "bo_inc_angle_deg",
           "declination", "bo_dec_angle_deg", "coil_diameter", "detect_coil", "noise_coil",
           "pulse_sequence", "pulse_width", "year",
           # provenance, not read by the stages
           "site_name", "projection", "latitude", "longitude", "F_nT", "detect_freq_hz",
           "n_pulse_moments", "live_channels", "mat_layout", "src_mat", "src_params", "src_agrf", "emi"]
    rows = Vector{Vector{String}}()
    for pf in sort(pfiles)
        pdir = dirname(pf)
        nested = lowercase(basename(pdir)) == "siteinfo"
        root = nested ? dirname(pdir) : pdir
        # delivery layout: the site is its folder; flat: its params file
        site = nested ? basename(root) : replace(basename(pf), r"_params\.txt$"i => "")
        want !== nothing && site ∉ want && continue
        af = one_match(filter(f -> dirname(f) == pdir, files), r"_agrf25\.txt$"i, "*_agrf25.txt beside $(basename(pf))", site)
        mats = filter(f -> startswith(f, root * "/") && occursin(r"\.mat$"i, f) &&
                           !(nested && startswith(f, pdir * "/")), files)
        isempty(folder) || (mats = filter(f -> folder in splitpath(relpath(f, root)), mats))
        mf, passed = pick_mat(mats, root, site; emi = emi, strict = strict)
        if mf === nothing
            println("  $site: SKIPPED, no $emi file (--strict): has ", join(relpath.(passed, root), ", "))
            push!(skipped, site); continue
        end
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
            relsrc(mf), relsrc(pf), relsrc(af),
            occursin("noemi", basename(mf)) ? "uncancelled" : "cancelled"]))
        @printf("%-10s zone %d%s  E %.0f N %.0f  side %.0f m  I %.2f D %.2f  az %.1f (rel %.2f)  %s %d q, %.0f ms, %.0f Hz, ch %s\n",
                fid, p.zone, p.hemisphere, p.easting, p.northing, p.side, g.I, g.D, az, rel,
                m.layout, m.nq, 1e3m.pulse_s, m.f0, join(m.live, ","))
    end
    isempty(rows) && error("no sites selected")
    sids = [r[1] for r in rows]
    length(unique(sids)) == length(sids) ||
        error("duplicate Field_ID: $(join(unique(filter(x -> count(==(x), sids) > 1, sids)), ", "))")
    isempty(skipped) || println("skipped (no $emi file): ", join(skipped, ", "))
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
