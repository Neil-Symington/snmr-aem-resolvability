# prep_aem_xcite.jl -- the nearest XCite sounding to each SMR site, as fieldrun AEM inputs.
#
#   julia --project=. fieldrun/prep/prep_aem_xcite.jl --config <toml> --dat <file>[,<file>...]
#         [--sites Mulan31,Mulan32]
#
# --dat takes the delivered ASEG .dat files, either as the delivery's .zip (one
# .dat inside; STREAMED with `unzip -p`, never extracted -- each is ~60 GB
# uncompressed) or as a plain .dat, such as a line subset cut earlier. The .dat
# is not sorted by line, so a whole block is one full pass: ~20 min (South) to
# ~40 min (Central). A subset of the right lines is seconds.
#
# Reads from the control file (TOML only: load_config validates [run] sites
# against sites.csv, which this script is what writes):
#   [paths] snmr_meta    the SMR sites (prep_smr.jl's snmr_metadata.csv)
#   [paths] aem_system   stem of the .des; the .dfn is found beside it
#   [aem]   nominal_unc, halt_noise_pv, max_sep_m, max_tilt_deg
# and writes:
#   [paths] aem_sounds/manifest.csv     one row per site, the columns aem_lib.jl's
#                                       manifest() reads, plus provenance
#   [paths] aem_sounds/<site>.csv       moment, gate, d_pV, halt_pV, sigma_pV
#   [paths] sites                       Field_ID + the SMR <-> AEM pairing
#
# ONE SOUNDING PER SITE (Carl, 2026-11-12): the nearest record that passes QC,
# not a stack. QC: all gates of dBdt_Z_F present, tilt sqrt(pitch^2 + roll^2)
# <= [aem] max_tilt_deg (where attitude was recorded; see below), and within
# [aem] max_sep_m. The nearest record on each
# OTHER line within max_sep_m is printed too: Mulan31 has a repeat pass.
#
# NOISE: sigma_i = sqrt((nominal_unc * d_i)^2 + halt_i^2), the HiQGA VTEM form
# (VTEM1DInversion.read_survey_files). Units pV/(A.m^4) throughout.
# See staging/data/aem/XCITE_DECISIONS.md.

using TOML, Printf, DelimitedFiles

const FIELDRUN = normpath(joinpath(@__DIR__, ".."))

function parse_args(a)
    d = Dict{String,String}(); i = 1
    while i <= length(a)
        startswith(a[i], "--") && i < length(a) || error("usage: --config <toml> --dat <file>[,...] [--sites A,B]")
        d[a[i][3:end]] = a[i+1]; i += 2
    end
    haskey(d, "config") && haskey(d, "dat") || error("usage: --config <toml> --dat <file>[,...] [--sites A,B]")
    d
end

resolve(p) = isabspath(p) ? p : normpath(joinpath(FIELDRUN, p))

"""
Field positions from the .dfn: name => first column (1-based) and width.
`24F13.5` is 24 columns, `F12.2` one. The leading RT=COMM record is skipped.
"""
function dfn_layout(dfn)
    lay = Dict{String,Tuple{Int,Int}}(); col = 1
    for l in readlines(dfn)
        m = match(r"^DEFN\s+\d+\s+ST=RECD,RT=;([A-Za-z0-9_]+):(\d*)[A-Z]", l)
        m === nothing && continue
        n = isempty(m[2]) ? 1 : parse(Int, m[2])
        lay[m[1]] = (col, n); col += n
    end
    isempty(lay) && error("$dfn: no field definitions parsed")
    lay, col - 1
end

function find_dfn(stem)
    d, b = dirname(stem), basename(stem)
    c = filter(f -> startswith(f, b) && endswith(lowercase(f), ".dfn"), readdir(d))
    length(c) == 1 || error("expected one .dfn beside $stem, found $(length(c)): $c")
    joinpath(d, only(c))
end

function read_sites(meta, want)
    raw, hdr = readdlm(meta, ',', String; header = true)
    ci = Dict(strip(h) => i for (i, h) in enumerate(vec(hdr)))
    S = [(id = raw[i, ci["Field_ID"]], e = parse(Float64, raw[i, ci["easting"]]),
          n = parse(Float64, raw[i, ci["northing"]]), zone = raw[i, ci["zone"]])
         for i in axes(raw, 1)]
    want === nothing ? S : filter(s -> s.id in want, S)
end

openstream(f) = endswith(lowercase(f), ".zip") ? open(`unzip -p $f`) : open(f)

function main(args)
    A = parse_args(args)
    C = TOML.parsefile(A["config"]); P = C["paths"]; AE = C["aem"]
    get(AE, "system", "skytem") == "xcite" || error("$(A["config"]): [aem] system is not \"xcite\"")
    stem  = resolve(P["aem_system"])
    lay, ncol = dfn_layout(find_dfn(stem))
    col(k) = (haskey(lay, k) || error(".dfn has no field $k"); lay[k][1])
    iL, iF, iFl, iD = col("Line"), col("Fiducial"), col("Flight"), col("Date")
    iE, iN = col("Easting"), col("Northing")
    iTx, iRx, iP, iR = col("Tx_Loop_Clearance"), col("Rx_Clearance"), col("Tx_Loop_Pitch"), col("Tx_Loop_Roll")
    iG, iPLM = col("Ground_Elevation"), col("PLM")
    z0, nz = lay["dBdt_Z_F"]
    halt = Float64.(AE["halt_noise_pv"]); mult = Float64(AE["nominal_unc"])
    nz == length(halt) || error("dBdt_Z_F has $nz gates but [aem] halt_noise_pv has $(length(halt))")
    gates = haskey(AE, "gates") ? collect(AE["gates"][1]:AE["gates"][2]) : collect(0:nz-1)
    length(gates) == nz || error("[aem] gates $(AE["gates"]) is not $nz gates")
    maxsep = Float64(AE["max_sep_m"]); maxtilt = Float64(get(AE, "max_tilt_deg", 10.0))
    NULL = -99999.0

    want = haskey(A, "sites") ? String.(split(A["sites"], ',')) : nothing
    sites = read_sites(resolve(P["snmr_meta"]), want)
    isempty(sites) && error("no sites selected")
    length(unique(s.zone for s in sites)) == 1 || error("SMR sites span UTM zones; the AEM is one zone")
    @printf("%d site(s); %s: %d fields per record, dBdt_Z_F at %d-%d; max sep %.0f m, max tilt %.0f deg\n",
            length(sites), basename(find_dfn(stem)), ncol, z0, z0 + nz - 1, maxsep, maxtilt)

    R2 = maxsep^2
    best = Dict(s.id => (Inf, String[], "", 0) for s in sites)           # (d, fields, source, rec)
    perline = Dict(s.id => Dict{Int,Tuple{Float64,Int}}() for s in sites) # line => (d, fid)
    nnear = Dict(s.id => 0 for s in sites); ntilt = Dict(s.id => 0 for s in sites)
    nnull = Dict(s.id => 0 for s in sites); nnoatt = Dict(s.id => 0 for s in sites)
    for f in split(A["dat"], ',')
        t0 = time(); nrec = 0
        io = openstream(f)
        for l in eachline(io)
            nrec += 1
            # cheap reject on E/N before splitting the whole record: almost
            # every record of a 60 GB block is far from every site
            h = split(l; limit = max(iE, iN) + 1)
            e = parse(Float64, h[iE]); n = parse(Float64, h[iN])
            any(s -> (e - s.e)^2 + (n - s.n)^2 <= R2, sites) || continue
            fs = split(l)
            length(fs) == ncol || error("$f record $nrec has $(length(fs)) fields, the .dfn says $ncol")
            for s in sites
                d2 = (e - s.e)^2 + (n - s.n)^2
                d2 <= R2 || continue
                nnear[s.id] += 1
                z = parse.(Float64, fs[z0:z0+nz-1])
                any(<=(NULL + 1), z) && (nnull[s.id] += 1; continue)
                pr = (parse(Float64, fs[iP]), parse(Float64, fs[iR]))
                # pitch = roll = 0.0 EXACTLY is "not recorded", not "level":
                # all 123,965 records of line 40052 (flight 64) read that way
                # while every neighbouring line carries real attitude. The tilt
                # check cannot be applied there, so the record is kept and
                # flagged rather than refused or passed as level.
                pr == (0.0, 0.0) && (nnoatt[s.id] += 1)
                pr != (0.0, 0.0) && hypot(pr...) > maxtilt && (ntilt[s.id] += 1; continue)
                d = sqrt(d2); ln = parse(Int, fs[iL])
                (!haskey(perline[s.id], ln) || d < perline[s.id][ln][1]) &&
                    (perline[s.id][ln] = (d, parse(Int, fs[iF])))
                d < best[s.id][1] && (best[s.id] = (d, String.(fs), basename(f), nrec))
            end
        end
        close(io)
        @printf("  %s: %d records in %.1f min\n", basename(f), nrec, (time() - t0) / 60)
    end

    snd = resolve(P["aem_sounds"]); mkpath(snd)
    mrows = String[]; srows = String[]
    for s in sites
        d, fs, src, rec = best[s.id]
        @printf("%-10s %d records within %.0f m (%d tilted > %.0f deg, %d with null gates, %d with no attitude)\n",
                s.id, nnear[s.id], maxsep, ntilt[s.id], maxtilt, nnull[s.id], nnoatt[s.id])
        if !isfinite(d)
            @warn "$(s.id): no QC-passing XCite sounding within $maxsep m -- not written"
            continue
        end
        for (ln, (dl, fid)) in sort(collect(perline[s.id]), by = x -> x[2][1])
            @printf("    line %d: nearest %.0f m (fiducial %d)%s\n", ln, dl, fid,
                    dl == d ? "  <- used" : "")
        end
        g(i) = parse(Float64, fs[i])
        z = parse.(Float64, fs[z0:z0+nz-1])
        σ = sqrt.((mult .* z).^2 .+ halt.^2)
        open(joinpath(snd, "$(s.id).csv"), "w") do io
            println(io, "moment,gate,d_pV,halt_pV,sigma_pV")
            for k in 1:nz
                @printf(io, "Z,%d,%.6g,%.6g,%.6g\n", gates[k], z[k], halt[k], σ[k])
            end
        end
        att = (g(iP), g(iR)) == (0.0, 0.0) ? "not_recorded" : "recorded"
        tilt = att == "recorded" ? hypot(g(iP), g(iR)) : NaN
        att == "recorded" || @warn "$(s.id): the chosen sounding has no attitude (pitch = roll = 0.0 " *
            "exactly), so the $(maxtilt)-degree tilt check could not be applied to it"
        push!(mrows, join([s.id, "$(s.id).csv", fs[iL], fs[iF], string(rec), fs[iTx], fs[iG],
                           fs[iE], fs[iN], @sprintf("%.1f", d), fs[iPLM], string(nz),
                           "none", "NaN", "NaN",
                           fs[iRx], fs[iP], fs[iR], @sprintf("%.1f", tilt), att, fs[iFl], fs[iD], src], ","))
        push!(srows, join([s.id, @sprintf("%.2f", s.e), @sprintf("%.2f", s.n), s.zone,
                           @sprintf("%.1f", d), fs[iL], fs[iF], fs[iTx], string(rec), src], ","))
        @printf("    used: line %s fiducial %s, %.0f m, Tx %s m, tilt %s, %s record %d\n",
                fs[iL], fs[iF], d, fs[iTx], isnan(tilt) ? "not recorded" : @sprintf("%.1f deg", tilt), src, rec)
    end
    isempty(mrows) && error("no site has a usable sounding")
    open(joinpath(snd, "manifest.csv"), "w") do io
        println(io, "field_id,csv,line,fiducial,rec,height_m,dtm_m,easting,northing,sep_m,plm,n_z," *
                    "bore,bore_sep_m,bore_td_m,rx_height_m,pitch_deg,roll_deg,tilt_deg,attitude,flight,date,source")
        foreach(r -> println(io, r), mrows)
    end
    sf = resolve(P["sites"]); mkpath(dirname(sf))
    open(sf, "w") do io
        println(io, "Field_ID,easting,northing,zone,aem_sep_m,aem_line,aem_fiducial,aem_tx_clearance_m,aem_rec,aem_source")
        foreach(r -> println(io, r), srows)
    end
    println("wrote $(length(mrows)) sounding(s) -> $snd, pairing -> $sf")
end

main(ARGS)
