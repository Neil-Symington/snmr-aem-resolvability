# udf_system.jl — the Upper Darling Floodplain SkyTEM 312 system, as delivered.
#
# Survey AUS_10087 (GA project 5011), SkyTEM Australia, flown March–July 2022.
# Sources, both of which agree:
#
#   * the .des named by [paths] aem_system -- shipped as
#     example/aem/system/UDF_AEM_Data.des, copied verbatim from the delivery
#   * AUS_10087_report.pdf, Tables 9 and 10 (gate times) and p.~19 (system)
#
# The gate windows and both waveforms are PARSED from the .des at load time, and
# the gate table is cross-checked against itself: the .des lists open, centre and
# close, and centre must be the mean of open and close. `check_against_config`
# enforces that, which is what catches a single-value edit.
#
# The GEOMETRY, though, IS transcribed: rTx, rRx, dzRx, lowpassfcs and both ramp
# holds come from the control file, not from the .des, and nothing compares them
# against it. An earlier version of this header claimed "nothing here is
# transcribed by hand", which a review showed was false for exactly those
# values. Changing tx_area_m2 or rx_radius_m in the config is accepted in
# silence. Treat them as declared, not verified. A transcription error in a gate table is invisible until it
# quietly biases every depth in the survey, and the Ord-Keep system file had to
# warn that HiQGA's bundled SkyTEM example is a different survey's system. The
# same trap applies here: this survey's LM gate 19 closes at 203.2 us where
# Ord-Keep's closes at 192.2 us.
#
# Gate times are referenced to the initiation of current switch-off. The report
# states a calibration correction shift has ALREADY been applied to the
# delivered gate times, so none is applied here.
#
# Gate centres are the GEOMETRIC mean of open and close, which is HiQGA's
# convention for SkyTEM. The .des also supplies an ARITHMETIC centre; the two
# differ by under 1% for the narrow early gates but by up to 4% for the wide
# late ones, so the choice is recorded rather than assumed. `report_centres()`
# returns the delivered arithmetic centres for comparison.

module UDFSystem

using HiQGA.transD_GP
using Statistics, DelimitedFiles
using ..FieldConfig

export LM_times, HM_times, LM_ramp, HM_ramp, rTx, rRx, dzRx, lowpassfcs,
       LM_GATES, HM_GATES, LM_WINDOW_US, HM_WINDOW_US,
       build_grid, build_operator, PV, MU0, DES, report_centres, describe

const MU0 = 4π * 1e-7
const PV  = 1e-12                 # delivered units: pV/(A.turns.m^4)

# The .des is named by the control file. This module PARSES it at load time --
# deliberately, so a transcription error cannot survive -- so FieldConfig must
# be loaded first. Including these the other way round gives the message below
# rather than an UndefVarError.
const DES = let c = FieldConfig.CFG[]
    c === nothing && error("udf_system.jl needs the control file: include " *
                           "common/fieldconfig.jl and call load_config first")
    c["paths"]["aem_system"] * ".des"
end
const _AEM = FieldConfig.CFG[]["aem"]

"Strip the leading COMM and return the .des body as lines."
desbody() = [replace(l, r"^COMM\s?" => "") for l in readlines(DES)]

"""
Parse one of the two gate tables.

The .des blocks are headed `LM GATE TIMES SUPPLIED` / `HM GATE TIMES SUPPLIED`
and carry columns: gate number, width, open, centre, close (microseconds).
Returns (gates, window_us, centre_us) with window = [open close].
"""
function parse_gates(tag)
    L = desbody()
    i = findfirst(l -> occursin("$tag GATE TIMES SUPPLIED", l), L)
    i === nothing && error("no $tag gate table in $DES")
    # wd is the WIDTH column. It was parsed by the regex and discarded, so
    # nothing could check it -- a reviewer changed gate 19's width from 50.57 to
    # 99.99 and every check passed.
    n, w, c, wd = Int[], Matrix{Float64}(undef, 0, 2), Float64[], Float64[]
    for l in L[i:end]
        occursin("---", l) && !isempty(n) && break
        m = match(r"^\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", l)
        m === nothing && continue
        push!(n, parse(Int, m[1]))
        push!(wd, parse(Float64, m[2]))
        w = vcat(w, [parse(Float64, m[3]) parse(Float64, m[5])])
        push!(c, parse(Float64, m[4]))
    end
    isempty(n) && error("$tag gate table parsed empty")
    (n, w, c, wd)
end

"""
Parse one of the two transmitter waveforms.

The .des lists `Time[s]  Amplitude` pairs for the LM and HM ramps. HiQGA wants
the waveform up to and including switch-off, with time measured from the same
zero as the gates — the top of the ramp-down, which is where the .des listing
has amplitude 1.0. Only the ramp DOWN is kept (t >= 0), prefixed by the
switch-off start, which is what `dBzdt` integrates over.
"""
function parse_waveform(tag)
    L = desbody()
    i = findfirst(l -> occursin("$tag Tx WAVEFORM", l) ||
                       occursin("$tag Tx Waveform", l), L)
    i === nothing && error("no $tag waveform in $DES")
    t, a = Float64[], Float64[]
    seen = false
    for l in L[i+1:end]
        m = match(r"^\s*(-?[\d.]+E[+-]\d+)\s+(-?[\d.]+E[+-]\d+)", l)
        if m === nothing
            seen && occursin("---", l) && break
            continue
        end
        seen = true
        push!(t, parse(Float64, m[1])); push!(a, parse(Float64, m[2]))
    end
    isempty(t) && error("$tag waveform parsed empty")
    (t, a)
end

"""
The ramp HiQGA needs: the whole final pulse, ON-TIME RAMP-UP INCLUDED, from the
last zero before the current peak through to the first zero after it.

The on-time limb is not optional. HiQGA's own SkyTEM example
(examples/SkyTEM1D/Menindee/.../electronics_halt.jl) gives `LM_ramp` starting
at -0.0008 s with amplitude 0, and ordk_system.jl follows it. Keeping only the
switch-off limb tells the forward that the current had been on forever, which
over-predicts the late-time response — invisible over the resistive Ord-Keep
ground, but on this floodplain (~0.5 S/m at 20-60 m) it over-predicted LM gate
26 by 52% and forced a 6-8x noise inflation at every site.

Times keep the listing's own base, which is the same zero the gates use.
"""
function rampof(tag; hold)
    t, a = parse_waveform(tag)
    ipk = findlast(x -> x >= 0.999, a)                  # top of the ramp down
    i0 = findlast(j -> a[j] <= 1e-12, 1:ipk)            # start of the ramp up
    i0 === nothing && error("$tag waveform has no zero before the peak")
    i1 = findfirst(j -> a[j] <= 1e-12, ipk:length(a))
    i1 === nothing && error("$tag waveform has no switch-off zero")
    i1 = ipk + i1 - 1
    R = [t[i0:i1] a[i0:i1]]
    vcat(R, [hold 0.0])                     # hold at zero to the next pulse
end

const _LM = parse_gates("LM")
const _HM = parse_gates("HM")
const LM_GATES, LM_WINDOW_US = _LM[1], _LM[2]
const HM_GATES, HM_WINDOW_US = _HM[1], _HM[2]

geomcentre(W) = sqrt.(W[:, 1] .* W[:, 2]) .* 1e-6
const LM_times = geomcentre(LM_WINDOW_US)
const HM_times = geomcentre(HM_WINDOW_US)

"The arithmetic gate centres the .des and the report supply, in seconds."
report_centres() = (lm = _LM[3] .* 1e-6, hm = _HM[3] .* 1e-6)

# Off-times: LM 1.018 ms, HM 15.0 ms. Also declared in the control file rather
# than parsed; the ramp SHAPES are parsed, the holds are not.
const LM_ramp = rampof("LM"; hold = _AEM["lm_ramp_hold"])
const HM_ramp = rampof("HM"; hold = _AEM["hm_ramp_hold"])

# ── geometry: DECLARED in the control file, not parsed ──────────────────────
# These values are from AUS_10087_report.pdf, which is not shipped. The .des
# does NOT contain them -- grep it for 342, 13.29 or 203.6 and you get
# nothing -- so there is no file here to check them against, and a wrong
# tx_area_m2 or rx_radius_m in the config is accepted in silence. An earlier
# version of this comment cited .des sections "INSTRUMENT POSITION" and
# "TRANSMITTER SPECIFICATIONS" that do not exist in the file.
# Tx loop area 342 m^2 (8-vertex polygon) -> equivalent circular radius.
const rTx = sqrt(_AEM["tx_area_m2"] / π)       # 10.43 m for 342 m^2
# Z-coil at x = -13.29 m, z = -2.00 m from the loop centre, in a frame whose z
# is POSITIVE BELOW the Tx loop, so the receiver sits 2 m above the Tx plane.
const rRx = _AEM["rx_radius_m"]
const dzRx = _AEM["rx_dz_m"]
# Rx coil cut-off 203.6 kHz (Z); Rx electronics low pass 300 kHz.
# (The X coil is 340 kHz; only the Z component is inverted here.)
const lowpassfcs = Float64.(_AEM["lowpass_hz"])

"""
    build_grid(; dz = 1.5, extendfrac = 1.03, n = 70)

Log-growing depth grid. The UDF ground is much more conductive than Ord-Keep —
the report caps SMR interpretation at ~90 m — so the grid is unchanged but the
useful depth will be shallower.
"""
function build_grid(; dz = 1.5, extendfrac = 1.03, n = 70)
    zall, znall, zbnd = transD_GP.setupz(0.0, extendfrac; dz = dz, n = n,
                                         showplot = false)
    z, ρ, nfixed = transD_GP.makezρ(zbnd; zfixed = [-1e5], ρfixed = [1e12])
    (zall = zall, znall = znall, zbnd = zbnd, z = z, ρ = ρ, nfixed = nfixed)
end

"""
    build_operator(d_lm, σ_lm, d_hm, σ_hm, height; ...)

`SkyTEM1DInversion.dBzdt` for one fiducial. `d_*`/`σ_*` are Z-component
dB/dt in pV/(A.turns.m^4); `height` is the laser height of the Tx loop centre.
"""
function build_operator(d_lm, σ_lm, d_hm, σ_hm, height::Real;
                        grid = build_grid(), calcjacobian = true,
                        σ_scale::Real = 1.0, ρ_start::Real = 30.0)
    zTx = -float(height)
    zRx = zTx + dzRx
    k = PV / MU0
    # The AIR layer must keep its 1e12 ohm.m: getfield! only overwrites layers
    # nfixed+1:end, and a conductive "air" halfspace makes the forward return
    # exactly zero, silently. See project_hiqga_air_layer_trap.
    ρ0 = fill(float(ρ_start), length(grid.z))
    ρ0[1:grid.nfixed] .= grid.ρ[1:grid.nfixed]
    transD_GP.SkyTEM1DInversion.dBzdt(;
        timeslow = LM_times, ramplow = LM_ramp, zRxlow = zRx, zTxlow = zTx,
        timeshigh = HM_times, ramphigh = HM_ramp, zRxhigh = zRx, zTxhigh = zTx,
        rRx = rRx, rTx = rTx,
        z = grid.z, ρ = ρ0,
        dlow = d_lm .* k, dhigh = d_hm .* k,
        σlow = σ_lm .* (σ_scale * k), σhigh = σ_hm .* (σ_scale * k),
        lowpassfcs = lowpassfcs, calcjacobian = calcjacobian)
end

"Print what was parsed, so the system can be eyeballed against the report."
function describe()
    rc = report_centres()
    println("UDF SkyTEM 312 (AUS_10087), parsed from $(basename(DES))")
    println("  LM gates $(first(LM_GATES))-$(last(LM_GATES)) ($(length(LM_times))), ",
            "HM gates $(first(HM_GATES))-$(last(HM_GATES)) ($(length(HM_times)))")
    for (n, t, c, W) in (("LM", LM_times, rc.lm, LM_WINDOW_US),
                         ("HM", HM_times, rc.hm, HM_WINDOW_US))
        d = maximum(abs.(t .- c) ./ c) * 100
        println("  $n window $(W[1,1]) - $(W[end,2]) us; geometric vs delivered ",
                "centres differ by at most $(round(d, digits=2))%")
    end
    println("  LM ramp $(size(LM_ramp,1)) points, 0 to $(LM_ramp[end-1,1]) s")
    println("  HM ramp $(size(HM_ramp,1)) points, 0 to $(HM_ramp[end-1,1]) s")
    println("  rTx $(round(rTx, digits=2)) m, rRx $rRx m, dzRx $dzRx m, ",
            "lowpass $(lowpassfcs) Hz")
end


"""
    check_against_config()

Assert the parsed .des agrees with what the control file declares. The point of
parsing rather than transcribing is lost if a config can quietly disagree with
the file, so [aem] lm_gates/hm_gates are treated as a declared cross-check and
not as a source of truth.
"""
function check_against_config()
    for (tag, gates, key) in (("LM", LM_GATES, "lm_gates"), ("HM", HM_GATES, "hm_gates"))
        lo, hi = _AEM[key]
        want = hi - lo + 1
        length(gates) == want || error(
            "[aem] $key declares gates $lo-$hi ($want gates) but the .des " *
            "supplies $(length(gates)) for $tag")
        (first(gates) == lo && last(gates) == hi) || error(
            "[aem] $key declares gates $lo-$hi but the .des $tag table runs " *
            "$(first(gates))-$(last(gates))")
    end
    issorted(LM_times) || error("LM gate centres are not increasing")
    issorted(HM_times) || error("HM gate centres are not increasing")

    # THE CHECK THIS FILE'S HEADER PROMISED AND DID NOT HAVE. The header says
    # the gates are "cross-checked against the report's own centre times";
    # nothing did that, and a reviewer changed LM gate 19's close time from
    # 203.2 to 192.2 us -- the exact transcription trap the header names, by
    # name -- and check_against_config accepted it.
    #
    # The .des gives open, centre and close for every gate, so the table can be
    # checked against ITSELF: centre must be the arithmetic mean of open and
    # close, and it holds to the printed digit on the shipped delivery. Any
    # single-value edit to one of the three columns breaks it.
    # A review measured this check catching ONE of eleven .des edits. The
    # shipped file satisfies every identity below to machine precision --
    # |mid - centre|/centre is 2.0e-16 (LM) and 1.3e-16 (HM), and
    # |(close-open) - width| is under 2e-12 -- so a 1e-3 tolerance left
    # 5e12 times more slack than needed, and three of the four columns were
    # unchecked. What got through: a mean-preserving edit (open -5, close +5),
    # a +0.3 us edit under the tolerance, the WIDTH column entirely, all three
    # columns shifted together so gates overlap, a duplicate gate number, and
    # any change to the waveform.
    #
    # Tolerances are now 1e-9 relative, which the delivery satisfies with eleven
    # orders of magnitude to spare, and all four columns plus the gate numbering
    # and the waveform are checked.
    TOL = 1e-9
    for (tag, N, W, C, Wd) in (("LM", LM_GATES, LM_WINDOW_US, _LM[3], _LM[4]),
                               ("HM", HM_GATES, HM_WINDOW_US, _HM[3], _HM[4]))
        # the gate NUMBER, not the row index: the old message printed argmax(d),
        # a row, labelled "gate", so a reader searching the .des landed on the
        # wrong line
        bad(d) = (i = argmax(d); (i, N[i], d[i]))

        # 1. centre is the mean of open and close
        i, g, dv = bad(abs.((W[:, 1] .+ W[:, 2]) ./ 2 .- C) ./ max.(abs.(C), eps()))
        dv <= TOL || error(
            "$tag gate $g: open $(W[i,1]) us and close $(W[i,2]) us give a centre " *
            "of $((W[i,1]+W[i,2])/2) us, but the .des lists $(C[i]) us " *
            "(relative difference $(round(dv, sigdigits=3))). One of the three " *
            "time columns has been altered.")

        # 2. width is close minus open -- a column nothing used to read
        i, g, dv = bad(abs.((W[:, 2] .- W[:, 1]) .- Wd) ./ max.(abs.(Wd), eps()))
        dv <= TOL || error(
            "$tag gate $g: close minus open is $(W[i,2]-W[i,1]) us but the .des " *
            "lists a width of $(Wd[i]) us (relative difference " *
            "$(round(dv, sigdigits=3)))")

        # 3. ordering, and NO OVERLAP -- the old comment claimed non-overlap and
        #    checked only that opens increase, so shifting one gate onto its
        #    neighbours passed
        all(W[:, 1] .< W[:, 2]) || error("$tag has a gate whose open is not before its close")
        for j in 1:size(W, 1) - 1
            W[j, 2] <= W[j+1, 1] + TOL * max(abs(W[j,2]), 1) || error(
                "$tag gates $(N[j]) and $(N[j+1]) overlap: $(N[j]) closes at " *
                "$(W[j,2]) us and $(N[j+1]) opens at $(W[j+1,1]) us")
        end

        # 4. gate numbers must be consecutive and unique
        N == collect(first(N):last(N)) || error(
            "$tag gate numbers are not consecutive and unique: $(N)")
    end

    # 5. the waveforms. Neither ramp was checked at all, so a reviewer doubled
    #    LM's switch-off time and changed its peak amplitude without complaint.
    for (tag, R, hold) in (("LM", LM_ramp, _AEM["lm_ramp_hold"]),
                           ("HM", HM_ramp, _AEM["hm_ramp_hold"]))
        size(R, 1) >= 3 || error("$tag ramp has only $(size(R,1)) points")
        issorted(R[:, 1]) || error("$tag ramp times are not increasing")
        R[1, 1] < 0 || error(
            "$tag ramp starts at $(R[1,1]) s, not before zero: the ON-TIME limb " *
            "is missing. Without it the forward is told the current had been on " *
            "forever, which over-predicted LM gate 26 by 52% on this floodplain " *
            "and is invisible in resistive ground.")
        maximum(R[:, 2]) ≈ 1.0 || error(
            "$tag ramp peak amplitude is $(maximum(R[:,2])), not 1.0: the " *
            "waveform is normalised in the .des and this one is not")
        R[end, 1] ≈ hold || error(
            "$tag ramp ends at $(R[end,1]) s but [aem] declares a hold of $hold s")
        # THE SWITCH-OFF TIME. A reviewer doubled LM's and the check passed, so
        # "both waveforms are checked" was true only of start, peak and end.
        #
        # The invariant is physical, not a range anyone picked: the current must
        # reach zero BEFORE the first gate opens, which is why the gates start
        # where they do. Measured on the shipped delivery, LM switches off in
        # 12.2 us against a first gate at 14.6 us, and HM in 304 us against
        # 427 us. Doubling either breaks it. A first attempt at this check
        # asserted "1 to 100 us" from assumption and refused the shipped HM
        # ramp on load -- the same mistake as every other claim in this file's
        # history, which is asserting instead of measuring.
        let ipk = findlast(x -> x >= 0.999, R[:, 2]),
            iz  = findfirst(j -> R[j, 2] <= 1e-12, ipk:size(R, 1))
            iz === nothing && error("$tag ramp has no zero after the peak")
            t_off = R[ipk + iz - 1, 1] - R[ipk, 1]
            t_off > 0 || error("$tag switch-off is not positive: $(t_off) s")
            g1 = (tag == "LM" ? LM_WINDOW_US[1, 1] : HM_WINDOW_US[1, 1]) * 1e-6
            t_off < g1 || error(
                "$tag switches off $(t_off) s after the peak but its first gate " *
                "opens at $(g1) s. The current must reach zero before the first " *
                "gate opens -- the .des waveform or the gate table has been " *
                "altered, and the forward would integrate current the gates " *
                "are supposed to be measuring after.")
        end
    end
    _AEM["lm_ramp_hold"] > 0 && _AEM["hm_ramp_hold"] > _AEM["lm_ramp_hold"] || error(
        "[aem] ramp holds look wrong: lm $(_AEM["lm_ramp_hold"]), hm $(_AEM["hm_ramp_hold"])")
    return (lm = length(LM_GATES), hm = length(HM_GATES))
end

# CALL IT ON LOAD. A reviewer grepped the whole branch and found
# check_against_config defined once and called nowhere -- not on load, not by
# FieldConfig, not by any test -- while the header and the README both said the
# gate table "is cross-checked". It was cross-checked only if a person thought
# to run the function. Parsing the .des at load time is pointless if the
# consistency of what was parsed is optional.
check_against_config()

end # module
