# xcite_system.jl -- the NRG XCite system, as delivered (Birrindudu, GA job 5031).
#
# The counterpart of udf_system.jl for a central-loop, single-moment system,
# inverted with HiQGA's VTEM1DInversion.dBzdt. Selected by [aem] system = "xcite".
# It exports the same API as UDFSystem (build_grid, build_operator, predicted,
# data_blocks, gate_times, describe, MU0, PV, rTx), so the AEM stages are
# written once.
#
# PARSED FROM THE .des AT LOAD TIME ([paths] aem_system * ".des"):
#   * the gate table ("Processed Data Window Times": gate, start, end, middle,
#     width in ms, measured from the START of the current pulse), and
#   * the Tx current waveform, column 5 ("Normalised Integrated Primary Field
#     Calibration Waveform").
#   Both are checked against themselves (middle = mean of start and end,
#   width = end - start) and the loop area is checked against the config.
#
# DECLARED IN THE CONFIG: tx_area_m2, rx_dz_m, lowpass_hz, nominal_unc,
# halt_noise_pv. The choices, and why (fieldrun/README.md, "AEM systems",
# has the summary):
#
#   WAVEFORM = column 5 + two fitted first-order low-pass corners, NOT column 6.
#   Column 6 is the same current already convolved with the Rx coil and
#   filters, which is the physically right signal path, but it has not decayed
#   when the first gates open (0.8% of peak at gate 0) and HiQGA's ramp
#   convolution does not converge for gates inside a ramp: gates 0-7 moved by up
#   to 65 noise sd with integration settings alone. Column 5 reaches zero at
#   5.4528 ms, before the first gate opens at 5.5136 ms, and the Rx path is then
#   applied analytically by `lowpassfcs` -- converged to 0.16 sd. The corners
#   ([aem] lowpass_hz, in the control file) were FITTED so that column 5 through
#   them reproduces column 6 to 0.03% rms of peak. The contractor's real filter
#   parameters are commercial-in-confidence, so this is an effective model, not
#   the filter design: a review item.
#
#   GATE TIMES = the .des arithmetic middles, on the same clock as the ramp, so
#   no shift is applied. They match the window_times_ms.txt of GA's own HiQGA
#   run exactly. HiQGA evaluates each gate at its centre (no window
#   integration). Gates 4-8 OVERLAP in the delivered table; that is the
#   contractor's windowing, kept and flagged, and is why the non-overlap check
#   udf_system.jl makes is deliberately absent here.
#
#   PITCH AND ROLL are not modelled (HiQGA's VTEM operator assumes a level
#   loop); prep_aem_xcite.jl refuses soundings tilted more than [aem]
#   max_tilt_deg instead.

module XCiteSystem

using HiQGA.transD_GP
using Statistics
using ..FieldConfig

export SYSTEM, TIMES, ramp, rTx, dzRx, lowpassfcs, GATES, WINDOW_MS, DES,
       build_grid, build_operator, predicted, data_blocks, gate_times,
       MU0, PV, describe, check_against_config

const SYSTEM = "xcite"
const MU0 = 4π * 1e-7
const PV  = 1e-12                 # delivered units: pV/(A.m^4)

const DES = let c = FieldConfig.CFG[]
    c === nothing && error("xcite_system.jl needs the control file: include " *
                           "common/fieldconfig.jl and call load_config first")
    c["paths"]["aem_system"] * ".des"
end
const _AEM = FieldConfig.CFG[]["aem"]

# Douglas-Peucker tolerance for the ramp, in units of normalised current.
# Measured (2026-11-12, 6 test earths): 1e-3 keeps 45 of 853 knots with a max
# forward error of 0.036 noise sd and cuts the forward from 11.6 to 2.8 ms.
const RAMP_TOL = 1e-3

isnum(s) = tryparse(Float64, s) !== nothing

"""
The gate table. Rows under "Processed Data Window Times" are
`gate start end middle width`, in ms from the start of the current pulse.
"""
function parse_gates()
    L = readlines(DES)
    i = findfirst(l -> occursin("Processed Data Window Times", l), L)
    i === nothing && error("no 'Processed Data Window Times' table in $DES")
    n, W, C, Wd = Int[], Matrix{Float64}(undef, 0, 2), Float64[], Float64[]
    for l in L[i+1:end]
        f = split(strip(l))
        if length(f) == 5 && all(isnum, f)
            push!(n, parse(Int, f[1]))
            W = vcat(W, [parse(Float64, f[2]) parse(Float64, f[3])])
            push!(C, parse(Float64, f[4])); push!(Wd, parse(Float64, f[5]))
        elseif !isempty(n) && !isempty(strip(l))
            break
        end
    end
    isempty(n) && error("gate table in $DES parsed empty")
    (n, W, C, Wd)
end

"""
The Tx current: column 5 of the "Transmitter Waveform" listing, sign-flipped
and normalised to peak 1, truncated at the first zero after the peak.
"""
function parse_current()
    L = readlines(DES)
    i = findfirst(l -> occursin("Transmitter Waveform", l), L)
    i === nothing && error("no 'Transmitter Waveform' listing in $DES")
    t, a = Float64[], Float64[]
    for l in L[i+1:end]
        f = split(strip(l))
        length(f) == 6 && all(isnum, f) || continue
        push!(t, parse(Float64, f[2])); push!(a, parse(Float64, f[5]))
    end
    isempty(t) && error("waveform listing in $DES parsed empty")
    pk = argmax(abs.(a)); I = a ./ a[pk]          # flips the sign, peak = +1
    iz = findfirst(j -> I[j] <= 0, pk:length(I))
    iz === nothing && error("the $DES current never returns to zero after its peak")
    iz += pk - 1
    I[iz] = 0.0
    (t[1:iz], I[1:iz], length(t))
end

"Douglas-Peucker: indices of the knots kept at tolerance `tol`."
function decimate(t, y, tol)
    keep = falses(length(t)); keep[1] = keep[end] = true
    st = [(1, length(t))]
    while !isempty(st)
        a, b = pop!(st); b - a < 2 && continue
        yl = y[a] .+ (y[b] - y[a]) .* (t[a+1:b-1] .- t[a]) ./ (t[b] - t[a])
        e, k = findmax(abs.(y[a+1:b-1] .- yl)); k += a
        e > tol && (keep[k] = true; push!(st, (a, k)); push!(st, (k, b)))
    end
    findall(keep)
end

"The loop area the .des states, m^2."
function des_loop_area()
    for l in readlines(DES)
        m = match(r"Loop Area:\s*([\d.]+)", l)
        m === nothing || return parse(Float64, m[1])
    end
    error("no 'Loop Area' in $DES")
end

const _G = parse_gates()
const GATES, WINDOW_MS = _G[1], _G[2]
const TIMES = _G[3] .* 1e-3                         # arithmetic middles, s
const _I = parse_current()
const ramp = let k = decimate(_I[1], _I[2], RAMP_TOL)
    [_I[1][k] _I[2][k]]
end
const N_RAMP_FULL = length(_I[1])

const rTx = sqrt(_AEM["tx_area_m2"] / π)            # equal-area radius, 9.20 m
const dzRx = Float64(_AEM["rx_dz_m"])               # z DOWN: -0.25 = Rx 0.25 m above Tx
const lowpassfcs = Float64.(_AEM["lowpass_hz"])

"Same log-growing grid as UDFSystem.build_grid, so the two surveys compare."
function build_grid(; dz = 1.5, extendfrac = 1.03, n = 70)
    zall, znall, zbnd = transD_GP.setupz(0.0, extendfrac; dz = dz, n = n,
                                         showplot = false)
    z, ρ, nfixed = transD_GP.makezρ(zbnd; zfixed = [-1e5], ρfixed = [1e12])
    (zall = zall, znall = znall, zbnd = zbnd, z = z, ρ = ρ, nfixed = nfixed)
end

"""
    build_operator(d, σ, height; grid, calcjacobian, σ_scale, ρ_start)

`VTEM1DInversion.dBzdt` for one sounding. `d`/`σ` are dBdt_Z_F in pV/(A.m^4);
`height` is Tx_Loop_Clearance. The receiver is concentric, `dzRx` above.
"""
function build_operator(d::AbstractVector, σ::AbstractVector, height::Real;
                        grid = build_grid(), calcjacobian = true,
                        σ_scale::Real = 1.0, ρ_start::Real = 30.0, useML::Bool = false)
    length(d) == length(TIMES) || error("$(length(d)) data for $(length(TIMES)) gates")
    zTx = -float(height)
    zRx = zTx + dzRx
    k = PV / MU0
    # the air layer keeps its 1e12 ohm.m; see udf_system.jl build_operator
    ρ0 = fill(float(ρ_start), length(grid.z))
    ρ0[1:grid.nfixed] .= grid.ρ[1:grid.nfixed]
    transD_GP.VTEM1DInversion.dBzdt(;
        times = TIMES, ramp = ramp, rTx = rTx, zTx = zTx, zRx = zRx,
        z = grid.z, ρ = ρ0, nfixed = grid.nfixed,
        d = d .* k, σ = σ .* (σ_scale * k),
        lowpassfcs = lowpassfcs, calcjacobian = calcjacobian, useML = useML)
end

# ── the system-independent API the stages call ─────────────────────────────
build_operator(g::NamedTuple, height::Real; kw...) = build_operator(g.d, g.σ, height; kw...)
"Forward response at the operator's current model, pV/(A.m^4), in data order."
predicted(aem) = aem.F.dBzdt[aem.select] .* (MU0 / PV)
"Named index ranges of the data vector, for per-block residual summaries."
data_blocks(g) = [("Z", 1:length(g.d))]
gate_times() = (Z_times = TIMES,)

function describe()
    println("NRG XCite (central loop, Z), parsed from $(basename(DES))")
    println("  gates $(first(GATES))-$(last(GATES)) ($(length(TIMES))), ",
            "$(round(1e3TIMES[1], digits=4))-$(round(1e3TIMES[end], digits=4)) ms from pulse start")
    println("  ramp: column 5, $(size(ramp,1)) of $N_RAMP_FULL knots (tol $RAMP_TOL), ",
            "off at $(round(1e3ramp[end,1], digits=4)) ms")
    println("  rTx $(round(rTx, digits=2)) m, dzRx $dzRx m, lowpass $(lowpassfcs) Hz (fitted, provisional)")
end

"""
    check_against_config()

The .des against itself and against the control file. Called on load.
"""
function check_against_config()
    N, W, C, Wd = _G
    TOL = 1e-9
    # 1. centre is the mean of start and end; 2. width is end - start. The
    #    table is printed to 4-5 decimals in ms, so compare at that precision.
    for i in eachindex(N)
        abs((W[i,1] + W[i,2]) / 2 - C[i]) <= 5e-5 || error(
            "gate $(N[i]): start $(W[i,1]) and end $(W[i,2]) ms give a middle of " *
            "$((W[i,1]+W[i,2])/2) but the .des lists $(C[i]): a time column has been altered")
        abs((W[i,2] - W[i,1]) - Wd[i]) <= 5e-5 || error(
            "gate $(N[i]): end - start is $(W[i,2]-W[i,1]) ms but the .des lists a width of $(Wd[i])")
    end
    all(W[:,1] .< W[:,2]) || error("a gate's start is not before its end")
    issorted(C) || error("gate middles are not increasing")
    N == collect(first(N):last(N)) || error("gate numbers are not consecutive and unique: $N")
    if haskey(_AEM, "gates")
        lo, hi = _AEM["gates"]
        (first(N) == lo && last(N) == hi) || error(
            "[aem] gates declares $lo-$hi but the .des table runs $(first(N))-$(last(N))")
    end
    # 3. the noise floor has one value per gate
    length(_AEM["halt_noise_pv"]) == length(N) || error(
        "[aem] halt_noise_pv has $(length(_AEM["halt_noise_pv"])) values for $(length(N)) gates")
    # 4. the loop area the config declares is the one the .des states
    a = des_loop_area()
    abs(a - _AEM["tx_area_m2"]) <= 1e-6 * a || error(
        "[aem] tx_area_m2 = $(_AEM["tx_area_m2"]) but the .des states a loop area of $a m^2")
    # 5. the ramp: starts at zero current, peaks at 1, and is OFF before the
    #    first gate opens -- the property that makes column 5 converge where
    #    column 6 does not
    issorted(ramp[:,1]) || error("ramp times are not increasing")
    abs(ramp[1,2]) <= 1e-12 || error("ramp does not start at zero current")
    maximum(ramp[:,2]) ≈ 1.0 || error("ramp peak is $(maximum(ramp[:,2])), not 1")
    ramp[end,2] == 0 || error("ramp does not end at zero current")
    ramp[end,1] < WINDOW_MS[1,1] * 1e-3 || error(
        "the current is still on at $(ramp[end,1]) s but gate $(first(N)) opens at " *
        "$(WINDOW_MS[1,1]*1e-3) s; HiQGA's ramp convolution is not converged for " *
        "gates inside the ramp")
    all(>(0), lowpassfcs) || error("[aem] lowpass_hz must be positive")
    (gates = length(N), ramp_knots = size(ramp, 1))
end

check_against_config()

end # module
