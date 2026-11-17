# smr_data.jl -- forward-model the SNMR data Stage 2 inverts.
#
#   julia --project=. tools/smr_data.jl --tag <TAG>
#
# FIXED 18 Nov 2026. THE GENERATING KERNEL IS BUILT ON THE TRUE CONDUCTIVITY.
#
# Until today this script took the true water content from the truth file -- which
# was right -- and then built the generating kernel on an AEM INVERSION's sigma(z),
# either the Occam profile or a Stage 1 trans-D p50 selected by --label. The
# measurements Stage 2 inverts were therefore forward-modelled through an
# inversion result. Measured against the true sigma_bulk at the same depths, the
# conductivity used departed by a median 55% (3P), 54% (2P), 45% (3S) and 52%
# (2S), with maxima of 400 to 550%, and was bit-identical to aem_occam<TAG>.jld2.
# The SMR kernel is conductivity-sensitive, so this was not cosmetic, and it hit
# hardest on 2P, the column whose aquifer is a resistor enclosed by conductors.
#
# Now: sigma for the generating kernel is T["sigma_bulk"] taken onto the KERNEL
# DEPTH GRID -- the same 65-layer log-growing grid src/aem_setup.jl defines and
# Stage 0 generates the AEM data on -- by the same piecewise_interp Stage 0 uses.
# The grid comes from the grid DEFINITION, not from any inversion output, and
# NOTHING FROM STAGE 1 ENTERS DATA GENERATION. Both datasets therefore see the
# same discretised ground.
#
# The truth's own 0.25 m grid was tried first and is not usable: SNMRForward's
# field computation overflows with 800 thin layers and returns amplitudes around
# 1e63 nV instead of tens of nV. Some discretisation is forced, so it is the same
# one the AEM data already carry.
# The data no longer depend on --label or on occam/transd, so there is ONE data
# file per column; both switches are now hard errors rather than being quietly
# accepted. The generating conductivity is stored as sigma_gen and asserted
# against the truth on write.
#
# The file name keeps its "T" suffix, which Stage 2 expects (it is invoked with
# --tag <TAG>T). The suffix no longer denotes a kernel source -- there is only
# one -- and is retained solely so Stage 2 needs no change.
#
# THIS SCRIPT WAS MISSING FROM THE RELEASE. Stage 2 loads
# data/smr_data_murray<TAG>.jld2 and nothing shipped wrote it, so an
# independent reviewer got through Stage 1 and then hit a wall with no way
# forward and no script to find. The chain is:
#
#   tools/aem_HYP.jl   -> data/aem<TAG><LABEL>/aem_posterior.jld2
#   tools/smr_data.jl  -> data/smr_data_murray<TAG>T.jld2   (transd branch)
#   stages/stage2_snmr.jl
#
# Note the T: the trans-D branch writes a "T" suffix so the Occam and trans-D
# kernels can never collide, which means Stage 2's --tag is <TAG>T, not <TAG>.
#
#
# Kernel is built from the AEM **p50** σ(z), because that is what the EK pre-run
# did and what condition C3 assumes. C4 rebuilds the kernel per AEM draw; that
# happens in run_conditions.jl, not here.
#
# Noise is additive and absolute, in two tiers taken from the EK pre-run's
# recovered noise population (effective √s·σ̂): "typical" = the survey median,
# "poor" = the p90. The per-q shape rises with pulse moment, as measured at EK
# (RESULTS_PRERUN §3: CV ≈ 0.60, rising ~2.5× across the moments).
#
# Amplitudes cannot be negative, so the noise draw is a TRUNCATED normal rather
# than a plain one — at the poor tier the deep, low-amplitude moments would
# otherwise take negative values and the amplitude-only likelihood would see
# data the forward cannot produce.
#
#   julia --project=ek_prerun synthetics/smr_data_murray.jl

using Random, Statistics, Printf, JLD2
using SNMRForward, SMRPInversion

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))   # tools/ -> repository root
const OUT  = joinpath(ROOT, "data")
include(joinpath(ROOT, "src", "prerun_config.jl")); const CFG = PrerunConfig
# For build_aem_grid and piecewise_interp: the kernel depth grid and the mapping
# onto it, from the one place that defines them. Loads cleanly beside
# SNMRForward; verified 18 Nov.
include(joinpath(ROOT, "src", "aem_setup.jl"))

# EFTF q-set, as synthetics/common.jl defines it.
const QGRID = collect(10.0 .^ range(log10(0.1), log10(15.0), length = 24))

# Effective noise from the EK pre-run's recovered s population, √s·σ̂ in nV.
# --aem transd uses the trans-D p50 instead of the Occam profile, and writes
# smr_data_murray2T.jld2 so the two never collide.
# Both switches chose which INVERSION generated the data. They are gone, and
# passing either is an error: an old invocation must fail loudly, not silently
# produce a file that looks the same and means something else.
"transd" in ARGS && error("""
    `transd` is no longer accepted. The generating kernel is built on the TRUE
    conductivity, so there is no occam/transd choice to make and one data file
    per column. Run: julia --project=. tools/smr_data.jl --tag <TAG>""")
findfirst(==("--label"), ARGS) === nothing || error("""
    `--label` is no longer accepted. It selected which Stage 1 ensemble's p50
    GENERATED the data, which was the fault fixed on 18 Nov. Nothing from Stage 1
    enters data generation now. Run: julia --project=. tools/smr_data.jl --tag <TAG>""")
# aem_HYP.jl's 4th argument writes data/aem<TAG><LABEL>/, so a labelled Stage 1
# run must be named here too or this looks in the wrong directory -- the same
# trap aem_draws.jl had.

# TAG selects the truth column: 2 = two-aquifer Murray, 3 = unconfined resistive.
# The suffix on the SMR data file is <TAG> for an Occam kernel and <TAG>T for a
# trans-D one, so the two kernels never collide for the same column.
const TAG = let i = findfirst(==("--tag"), ARGS)
    i === nothing ? "2" : ARGS[i+1]
end
# Retained so Stage 2, which is invoked with --tag <TAG>T, finds the file.
const SUF = TAG*"T"
const TIERS = Dict("typical" => 4.63e-9, "poor" => 11.73e-9)
const SHAPE_RISE = 2.5     # σ̂ rises this much across the moments (EK-measured)
const SHAPE_CV   = 0.60    # and has this CV overall (EK-measured)
const SEED       = 20261012
# --noiseseed varies the NOISE DRAW only. The per-q noise SHAPE (sigma_shape) is
# a survey property measured at EK, so it stays on SEED for every realisation;
# only add_noise is re-drawn. Without this every result in the project rests on
# one noise realisation, which makes a statement like "the truth sits at
# posterior quantile 0.55" an anecdote rather than a calibration.
const NOISESEED = let i = findfirst(==("--noiseseed"), ARGS)
    i === nothing ? SEED : parse(Int, ARGS[i+1])
end

"""
Per-q noise standard deviations: a monotone rise across the moments times
lognormal scatter, scaled so the MEDIAN equals `level`.

A pure linear ramp cannot reach CV 0.60 — its limit as the ratio grows is
2/√12 ≈ 0.577 — so the measured CV needs scatter on top of the rise. The
scatter width is solved for numerically rather than guessed.
"""
function sigma_shape(nq, level; rise = SHAPE_RISE, cv = SHAPE_CV, seed = SEED)
    ramp = range(1.0, rise, length = nq) |> collect
    f(w) = begin
        Random.seed!(seed)
        s = ramp .* exp.(w .* randn(nq))
        std(s) / mean(s)
    end
    lo, hi = 0.0, 2.0
    for _ in 1:60                       # bisection on the scatter width
        mid = (lo + hi) / 2
        f(mid) < cv ? (lo = mid) : (hi = mid)
    end
    Random.seed!(seed)
    s = ramp .* exp.(((lo + hi) / 2) .* randn(nq))
    s .* (level / median(s))
end

"Truncated-normal noise: redraw any sample that would make an amplitude negative."
function add_noise(V, σ̂; seed = SEED, maxtries = 1000)
    Random.seed!(seed)
    out = similar(V)
    for j in eachindex(V)
        for t in 1:maxtries
            v = V[j] + σ̂[j] * randn()
            if v > 0
                out[j] = v; break
            end
            t == maxtries && (out[j] = V[j] * 1e-3)   # pathological; logged below
        end
    end
    out
end

"AEM σ on the SMR grid, interpolated by LAYER BOUNDARIES (not midpoints)."
function to_smr_grid(σ_aem, z_aem, zsmr)
    [σ_aem[clamp(searchsortedlast(z_aem, zs), 1, length(σ_aem))] for zs in zsmr]
end

function main()
    T = load(joinpath(OUT, "truth_murray$(TAG).jld2"))

    zsmr  = T["zsmr"]
    w     = T["θ_smrgrid"]                    # truth θ_SMR on the SMR grid
    Be, ϕ = T["BE_T"], deg2rad(T["INC_DEG"])
    R     = T["LOOP_SIDE"] / 2                # MRSForward_square takes the HALF-side

    # Conductivity model for the generating kernel: the TRUE sigma_bulk, on the
    # truth's own grid. No interpolation, so no discretisation is introduced and
    # the kernel sees the ground the water content came from. Nothing from any
    # inversion, and nothing from Stage 1, enters here.
    #
    # ConductivityModel takes n conductivities and n-1 thicknesses, the last
    # layer being the half-space, which is why t_gen is a diff.
    zgen = collect(Float64, build_aem_grid().zall)
    σgen = piecewise_interp(zgen, collect(Float64, T["z"]), collect(Float64, T["σ_bulk"]))
    tgen = diff(zgen)
    length(σgen) == length(tgen) + 1 ||
        error("sigma_gen has $(length(σgen)) layers against $(length(tgen)) thicknesses")

    @printf("kernel from the TRUE sigma_bulk on the kernel grid: %d layers, σ %.4f–%.4f S/m, median %.4f\n",
            length(σgen), extrema(σgen)..., median(σgen))
    @printf("AGRF Be = %.0f nT, inc = %.1f°, loop %.1f m square, %d moments\n",
            Be*1e9, T["INC_DEG"], T["LOOP_SIDE"], length(QGRID))

    ct = SNMRForward.ConductivityModel(σgen, tgen)
    F  = SNMRForward.MRSForward_square(R, zsmr, QGRID, ϕ, 0.0, Be, ct)
    # SNMRForward.forward does trapezoidal integration over depth with F.dz and
    # the volume element already folded into the stored kernel — it is NOT a
    # plain kernel*w matrix multiply, and treating it as one silently drops the
    # depth weighting. Amplitude-only means abs of the complex response.
    V0 = abs.(SNMRForward.forward(F, w))

    @printf("forward V0: %.2f–%.2f nV\n", 1e9*minimum(V0), 1e9*maximum(V0))

    data = Dict{String,Any}()
    for (tier, level) in TIERS
        σ̂ = sigma_shape(length(QGRID), level; seed = SEED)
        Vobs = add_noise(V0, σ̂; seed = NOISESEED + (tier == "poor" ? 7 : 0))
        snr = V0 ./ σ̂
        @printf("\ntier %-8s σ̂ median %.2f nV (CV %.2f, maxmin %.1f×)  SNR %.1f–%.1f, %d moments < 2\n",
                tier, 1e9*median(σ̂), std(σ̂)/mean(σ̂), maximum(σ̂)/minimum(σ̂),
                minimum(snr), maximum(snr), count(<(2), snr))
        data[tier] = (V0_noisy = Vobs, σ̂ = σ̂, snr = snr)
    end

    # ASSERT BEFORE THE WRITE that what generated the data is the truth and
    # nothing else. This is the check whose absence let an inversion's sigma
    # generate the measurements for months.
    σcheck = piecewise_interp(zgen, collect(Float64, T["z"]), collect(Float64, T["σ_bulk"]))
    σgen == σcheck || error("""
        sigma_gen is not the truth's sigma_bulk on the kernel grid.
        max |difference| = $(maximum(abs.(σgen .- σcheck))) S/m over $(length(σgen)) layers.
        The generating kernel must be built on the truth. Refusing to write.""")
    # and that it is nowhere near any inversion's sigma, which is the failure mode
    aeminv = joinpath(OUT, "aem_occam$(TAG).jld2")
    if isfile(aeminv)
        so = load(aeminv, "σ_occam")
        length(so) >= length(σgen) && maximum(abs.(σgen .- so[1:length(σgen)])) < 1e-12 &&
            error("sigma_gen is bit-identical to $aeminv -- that is the 18 Nov fault. Refusing to write.")
    end
    @printf("assert: sigma_gen is the truth's sigma_bulk on the kernel grid over %d layers, exactly\n",
            length(σgen))
    @printf("        and is not the Occam inversion's sigma\n")

    jldsave(joinpath(OUT, "smr_data_murray$(SUF).jld2");
        zsmr = zsmr, w_true = w, qgrid = QGRID, V0_clean = V0,
        kernel = Matrix(F.kernel), kernel_dz = collect(F.dz),
        sigma_gen = σgen, z_gen = zgen, t_gen = tgen,
        generating_sigma_source = "truth sigma_bulk, truth grid, no interpolation",
        σ_kernel = σgen, z_kernel = zgen, t_kernel = tgen,
        Be = Be, inc_deg = T["INC_DEG"], loop_side = T["LOOP_SIDE"], R = R,
        typical_V = data["typical"].V0_noisy, typical_σ̂ = data["typical"].σ̂,
        poor_V    = data["poor"].V0_noisy,    poor_σ̂    = data["poor"].σ̂,
        tiers = TIERS, seed = SEED, noiseseed = NOISESEED,
        shape_rise = SHAPE_RISE, shape_cv = SHAPE_CV)
    println("\nwrote ", joinpath(OUT, "smr_data_murray$(SUF).jld2"))
end

main()
