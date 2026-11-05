# smr_data.jl -- forward-model the SNMR data Stage 2 inverts.
#
#   julia --project=. tools/smr_data.jl transd --tag <TAG> [--label <L>]
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

# EFTF q-set, as synthetics/common.jl defines it.
const QGRID = collect(10.0 .^ range(log10(0.1), log10(15.0), length = 24))

# Effective noise from the EK pre-run's recovered s population, √s·σ̂ in nV.
# --aem transd uses the trans-D p50 instead of the Occam profile, and writes
# smr_data_murray2T.jld2 so the two never collide.
const USE_TRANSD = "transd" in ARGS
# aem_HYP.jl's 4th argument writes data/aem<TAG><LABEL>/, so a labelled Stage 1
# run must be named here too or this looks in the wrong directory -- the same
# trap aem_draws.jl had.
const LABEL = let i = findfirst(==("--label"), ARGS)
    i === nothing ? "" : ARGS[i+1]
end
# TAG selects the truth column: 2 = two-aquifer Murray, 3 = unconfined resistive.
# The suffix on the SMR data file is <TAG> for an Occam kernel and <TAG>T for a
# trans-D one, so the two kernels never collide for the same column.
const TAG = let i = findfirst(==("--tag"), ARGS)
    i === nothing ? "2" : ARGS[i+1]
end
const SUF = USE_TRANSD ? TAG*"T" : TAG
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
    af = USE_TRANSD ? joinpath(OUT, "aem$(TAG)$(LABEL)", "aem_posterior.jld2") :
                      joinpath(OUT, TAG == "2" ? "aem_occam.jld2" : "aem_occam$(TAG).jld2")
    isfile(af) || error("missing $af")
    A = load(af)

    zsmr  = T["zsmr"]
    w     = T["θ_smrgrid"]                    # truth θ_SMR on the SMR grid
    Be, ϕ = T["BE_T"], deg2rad(T["INC_DEG"])
    R     = T["LOOP_SIDE"] / 2                # MRSForward_square takes the HALF-side

    # Conductivity model for the kernel: the AEM p50, on the AEM grid's own
    # layering. Interpolating it onto the coarse SMR grid first would throw away
    # resolution the kernel can use.
    # Occam gives ONE sigma(z), no posterior. The kernel is built from it and
    # every sigma-dependent predicate conditions on it as a point estimate.
    zA = A["zall_aem"]
    keep = findall(z -> z > 0, zA)
    σ50  = (USE_TRANSD ? A["σ50"] : A["σ_occam"])[keep]; zA = zA[keep]
    tA   = diff(zA)

    @printf("kernel from AEM %s: %d layers, σ %.4f–%.4f S/m\n",
            USE_TRANSD ? "TRANS-D p50" : "OCCAM", length(σ50), extrema(σ50)...)
    @printf("AGRF Be = %.0f nT, inc = %.1f°, loop %.1f m square, %d moments\n",
            Be*1e9, T["INC_DEG"], T["LOOP_SIDE"], length(QGRID))

    ct = SNMRForward.ConductivityModel(σ50, tA)
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

    jldsave(joinpath(OUT, "smr_data_murray$(SUF).jld2");
        zsmr = zsmr, w_true = w, qgrid = QGRID, V0_clean = V0,
        kernel = Matrix(F.kernel), kernel_dz = collect(F.dz),
        σ_kernel = σ50, z_kernel = zA, t_kernel = tA,
        Be = Be, inc_deg = T["INC_DEG"], loop_side = T["LOOP_SIDE"], R = R,
        typical_V = data["typical"].V0_noisy, typical_σ̂ = data["typical"].σ̂,
        poor_V    = data["poor"].V0_noisy,    poor_σ̂    = data["poor"].σ̂,
        tiers = TIERS, seed = SEED, noiseseed = NOISESEED,
        shape_rise = SHAPE_RISE, shape_cv = SHAPE_CV)
    println("\nwrote ", joinpath(OUT, "smr_data_murray$(SUF).jld2"))
end

main()
