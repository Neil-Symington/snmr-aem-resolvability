# quicklook.jl — do the two inversions run, and do they fit? Nothing else.
#
#   julia --project=. tools/quicklook.jl <TAG>
#   julia --project=. tools/quicklook.jl 2P --have-smrdata    # reuse smr_data
#
# WHAT THIS IS FOR. Stages 1 and 2 are trans-dimensional MCMC and take hours, so
# the first thing you want from a newly built column is cheaper than that: did
# the AEM inversion converge and fit its 42 gates, did the SNMR inversion
# converge and fit its pulse moments, and do the two recovered profiles look
# like the truth they were built from. That is a check on the inversions. It
# runs in a few minutes because an Occam solve is a handful of Gauss-Newton
# steps and needs no workers.
#
# WHAT THIS IS NOT. It is not a small version of the real answer, and no number
# it prints may be quoted as a result.
#
#   * There is no posterior here. Occam gives ONE model, and the linearised
#     covariance around it is not a credible interval: at East Kimberley
#     OK18_38 the published inversion of the same sounding sat outside the
#     bootstrap 5-95 band at 38 of 41 depths while the band claimed +/-4%.
#   * It cannot do detection. detect_aquifer.jl takes the smaller of
#     detectability and posterior-responding, and the responding criterion
#     needs a PRIOR width that a deterministic solve does not have. Run without
#     it, zlim fell to the bottom of the parameterisation and the selected
#     screen pinned at 137-147 m with 37 candidate windows within 5% of the
#     best. So: no W, no P(W >= w), no screen, no gate, no verdict.
#   * The lambda-squared it lands on is a damping choice, not a prior.
#
# It also closes a producer gap. tools/smr_data.jl and stages/stage2_snmr.jl
# both READ data/aem_occam<TAG>.jld2, and nothing in the repository wrote it --
# which is why [inversion.aem] engine = "occam" had no producer. This writes it,
# so an Occam kernel is reachable on the synthetic side.
using JLD2, Statistics, Printf, LinearAlgebra
using SNMRForward, SMRPInversion

const HERE = normpath(@__DIR__)
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = joinpath(ROOT, "data")
const TAG  = let a = filter(x -> !startswith(x, "--"), ARGS)
    isempty(a) && error("usage: quicklook.jl <TAG> [--have-smrdata]\n" *
                        "  TAG names data/truth_murray<TAG>.jld2")
    a[1]
end
const HAVE_SMRDATA = "--have-smrdata" in ARGS

# Pinned to HiQGA's makenoisydata! default, as every other synthetic stage is:
# the caller's seed was silently ignored until build_skytem_operator passed
# rseed through, so 123 is the realisation behind every published number.
const NOISE_RSEED = 123

# The Occam convention this project uses, from src/prerun_config.jl §5. beta2
# is NOT zero: a difference operator annihilates a constant reference, and
# SMRPInversion applies wref only when beta2 > 0, so at zero the reference does
# nothing at all.
const SMR_W0     = 0.02
const SMR_REG    = :R2
const SMR_LAM2   = (-1.0, 6.0)
const SMR_BETA2  = 0.01

const GADI_SRC = joinpath(ROOT, "src")
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "occam_core.jl")); using .OccamCore
using HiQGA.transD_GP
# aem_jac.jl, NOT build_skytem_operator. The shared builder does not set
# calcjacobian, so the Gauss-Newton steps in occam_core.jl cannot use it: the
# first version of this file did, and the forward model died with a BoundsError
# deep inside AEM_VMD_HMD rather than saying the operator was wrong. Same
# builder and the same seed as test/test_gate_solver.jl, so the chi2 here is
# directly comparable with the numbers that suite reports.
include(joinpath(ROOT, "src", "aem_jac.jl"))

rms(a, b) = sqrt(mean(abs2, a .- b))
band(z, lo, hi) = findall(zz -> lo <= zz <= hi, z)

"A chi2/N far from 1 means the fit failed, in one direction or the other."
function verdict(name, χ2N)
    if χ2N > 5.0
        @printf("  FAIL  %s did not fit: chi2/N %.2f\n", name, χ2N); return false
    elseif χ2N < 0.2
        @printf("  FAIL  %s fits far better than the noise allows: chi2/N %.2f\n", name, χ2N)
        @printf("        that is over-fitting or a mis-stated noise level, not a good result\n")
        return false
    end
    true
end

# ── AEM ─────────────────────────────────────────────────────────────────────
function aem_quicklook(T)
    B = build_aem_jac(T, NOISE_RSEED)
    aem, grid = B.aem, B.grid
    n = length(aem.ρ) - 1                 # ρ[1] is AIR: overwrite it and the
    ndata = length(aem.res)               # SkyTEM forward returns exactly zero
    zall = grid.zall
    @assert length(zall) == n "grid.zall ($(length(zall))) does not match the $n model cells"

    R = makeR1(n); m0 = fill(-2.0, n)     # model is log10 conductivity
    @printf("AEM Occam        %d gates, %d cells to %.1f m\n", ndata, n, zall[end])

    t0 = time()
    # :minimise first, as stage4_gate_v2.jl does, to see the best fit the
    # forward model can reach at all. Then :target from there for the Occam
    # answer: the smoothest model that still reaches chi2 <= N.
    mmin, χmin = occam(aem, m0, m0, R, ndata; mode = :minimise)
    mtar, χtar = occam(aem, mmin, m0, R, ndata; mode = :target)
    secs = time() - t0

    @printf("  best achievable   chi2 %8.2f   chi2/N %6.3f   (:minimise)\n", χmin, χmin/ndata)
    @printf("  Occam answer      chi2 %8.2f   chi2/N %6.3f   (:target)\n",   χtar, χtar/ndata)
    ok = verdict("the AEM inversion", χtar/ndata)

    σ    = 10 .^ mtar
    σtru = B.σt[2:end]                    # [1] is the air layer, set to 1e-12
    @assert length(σtru) == n
    d = log10.(σ) .- log10.(σtru)
    @printf("  recovered sigma vs truth, log10 RMS: %.3f dex over all %d cells\n", rms(log10.(σ), log10.(σtru)), n)
    for (lo, hi) in ((0.0, 50.0), (50.0, 100.0), (100.0, zall[end]))
        k = band(zall, lo, hi)
        isempty(k) && continue
        @printf("      %5.0f-%-5.0f m  %.3f dex   (bias %+.3f)\n", lo, hi, rms(d[k], 0 .* d[k]), mean(d[k]))
    end
    @printf("  %.1f s\n", secs)

    f = joinpath(OUT, TAG == "2" ? "aem_occam.jld2" : "aem_occam$(TAG).jld2")
    jldsave(f; σ_occam = σ, zall_aem = zall, zbnd_aem = grid.zbnd,
            σ_aem_truth = B.σt, chi2 = χtar, chi2_minimise = χmin,
            ndata = ndata, mode = "target", rseed_noise = NOISE_RSEED, tag = TAG,
            # Named so nothing downstream can mistake this for a posterior.
            note = "deterministic Occam, one model, no posterior; quicklook.jl")
    println("  wrote ", relpath(f, ROOT))
    ok, σ, zall
end

# ── SNMR ────────────────────────────────────────────────────────────────────
function smr_quicklook()
    df = joinpath(OUT, "smr_data_murray$(TAG).jld2")
    if !isfile(df) || !HAVE_SMRDATA
        println("\nbuilding the SNMR data on the Occam kernel")
        run(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(HERE, "smr_data.jl")) --tag $TAG`)
    end
    isfile(df) || error("smr_data.jl did not write $df")
    D = load(df)
    zsmr = D["zsmr"]; V = D["typical_V"]; σ̂ = D["typical_σ̂"]
    wtru = D["w_true"]; qgrid = D["qgrid"]
    Be = D["Be"]; ϕ = deg2rad(D["inc_deg"]); Rloop = D["R"]
    σk = D["σ_kernel"]; zk = D["z_kernel"]

    # The same kernel construction stage2_snmr.jl uses, on the Occam sigma.
    ct = SNMRForward.ConductivityModel(σk, diff(zk))
    F  = SNMRForward.MRSForward_square(Rloop, zsmr, qgrid, ϕ, 0.0, Be, ct)

    nz = length(zsmr)
    @printf("\nSNMR Occam       %d pulse moments, %d cells to %.1f m\n", length(V), nz, zsmr[end])
    # Known noise, not hierarchical: this is the deterministic branch, so sigma
    # is taken as stated rather than sampled. amponly means the phase is unused,
    # but the constructor still wants a sigma of matching length for it.
    S = SMRPInversion.newSMRSounding(V, zeros(length(V)), F;
                                     σ_V0 = σ̂, σ_ϕ = ones(length(V)))
    t0 = time()
    r = SMRPInversion.deterministic_invert(S;
            wstart = fill(SMR_W0, nz), wref = fill(SMR_W0, nz),
            regtype = SMR_REG, λ²min = SMR_LAM2[1], λ²max = SMR_LAM2[2],
            β² = SMR_BETA2, uncertainty = :none)
    secs = time() - t0
    nq = length(V)
    # gradientinv returns the whole path, not scalars: one entry per iteration,
    # each a vector over the lambda-squared ladder, with oidx naming the one it
    # kept. The answer is the last iteration's selection. Printing r.χ² itself
    # throws, which is how this was found.
    it = lastindex(r.oidx); j = r.oidx[it]
    w  = r.mnew[it][j]
    χ2 = r.χ²[it][j]
    λ2 = let v = r.λ²[it][j]; v isa AbstractVector ? first(v) : v end
    @printf("  chi2 %8.2f   chi2/N %6.3f   lambda2 %.3g   (%d iterations)\n",
            χ2, χ2/nq, λ2, it)
    ok = verdict("the SNMR inversion", χ2/nq)

    # Save it. This used to print the profile and keep nothing, which left the
    # deterministic SNMR branch with no producer -- the same gap aem_occam had
    # before this tool wrote it. Stages 3 and beyond need the profile, not the
    # RMS, so it goes to disk under its own name and cannot collide with a
    # cond_* ensemble file.
    sf = joinpath(OUT, "smr_occam$(TAG).jld2")
    jldsave(sf; θ_occam = w, zsmr = zsmr, chi2 = χ2, ndata = nq, λ² = λ2,
            θ_true = wtru, σ_kernel = σk, z_kernel = zk,
            regtype = String(SMR_REG), β² = SMR_BETA2, w0 = SMR_W0,
            tag = TAG, mode = "target",
            note = "deterministic Occam, one model, no posterior; quicklook.jl")
    println("  wrote ", relpath(sf, ROOT))

    @printf("  recovered theta vs truth, RMS: %.4f v/v over %d cells\n", rms(w, wtru), nz)
    for (lo, hi) in ((0.0, 40.0), (40.0, 80.0), (80.0, zsmr[end]))
        k = band(zsmr, lo, hi)
        isempty(k) && continue
        @printf("      %5.0f-%-5.0f m  %.4f v/v   (bias %+.4f)\n", lo, hi,
                rms(w[k], wtru[k]), mean(w[k] .- wtru[k]))
    end
    it, iw = argmax(wtru), argmax(w)
    @printf("  peak theta: truth %.3f at %.1f m,  recovered %.3f at %.1f m\n",
            wtru[it], zsmr[it], w[iw], zsmr[iw])
    @printf("  %.1f s\n", secs)
    ok
end

function main()
    tf = joinpath(OUT, "truth_murray$(TAG).jld2")
    isfile(tf) || error("""
    no truth at $tf

      Build the column first:
        julia --project=. synthetics/scripts/build_truth.jl synthetics/config/<scenario>.toml
        cp synthetics/output/<name>_truth.jld2 $tf
      or run ./run_smoke.sh with STAGES="truth".
    """)
    println("="^72)
    @printf("quicklook — tag %s — DO THE INVERSIONS RUN AND FIT?\n", TAG)
    println("this checks the inversions only. no posterior, no detection, no verdict.")
    println("="^72)
    T = load(tf)
    okA, _, _ = aem_quicklook(T)
    okS = smr_quicklook()
    println("\n", "="^72)
    if okA && okS
        println("BOTH INVERSIONS RAN AND FIT. Not an inference result: for W, for")
        println("P(W >= w) or for any verdict, run Stages 1 and 2 and then run_site.jl.")
    else
        println("AT LEAST ONE INVERSION DID NOT FIT. See the FAIL lines above.")
        println("A production run on this column would be hours spent on the same problem.")
        exit(1)
    end
end

main()
