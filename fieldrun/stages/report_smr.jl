# report_smr.jl -- standard SMR inversion figures and tables, per site.
#
#   julia --project=. fieldrun/stages/report_smr.jl --config <toml> [--site S | --sites A,B]
#
# Reads what run_smr_transd.jl left on disk -- smr_transd_<site>.jld2 (the
# pooled posterior) and smr_mcmc/<site>/ (every ensemble's chains) -- and runs
# no MCMC, so it can be re-run at will (smrp_invert.jl's --replot). Output,
# under <outputs>/report/<site>/:
#
#   <site>_post.png             pooled water-content posterior, full depth
#   <site>_post_100m.png        the same, top 100 m
#   <site>_sounding.png         V0 and phase vs pulse moment; bars = sigma-hat
#                               and sigma-hat x sqrt(posterior-median s)
#   <site>_himage_results.csv   posterior histogram, one row per depth cell
#   <site>_CI_results.csv       depth_from, depth_to, p05, p50, p95
#   misfit/<site>_e<k>_misfit_chi2.png   per ensemble: chi2 and k traces
#   swarm/<site>_e<k>_swarm.png          per ensemble: model swarm + fits
#
# These follow staging/scripts/smrp_invert.jl, with three departures, each
# deliberate:
#   * POSTERIOR AND TABLES ARE POOLED over the [kernel] n_draws ensembles: that
#     pooled posterior is the result every later stage uses.
#   * MISFIT AND SWARM ARE PER ENSEMBLE, because each ensemble has its own
#     kernel and its own chains; a pooled swarm would plot each model against
#     the wrong forward operator.
#   * NO CONDUCTIVITY OVERLAY on the posterior. It was the AEM p50 used as the
#     fixed kernel conductivity; here the kernel comes from AEM draws, and a
#     second posterior on the same axes would crowd it.
# The CSVs carry a header row, which smrp_invert.jl's did not.

ENV["MPLBACKEND"] = "Agg"   # headless: no windows, works on a Gadi node

using Distributed, Printf
include(joinpath(@__DIR__, "..", "common", "stage_prelude.jl"))
include(joinpath(COMMON, "prerun_data.jl"))
include(joinpath(COMMON, "prerun_lib.jl"))
include(joinpath(COMMON, "udf_smr.jl"))
using .PrerunLib
using HiQGA.transD_GP, JLD2, Statistics, SMRPInversion, DelimitedFiles, PyPlot
include(joinpath(COMMON, "smr_setup.jl"))

const REPORT     = joinpath(OUT, "report")
const BURNIN     = CFG.BURNINFRAC
const EXTENDFRAC = Float64(FCFG["grid"]["extendfrac"])
const NBINS, QP1, QP2 = 50, 0.05, 0.95

"Cell boundaries for pcolormesh, from the depths the posterior is indexed by."
meshof(z) = [transD_GP.CommonToAll.zcentertoboundary(z); z[end]]

"depth_from / depth_to as smrp_invert.jl wrote them."
function depth_table(zbnd)
    from = collect(Float64, zbnd)
    to   = [from[2:end]; from[end] + (from[end] - from[end-1]) * EXTENDFRAC]
    from, to
end

"The pooled water-content posterior image, with its 5/50/95 % lines."
function plot_pooled(zbnd, CI, himage, edges, title)
    f, ax = plt.subplots(1, 1, figsize = (5, 5))
    ax.pcolormesh(edges[:], meshof(zbnd), himage, cmap = "inferno")
    ax.plot(CI, zbnd, linewidth = 1, color = "w")
    ax.plot(CI, zbnd, linewidth = 1, color = "k", linestyle = "--")
    ax.set_ylim(zbnd[end], zbnd[1])
    ax.set_xlim(extrema(edges)...)
    ax.set_xlabel("fractional water content")
    ax.set_ylabel("depth (m)")
    ax.set_title(title, fontsize = 9)
    f.tight_layout()
    f, ax
end

"V0 and phase against pulse moment, with sigma-hat and the s-scaled noise."
function plot_sounding(D, fid)
    q, V0, ϕ = D["q"], D["V0"] .* 1e9, rem2pi.(D["phi"], RoundNearest)
    σ̂ = D["sigma_hat"] .* 1e9
    # s scales the VARIANCE (SMRPInversion HierarchicalNoise: sigma^2 = s sigma-hat^2),
    # so the noise the posterior settled on is sqrt(s) sigma-hat, not s sigma-hat.
    s50 = isempty(D["log_s"]) ? NaN : median(exp.(D["log_s"]))
    f, ax = plt.subplots(1, 2, sharex = true, figsize = (9, 3.5))
    isfinite(s50) && ax[1].errorbar(q, V0, sqrt(s50) .* σ̂, linestyle = "none",
        color = "0.65", elinewidth = 3,
        label = @sprintf("σ̂·√s  (s p50 %.0f, √s %.1f)", s50, sqrt(s50)))
    ax[1].errorbar(q, V0, σ̂, linestyle = "none", marker = ".", color = "k",
        elinewidth = 1, capsize = 2, label = "σ̂")
    ax[1].set_ylabel("V0 (nV)"); ax[1].legend(fontsize = 8)
    ax[2].plot(q, ϕ, linestyle = "none", marker = ".", color = "k")
    ax[2].set_ylabel("phase (rad)"); ax[2].set_ylim(-π, π)
    for a in ax
        a.set_xscale("log"); a.set_xlabel("pulse moment (A·s)")
        a.grid(true, which = "both", alpha = 0.5)
    end
    f.suptitle("$fid  ($(D["filename"]), $(D["nq"]) moments)", fontsize = 9)
    f.tight_layout()
    f
end

function report_site(fid)
    jf = joinpath(OUT, "smr_transd_$(fid)$(LABEL).jld2")
    isfile(jf) || error("$fid: no $(basename(jf)) -- run run_smr_transd.jl first")
    D = load(jf)
    zbnd, Θ = D["zbnd"], D["theta"]
    size(Θ, 2) == 0 && error("$fid: the pooled posterior holds no draws")
    od = joinpath(REPORT, fid); mkpath(joinpath(od, "swarm")); mkpath(joinpath(od, "misfit"))
    nens = Int(D["nensembles"])

    # Each ensemble's sounding is rebuilt on ITS kernel (udf_kernel, ens = e),
    # exactly as run_smr_transd.jl built it, so the swarm's forward responses
    # are the ones that ensemble was sampled against.
    r = udf_meta(fid)
    path, fname = udf_matfile(fid)
    fname == D["filename"] ||
        @warn "$fid: inversion used $(D["filename"]) but $fname is now preferred; plotting against $fname"
    job = merge(r, (filename = fname, matpath = path))
    b0  = build_udf_sounding(job)
    b0.ok || error("$fid: $(b0.why)")
    S1 = hier(udf_kernel(b0, job, fid, OUT; ens = 1, occam = OCCAMK))

    # ---- pooled posterior + tables ------------------------------------------
    M = [Θ[:, j] for j in axes(Θ, 2)]
    opt0, _ = make_opts("$(D["sounding_id"])$(LABEL)_e1_", zbnd)
    lo, hi = extrema(D["fbounds"])
    # gethimage consults the operator only for an affine stretch, which an SMR
    # sounding never has; ensemble 1's sounding, built above, serves.
    himage, edges, CI, = transD_GP.gethimage(S1, M, opt0; nbins = NBINS, qp1 = QP1,
        qp2 = QP2, rhomin = lo, rhomax = hi, islscale = false, pdfnormalize = false)
    title = @sprintf("%s  pooled: %d ensembles, %d draws", fid, nens, size(Θ, 2))
    f, ax = plot_pooled(zbnd, CI, himage, edges, title)
    savefig(joinpath(od, "$(fid)_post.png"), dpi = 150)
    ax.set_ylim(100, zbnd[1])
    savefig(joinpath(od, "$(fid)_post_100m.png"), dpi = 150)
    close(f)

    from, to = depth_table(zbnd)
    open(joinpath(od, "$(fid)_CI_results.csv"), "w") do io
        println(io, "depth_from_m,depth_to_m,p05,p50,p95")
        writedlm(io, hcat(from, to, CI), ',')
    end
    open(joinpath(od, "$(fid)_himage_results.csv"), "w") do io
        # one row per depth cell; columns are the bin centres
        cen = 0.5 .* (edges[1:end-1] .+ edges[2:end])
        println(io, "depth_from_m,depth_to_m," * join((@sprintf("%.4f", c) for c in cen), ","))
        writedlm(io, hcat(from, to, himage), ',')
    end

    f = plot_sounding(D, fid)
    savefig(joinpath(od, "$(fid)_sounding.png"), dpi = 150); close(f)

    # ---- per ensemble: misfit traces and swarm --------------------------------
    rundir = joinpath(RUNDIR, fid * LABEL)
    for e in 1:nens
        prefix = "$(D["sounding_id"])$(LABEL)_e$(e)_"
        be = udf_kernel(b0, job, fid, OUT; ens = e, occam = OCCAMK)
        S  = hier(be)
        cd(rundir) do
            opt, optn = make_opts(prefix, zbnd)
            set_run_filenames!(opt, optn)
            transD_GP.getchi2forall(opt)
            gcf().suptitle("$fid ensemble $e", fontsize = 9)
            savefig(joinpath(od, "misfit", "$(fid)_e$(e)_misfit_chi2.png"), dpi = 150)
            close("all")
            SMRPInversion.plot_model_field(S, opt, optn; decfactor = 10, lcolor = "k",
                                           modelalpha = 0.08, burninfrac = BURNIN)
            gcf().suptitle("$fid ensemble $e  ($(be.cond_source))", fontsize = 9)
            savefig(joinpath(od, "swarm", "$(fid)_e$(e)_swarm.png"), dpi = 150)
            close("all")
        end
    end
    @printf("  %s: pooled %d draws over %d ensembles -> %s\n", fid, size(Θ, 2), nens, od)
end

function main()
    sites = !isempty(ONLY) ? [ONLY] :
            haskey(ARG, "sites") ? String.(split(ARG["sites"], ",")) : RUN_SITES
    failed = String[]
    for s in sites
        try
            report_site(s)
        catch e
            @warn "report failed for $s" exception = (e, catch_backtrace())
            push!(failed, s)
        end
    end
    isempty(failed) || (@printf("FAILED: %s\n", join(failed, ", ")); exit(1))
    @printf("all %d sites reported\n", length(sites))
end

main()
