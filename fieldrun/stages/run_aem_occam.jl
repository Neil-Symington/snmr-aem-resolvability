# Ported from udf_realdata/ into the repository. The include paths and the
# project are the only things changed from the working version: COMMON is
# fieldrun/common/ and the project is the repository root, not ../ek_prerun.
# The control file first: support modules read it at load time.
include(joinpath(@__DIR__, "..", "common", "stage_prelude.jl"))
# run_aem_occam.jl — smooth AEM inversion of the nearest fiducial at each of the
# SMR sites, with a per-sounding noise calibration. System-independent: [aem]
# system picks SkyTEM (UDF) or XCite (Birrindudu) through common/aem_lib.jl.
#
#   julia --project=. run_aem_occam.jl
#
# Writes outputs/aem_occam_<site>.jld2 and prints a summary table. The resulting
# sigma(z) is used twice downstream: as the conductive-earth model for the SNMR
# forward (the SMR kernel is conductivity-sensitive), and as the starting point
# for the forced-conductivity profile that the gate reads.

using JLD2, Printf, Statistics

include(joinpath(COMMON, "aem_lib.jl"))

# [aem] noise_calibration: "inflate" (default, the scalar floored at 1 that
# calibrate_noise has always applied) or "two_term" (XCite: the nominal_unc
# convention rescaled either way, the measured halt floor kept, solved for
# chi2_min = ndata - p_eff; aem_lib.jl calibrate_noise_twoterm). Read here,
# not merely validated -- the [detect] section was validated and never read.
const NOISE_CAL = String(get(FCFG["aem"], "noise_calibration", "inflate"))

# --site NAME restricts the run to one site, so the five can run as five
# independent processes; they share nothing but the manifest.

function main()
    mkpath(OUT)
    grid = build_grid()
    @printf("AEM depth grid: %d cells, 0 to %.1f m (last cell %.1f m thick)\n",
            length(grid.zall), grid.zbnd[end], grid.zbnd[end] - grid.zbnd[end-1])
    describe(); println()     # the system module says what it parsed

    rows = []
    for s in manifest()
        (!isempty(ONLY) && s.field_id != ONLY) && continue
        g = load_gates(s.field_id)
        @printf("%s  line %d fid %.1f  sep %.1f m  height %.1f m  (bore %s)\n",
                s.field_id, s.line, s.fiducial, s.sep_m, s.height_m, s.bore)
        t0 = time()
        local cal, calrec
        if NOISE_CAL == "two_term"
            println("    two-term noise calibration (nominal_unc rescaled, halt floor kept):")
            tc = calibrate_noise_twoterm(g, s.height_m; grid = grid)
            h = tc.history[end]
            cal = (scale = NaN, fit = tc.fit)
            calrec = (noise_calibration = "two_term", sigma_cal = tc.σ, noise_conv_factor = tc.f,
                      p_eff = tc.p_eff, lambda2 = h.l2, roughness = h.roughness, n_on_bound = h.n_bound,
                      cal_history = [collect(values(x)) for x in tc.history],
                      cal_history_cols = collect(String.(keys(tc.history[1]))))
            @printf("    convention x %.4f  chi2 %.3f = ndata %d - p_eff %.2f (ratio %.4f)  lambda2 %.3g  roughness %.4g  %d cells on a bound  (%.1f s)\n",
                    tc.f, tc.fit.chi2, tc.fit.ndata, tc.p_eff, h.ratio, h.l2, h.roughness, h.n_bound, time() - t0)
        else
            cal = calibrate_noise(g, s.height_m; grid = grid)
            calrec = (noise_calibration = "inflate",)
        end
        f = cal.fit
        σz = 10 .^ f.m
        @printf("    noise scale %.3f  chi2/N %.4f  ndata %d  (%.1f s)\n",
                cal.scale, f.chi2 / f.ndata, f.ndata, time() - t0)
        @printf("    sigma(z): 0-20 m %.4f   20-60 m %.4f   60-120 m %.4f S/m\n",
                unit_geomean(f.m, grid.zall, 0, 20),
                unit_geomean(f.m, grid.zall, 20, 60),
                unit_geomean(f.m, grid.zall, 60, 120))
        # Residual structure, per moment. aem.res after occam() holds whatever
        # the LAST candidate of the line search left there, not the accepted
        # model, so the forward is re-run at the accepted model first.
        chi2(f.aem, f.m)
        rn = (f.aem.W * f.aem.res)
        println("    normalised residual rms: ",
                join([@sprintf("%s %.2f", b, sqrt(mean(rn[ix].^2))) for (b, ix) in data_blocks(g)], "  "), "\n")
        pred = predicted(f.aem)
        # SkyTEM's per-moment keys, kept so files written before the system
        # switch and after it read the same way. XCite has one moment.
        legacy = AEM_SYSTEM == "skytem" ?
            (d_lm = g.d_lm, sigma_lm = g.σ_lm, g_lm = g.g_lm,
             d_hm = g.d_hm, sigma_hm = g.σ_hm, g_hm = g.g_hm,
             pred_lm = pred[1:length(g.d_lm)], pred_hm = pred[length(g.d_lm)+1:end]) : (;)

        jldsave(joinpath(OUT, "aem_occam_$(s.field_id).jld2");
            field_id = s.field_id, zall = grid.zall, zbnd = grid.zbnd,
            m_log10sigma = f.m, sigma = σz,
            chi2 = f.chi2, ndata = f.ndata, chi2n = f.chi2 / f.ndata,
            noise_scale = cal.scale, aem_system = AEM_SYSTEM, calrec...,
            d = g.d, sigma_d = g.σ, gate = g.gate, moment = g.moment, pred = pred,
            legacy...,
            height_m = s.height_m,
            line = s.line, fiducial = s.fiducial, sep_m = s.sep_m,
            easting = s.easting, northing = s.northing, dtm_m = s.dtm_m,
            plm = s.plm, bore = s.bore, bore_sep_m = s.bore_sep_m,
            gate_times()...)
        push!(rows, (s.field_id, NOISE_CAL == "two_term" ? calrec.noise_conv_factor : cal.scale, f.chi2 / f.ndata, f.ndata,
                     unit_geomean(f.m, grid.zall, 0, 20),
                     unit_geomean(f.m, grid.zall, 20, 60),
                     unit_geomean(f.m, grid.zall, 60, 120)))
    end

    println("="^96)
    @printf("%-14s %8s %9s %6s %10s %10s %10s\n", "site", NOISE_CAL == "two_term" ? "conv x" : "sig x", "chi2/N",
            "ndata", "0-20 m", "20-60 m", "60-120 m")
    for r in rows
        @printf("%-14s %8.3f %9.4f %6d %10.4f %10.4f %10.4f\n", r...)
    end
    println("="^96)
    println("sigma(z) in S/m, geometric mean over the window")
end

main()
