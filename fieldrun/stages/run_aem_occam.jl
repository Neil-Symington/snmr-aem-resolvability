# Ported from udf_realdata/ into the repository. The include paths and the
# project are the only things changed from the working version: COMMON is
# fieldrun/common/ and the project is the repository root, not ../ek_prerun.
# The control file first: support modules read it at load time.
include(joinpath(@__DIR__, "..", "common", "stage_prelude.jl"))
# run_aem_occam.jl — smooth AEM inversion of the nearest fiducial at each of the
# four chosen UDF SMR sites, with a per-sounding noise calibration.
#
#   julia --project=. run_aem_occam.jl
#
# Writes outputs/aem_occam_<site>.jld2 and prints a summary table. The resulting
# sigma(z) is used twice downstream: as the conductive-earth model for the SNMR
# forward (the SMR kernel is conductivity-sensitive), and as the starting point
# for the forced-conductivity profile that the gate reads.

using JLD2, Printf, Statistics

include(joinpath(COMMON, "aem_lib.jl"))

# --site NAME restricts the run to one site, so the five can run as five
# independent processes; they share nothing but the manifest.

function main()
    mkpath(OUT)
    grid = build_grid()
    @printf("AEM depth grid: %d cells, 0 to %.1f m (last cell %.1f m thick)\n",
            length(grid.zall), grid.zbnd[end], grid.zbnd[end] - grid.zbnd[end-1])
    @printf("system: rTx %.2f m, rRx %.2f m, dzRx %.1f m, %d LM + %d HM gates\n",
            rTx, rRx, dzRx, length(LM_times), length(HM_times))
    @printf("        LM %.2f-%.1f us, HM %.1f-%.1f us\n\n",
            LM_times[1]*1e6, LM_times[end]*1e6, HM_times[1]*1e6, HM_times[end]*1e6)

    rows = []
    for s in manifest()
        (!isempty(ONLY) && s.field_id != ONLY) && continue
        g = load_gates(s.field_id)
        @printf("%s  line %d fid %.1f  sep %.1f m  height %.1f m  (bore %s)\n",
                s.field_id, s.line, s.fiducial, s.sep_m, s.height_m, s.bore)
        t0 = time()
        cal = calibrate_noise(g, s.height_m; grid = grid)
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
        nlm = length(g.d_lm)
        rn = (f.aem.W * f.aem.res)
        @printf("    normalised residual rms: LM %.2f  HM %.2f\n\n",
                sqrt(mean(rn[1:nlm].^2)), sqrt(mean(rn[nlm+1:end].^2)))

        jldsave(joinpath(OUT, "aem_occam_$(s.field_id).jld2");
            field_id = s.field_id, zall = grid.zall, zbnd = grid.zbnd,
            m_log10sigma = f.m, sigma = σz,
            chi2 = f.chi2, ndata = f.ndata, chi2n = f.chi2 / f.ndata,
            noise_scale = cal.scale,
            d_lm = g.d_lm, sigma_lm = g.σ_lm, g_lm = g.g_lm,
            d_hm = g.d_hm, sigma_hm = g.σ_hm, g_hm = g.g_hm,
            pred_lm = (f.aem.Flow.dBzdt[f.aem.selectlow]) .* (MU0 / PV),
            pred_hm = (f.aem.Fhigh.dBzdt[f.aem.selecthigh]) .* (MU0 / PV),
            height_m = s.height_m,
            line = s.line, fiducial = s.fiducial, sep_m = s.sep_m,
            easting = s.easting, northing = s.northing, dtm_m = s.dtm_m,
            plm = s.plm, bore = s.bore, bore_sep_m = s.bore_sep_m,
            LM_times = LM_times, HM_times = HM_times)
        push!(rows, (s.field_id, cal.scale, f.chi2 / f.ndata, f.ndata,
                     unit_geomean(f.m, grid.zall, 0, 20),
                     unit_geomean(f.m, grid.zall, 20, 60),
                     unit_geomean(f.m, grid.zall, 60, 120)))
    end

    println("="^96)
    @printf("%-14s %8s %9s %6s %10s %10s %10s\n", "site", "sig x", "chi2/N",
            "ndata", "0-20 m", "20-60 m", "60-120 m")
    for r in rows
        @printf("%-14s %8.3f %9.4f %6d %10.4f %10.4f %10.4f\n", r...)
    end
    println("="^96)
    println("sigma(z) in S/m, geometric mean over the window")
end

main()
