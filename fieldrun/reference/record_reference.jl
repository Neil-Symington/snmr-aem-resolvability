# record_reference.jl -- emit the [deterministic] section of a reference file
# from a finished run, so the stored numbers are never hand-typed.
#
#   julia --project=. fieldrun/reference/record_reference.jl \
#       --config fieldrun/config/udf_example.toml >> fieldrun/reference/udf_example.toml
#
# It prints TOML to stdout and nothing else, and it REFUSES to emit a site whose
# Occam solve is missing or non-finite -- a reference recorded from a broken run
# would make every later run agree with the breakage.
#
# Only the Occam numbers go in. They are deterministic: no Monte Carlo, so they
# reproduce to the printed digits anywhere. Posterior quantiles are deliberately
# not recorded; see README.md.

const COMMON = normpath(joinpath(@__DIR__, "..", "common"))
include(joinpath(COMMON, "fieldconfig.jl"))
using .FieldConfig
using JLD2, Printf, Statistics

const CFGPATH = let i = findfirst(==("--config"), ARGS)
    i === nothing && error("record_reference.jl needs --config <control file>")
    ARGS[i+1]
end
const FCFG = FieldConfig.load_config(CFGPATH)
const OUT  = FieldConfig.cfgpath(FCFG, "outputs")
const SITES = String.(FCFG["run"]["sites"])

# The windows run_aem_occam.jl itself reports.
const WINDOWS = [(0.0, 20.0), (20.0, 60.0), (60.0, 120.0)]

geomean_window(m, z, za, zb) = let k = findall(zz -> za <= zz <= zb, z)
    isempty(k) && error("window $za-$zb m contains no cells")
    10^mean(m[k])
end

println("[deterministic]")
for s in SITES
    f = joinpath(OUT, "aem_occam_$(s).jld2")
    isfile(f) || error("no Occam output for $s at $f -- run the batch first")
    A = load(f)
    all(isfinite, A["m_log10sigma"]) ||
        error("$s Occam model has non-finite values; refusing to record it")
    z = Float64.(A["zall"]); m = Float64.(A["m_log10sigma"])
    gm = [geomean_window(m, z, w...) for w in WINDOWS]
    println()
    println("[deterministic.$s]")
    @printf("chi2n = %.12g\n", Float64(A["chi2n"]))
    @printf("ndata = %d\n", Int(A["ndata"]))
    @printf("sigma_windows_m = [[%.0f, %.0f], [%.0f, %.0f], [%.0f, %.0f]]\n",
            WINDOWS[1]..., WINDOWS[2]..., WINDOWS[3]...)
    @printf("sigma_geomean = [%.12g, %.12g, %.12g]\n", gm...)
    # 1e-4 relative is far looser than the run-to-run spread of a deterministic
    # solve and far tighter than any real change in the forward model or the
    # noise specification would be.
    println("tol = 1e-4")
end
