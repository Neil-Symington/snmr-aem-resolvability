# test_fieldconfig.jl -- the control-file reader's failure paths.
#
# fieldrun/README.md claimed these were "each negative-tested". They were, by
# reviewers, on their own machines -- nothing was committed and no CI step
# touched this file, which a review pointed out. This makes the claim true.
#
# Each case below is a defect that actually got through: three rounds of review
# found 36 escapes in two directions, and the WRONGLY-REFUSED ones at the end
# matter more than the accepted ones, because they break working configs.
using Test, TOML

const FR   = normpath(joinpath(@__DIR__, ".."))
const BASE = TOML.parsefile(joinpath(FR, "config", "udf_example.toml"))

refuses(mut) = begin
    c = deepcopy(BASE); mut(c)
    p = tempname() * ".toml"
    open(p, "w") do f; TOML.print(f, c); end
    r = try (FieldConfig.load_config(p); false) catch; true end
    rm(p, force = true)
    r
end

@testset "the shipped control file loads" begin
    c = FieldConfig.load_config(joinpath(FR, "config", "udf_example.toml"))
    @test c["run"]["name"] == "udf_example"
    @test length(c["run"]["mode_overrides"]) == 4      # smoke overlay recorded
    @test c["mcmc"]["nsamples_aem"]  == 2001           # and applied
    @test c["mcmc"]["nsamples_snmr"] == 2001           # both of them
    @test !haskey(c["mcmc"], "smoke")                  # and consumed
end

@testset "unknown sections and keys" begin
    @test refuses(c -> c["wibble"] = Dict("a" => 1))
    @test refuses(c -> c["grid"]["nzz"] = 60)
    @test refuses(c -> c["mcmc"]["nsample"] = 1000)     # near miss
    @test refuses(c -> c["grid"]["nchains"] = 4)        # right key, wrong section
    @test refuses(c -> c["run"]["mode_overrides"] = ["x"])   # loader-written
    @test refuses(c -> c["run"]["config_file"] = "x")        # loader-written
end

@testset "missing sections the stages need" begin
    for sec in ("aem", "snmr", "grid", "mcmc", "gate", "detect")
        @test refuses(c -> delete!(c, sec))
    end
end

@testset "types" begin
    @test refuses(c -> c["grid"]["nz"] = true)          # true isa Integer in Julia
    @test refuses(c -> c["grid"]["nz"] = 2.5)
    @test refuses(c -> c["grid"]["nz"] = "60")          # was a raw MethodError
    @test refuses(c -> c["grid"]["dz"] = Inf)
    @test refuses(c -> c["mcmc"]["seed_base"] = "abc")
    @test refuses(c -> c["run"]["name"] = 5)            # was a raw MethodError
    @test refuses(c -> c["gate"]["thresholds_mgL"] = ["a"])   # was a raw MethodError
    @test refuses(c -> c["aem"]["lm_gates"] = ["a", "b"])     # was a raw MethodError
    @test refuses(c -> c["gate"]["thresholds_mgL"] = 1000)    # scalar, not a list
    @test refuses(c -> c["snmr"]["mat_preference"] = [1, 2])
end

@testset "values" begin
    @test refuses(c -> c["run"]["npar"] = -1)
    @test refuses(c -> c["detect"]["theta"] = 1.8)
    @test refuses(c -> c["mcmc"]["tmax"] = 1.0)                    # no ladder
    @test refuses(c -> c["mcmc"]["fbounds"] = [0.45, 0.001])       # reversed
    @test refuses(c -> c["aem"]["lm_gates"] = [26, 9])             # reversed
    @test refuses(c -> c["aem"]["lm_gates"] = [9.5, 26])           # fractional gate
    @test refuses(c -> (c["aem"]["lm_ramp_hold"] = 15e-3; c["aem"]["hm_ramp_hold"] = 1e-3))
    @test refuses(c -> c["gate"]["thresholds_mgL"] = [3000, 1000]) # unsorted
    @test refuses(c -> c["gate"]["thresholds_mgL"] = [1000, 1000]) # duplicate
    @test refuses(c -> c["detect"]["metric"] = "pointwise")
    @test refuses(c -> c["snmr"]["loop_shape"] = "hexagon")
    @test refuses(c -> c["grid"]["extendfrac"] = 0.9)
    @test refuses(c -> c["mcmc"]["lambda_gp"] = [-3.0])
    @test refuses(c -> c["mcmc"]["delta"] = -1.0)
    @test refuses(c -> c["noise"]["s_sdev"] = [-0.3])
    @test refuses(c -> c["noise"]["s_prior_tau"] = 0.0)
    @test refuses(c -> c["occam"]["beta2"] = -1.0)
    @test refuses(c -> c["occam"]["regtype"] = "R9")
    @test refuses(c -> c["occam"]["uncertainty"] = "lapalce")
    @test refuses(c -> c["noise"]["sigma_hat_branch"] = "nonsense")
    @test refuses(c -> c["detect"]["min_thick_m"] = -5.0)
    @test refuses(c -> c["detect"]["occam_chi2_max"] = -1.0)
    @test refuses(c -> c["detect"]["doi_widthfrac_alt"] = [0.9, 0.1])
    @test refuses(c -> c["convergence"]["rhat"] = 0.5)
    @test refuses(c -> c["convergence"]["depths_m"] = [-5.0, 15.0])
    @test refuses(c -> c["aem"]["rx_radius_m"] = -13.29)
    @test refuses(c -> c["aem"]["lowpass_hz"] = [-1.0, 3e5])
    @test refuses(c -> c["run"]["name"] = "")
    @test refuses(c -> c["run"]["name"] = "../../x")    # becomes a directory
end

@testset "sites and stage order" begin
    @test refuses(c -> c["run"]["sites"] = ["SMR_UD21", "SMR_UD99"])   # not in sites.csv
    @test refuses(c -> c["run"]["sites"] = String[])
    @test refuses(c -> c["run"]["sites"] = ["SMR_UD21", "SMR_UD21"])   # duplicate
    @test refuses(c -> c["run"]["stages"] = String[])
    @test refuses(c -> c["run"]["stages"] = ["aem_occam", "aem_occam"]) # duplicate
    @test refuses(c -> c["run"]["stages"] = ["profile", "detect"])      # out of order
    @test refuses(c -> c["run"]["stages"] = ["aem_occam", "profile"])   # profile needs detect
    @test refuses(c -> c["run"]["stages"] = ["aem_occam", "wibble"])
end

@testset "the smoke overlay" begin
    # check_values used to run BEFORE the overlay, so these went straight through
    @test refuses(c -> c["mcmc"]["smoke"]["nsamples_aem"] = -5)
    @test refuses(c -> c["mcmc"]["smoke"]["nchains"] = 0)
    @test refuses(c -> c["mcmc"]["smoke"]["kernel_draws"] = -3)
    @test refuses(c -> c["mcmc"]["smoke"]["nsamples_snmr"] = "many")
    @test refuses(c -> c["mcmc"]["smoke"]["nsampels"] = 10)   # overrides nothing
    # the two counts are separate keys, so a smoke table naming the old
    # shared one overrides nothing and must be refused
    @test refuses(c -> c["mcmc"]["smoke"]["nsamples"] = 2001)
    # and the table must be well formed even when the mode does not use it
    @test refuses(c -> (c["run"]["mode"] = "full"; c["mcmc"]["smoke"]["nsamples_aem"] = -5))
    @test refuses(c -> c["run"]["mode"] = "quick")
end

@testset "valid configs must be ACCEPTED" begin
    # The direction that matters more. Each of these was wrongly refused by the
    # first version of the type and range checks, and a review caught them.
    @test !refuses(c -> c["run"]["minfree_gb"] = 6.5)   # GB, fractional is natural
    @test !refuses(c -> c["grid"]["zstart"] = 0.0)      # build_grid calls setupz(0.0, ...)
    @test !refuses(c -> c["mcmc"]["burninfrac"] = 0.0)  # zero burn-in is meaningful
    @test !refuses(c -> (c["mcmc"]["nmin"] = 40; c["mcmc"]["nmax"] = 40))  # fixed k
    @test !refuses(c -> c["mcmc"]["nsamples_aem"] = 2e5)  # TOML float, integral
    @test !refuses(c -> nothing)                        # the shipped file itself
end

# ─────────────────────────────────────────────────────────────────────────────
# The control file is not the only declaration of these parameters.
#
# `src/prerun_config.jl` (module PrerunConfig) hardcodes the survey pre-run
# spec, and `fieldrun/common/prerun_data.jl` aliases it to `CFG`. So a stage
# writing `CFG.SAVE_FREQ` reads THAT FILE, not the control file, and nothing at
# the call site says which source won. Thirty-two values are declared in both
# places. Thirty agree, which is exactly why nobody noticed; where they
# disagree the code wins silently, and that is how `s_prior_tau` shipped at 0.5
# in the control file while the stage ran 2.0 (fixed, in the other direction).
#
# These tests pin the two copies together. Edit one and this fails, naming the
# parameter. They are deliberately TEXTUAL -- prerun_config.jl is parsed as a
# string rather than loaded -- so the test asserts what a reader of the file
# sees and costs no HiQGA load.
const PCSRC = read(joinpath(FR, "..", "src", "prerun_config.jl"), String)

"The literal right-hand side of `const NAME = ...`, comment stripped."
pcraw(name) = begin
    m = match(Regex("^const\\s+" * name * "\\s*=\\s*([^#\\n]+)", "m"), PCSRC)
    m === nothing && error("prerun_config.jl declares no const $name")
    strip(m.captures[1])
end

"Every number in that right-hand side, so `[0.5, 0.7]` and `0.5` compare alike."
pcnums(name) = [parse(Float64, replace(m.match, "_" => "")) for m in
                eachmatch(r"-?\d+(?:_\d+)*(?:\.\d+)?(?:[eE][-+]?\d+)?", pcraw(name))]

tonums(v) = v isa AbstractArray ? Float64.(collect(Iterators.flatten(v))) : [Float64(v)]

@testset "control file and prerun_config.jl declare the same values" begin
    # BASE, not the loaded config: load_config applies the smoke overlay, which
    # would make nsamples_snmr 2001 against the module's 150001 and report a
    # disagreement that is not one.
    for (sec, key, name) in (
            ("grid", "nz",             "NZ"),
            ("grid", "zstart",         "ZSTART"),
            ("grid", "dz",             "DZ"),
            ("grid", "extendfrac",     "EXTENDFRAC"),
            ("mcmc", "nchains",        "NCHAINS"),
            ("mcmc", "tmax",           "TMAX"),
            ("mcmc", "nsamples_snmr",  "NSAMPLES"),
            ("mcmc", "burninfrac",     "BURNINFRAC"),
            ("mcmc", "save_freq",      "SAVE_FREQ"),
            ("mcmc", "seed_base",      "SEED_BASE"),
            ("mcmc", "lambda_gp",      "LAMBDA_GP"),
            ("mcmc", "delta",          "DELTA"),
            ("mcmc", "fbounds",        "FBOUNDS"),
            ("noise", "s_prior_center", "S_PRIOR_CENTER"),
            ("noise", "s_sdev",         "S_SDEV"),
            ("detect", "theta",             "DETECT_THETA"),
            ("detect", "min_thick_m",       "DETECT_MINTHICK"),
            ("detect", "occam_chi2_max",    "DETECT_OCCAM_CHI2"),
            ("detect", "doi_widthfrac",     "DOI_WIDTHFRAC"),
            ("detect", "doi_widthfrac_alt", "DOI_WIDTHFRAC_ALT"),
            ("occam", "w0",    "OCCAM_W0"),
            ("occam", "beta2", "OCCAM_BETA2"),
            ("convergence", "rhat",         "GATE_RHAT"),
            ("convergence", "min_ess",      "GATE_MIN_ESS"),
            ("convergence", "splithalf_dw", "GATE_SPLITHALF_DW"),
            ("convergence", "trend_sigma",  "TREND_SIGMA"),
            ("convergence", "depths_m",     "GATE_DEPTHS"),
        )
        @test tonums(BASE[sec][key]) ≈ pcnums(name)
    end
    # declared as one line, `const NMIN, NMAX = 2, 40`
    @test [Float64(BASE["mcmc"]["nmin"]), Float64(BASE["mcmc"]["nmax"])] ≈
          pcnums("NMIN,\\s*NMAX")
    # and as two separate consts
    @test tonums(BASE["occam"]["lam2_log10"]) ≈
          vcat(pcnums("OCCAM_LAM2MIN"), pcnums("OCCAM_LAM2MAX"))

    # ── KNOWN STALE COPIES ──────────────────────────────────────────────────
    # These two disagree, and the CONTROL FILE is the one that runs:
    # run_smr_transd.jl:127-133 reads both from [noise], and its header records
    # why tau = 2.0 and an upper bound of 1e4 are the right values. The module's
    # copies are left because the East Kimberley pre-run was produced with them.
    #
    # @test_broken, so that deleting or correcting a stale copy shows up as an
    # "Unexpected Pass" and asks for its line here to be removed rather than
    # quietly agreeing.
    @test_broken tonums(BASE["noise"]["s_prior_tau"]) ≈ pcnums("S_PRIOR_TAU")
    # log10 in the control file, natural log in the module, so compare the
    # bounds themselves: 1e-2..1e4 against 1e-2..1e2.
    @test_broken (10.0 .^ tonums(BASE["noise"]["s_bounds_log10"])) ≈ pcnums("S_BOUNDS")
end
