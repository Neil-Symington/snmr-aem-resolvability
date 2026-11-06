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
    @test length(c["run"]["mode_overrides"]) == 3      # smoke overlay recorded
    @test c["mcmc"]["nsamples"] == 2001                # and applied
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
    @test refuses(c -> c["mcmc"]["smoke"]["nsamples"] = -5)
    @test refuses(c -> c["mcmc"]["smoke"]["nchains"] = 0)
    @test refuses(c -> c["mcmc"]["smoke"]["kernel_draws"] = -3)
    @test refuses(c -> c["mcmc"]["smoke"]["nsamples"] = "many")
    @test refuses(c -> c["mcmc"]["smoke"]["nsampels"] = 10)   # overrides nothing
    # and the table must be well formed even when the mode does not use it
    @test refuses(c -> (c["run"]["mode"] = "full"; c["mcmc"]["smoke"]["nsamples"] = -5))
    @test refuses(c -> c["run"]["mode"] = "quick")
end

@testset "valid configs must be ACCEPTED" begin
    # The direction that matters more. Each of these was wrongly refused by the
    # first version of the type and range checks, and a review caught them.
    @test !refuses(c -> c["run"]["minfree_gb"] = 6.5)   # GB, fractional is natural
    @test !refuses(c -> c["grid"]["zstart"] = 0.0)      # build_grid calls setupz(0.0, ...)
    @test !refuses(c -> c["mcmc"]["burninfrac"] = 0.0)  # zero burn-in is meaningful
    @test !refuses(c -> (c["mcmc"]["nmin"] = 40; c["mcmc"]["nmax"] = 40))  # fixed k
    @test !refuses(c -> c["mcmc"]["nsamples"] = 2e5)    # TOML float, integral
    @test !refuses(c -> nothing)                        # the shipped file itself
end
