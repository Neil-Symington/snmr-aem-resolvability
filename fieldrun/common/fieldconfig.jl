# fieldconfig.jl — the control file, read once, validated, and shared.
#
# Every ported stage reads its parameters from here instead of from constants
# scattered across four files. Two things this must do that a bare TOML parse
# does not:
#
#   1. RESOLVE PATHS against fieldrun/, so a config is portable and a run
#      started from any directory reads the same inputs.
#   2. APPLY THE MODE. `mode = "smoke"` overlays [mcmc.smoke] onto [mcmc], and
#      the overlay is RECORDED in the returned config, so a smoke result can
#      never be mistaken for a production one downstream.
#
# LOAD ORDER MATTERS. udf_system.jl parses the SkyTEM .des at load time — on
# purpose, so a transcription error cannot survive — which means it needs the
# path before it is included. So this file is included FIRST and sets CFG, and
# udf_system.jl reads CFG[]. Including them the other way round fails with a
# message that says so rather than with an UndefVarError.

module FieldConfig

using TOML, Printf, DelimitedFiles

export load_config, CFG, cfgpath, is_smoke

const HERE = normpath(joinpath(@__DIR__, ".."))

# Set by load_config so the load-time consumers (udf_system.jl) can see it.
const CFG = Ref{Union{Nothing,Dict{String,Any}}}(nothing)

const REQUIRED = [
    ["run", "name"], ["run", "sites"], ["run", "stages"], ["run", "mode"],
    ["paths", "sites"], ["paths", "snmr_meta"], ["paths", "snmr_dir"],
    ["paths", "aem_sounds"], ["paths", "aem_system"], ["paths", "outputs"],
    ["grid", "nz"], ["grid", "zstart"], ["grid", "dz"], ["grid", "extendfrac"],
    ["mcmc", "nchains"], ["mcmc", "nsamples"], ["mcmc", "seed_base"],
    ["noise", "sigma_hat_branch"], ["kernel", "n_draws"],
    ["detect", "theta"], ["gate", "thresholds_mgL"],
]

# Relative [paths] entries are resolved against fieldrun/. These are stems or
# directories rather than single files, so existence is checked per use.
const PATH_KEYS = ("sites", "snmr_meta", "snmr_dir", "aem_sounds", "aem_nav",
                   "aem_system", "outputs", "logs")

getin(d, ks) = foldl((a, k) -> a === nothing ? nothing : get(a, k, nothing), ks; init = d)

# Every section and key the schema knows. A review found that an unknown
# section, an unknown key and an out-of-range VALUE were all accepted silently,
# which for a file whose only job is to state parameters is the worst available
# failure -- the same defect the truth generator had.
const KNOWN = Dict(
    "run"         => ["name", "sites", "stages", "mode", "npar", "minfree_gb",
                      "mode_overrides", "config_file"],
    "paths"       => ["sites", "snmr_meta", "snmr_dir", "aem_sounds", "aem_nav",
                      "aem_system", "outputs", "logs"],
    "snmr"        => ["mat_preference", "min_pulse_moments", "loop_shape"],
    "aem"         => ["nominal_unc", "lm_gates", "hm_gates", "lm_ramp_hold",
                      "hm_ramp_hold", "lowpass_hz", "tx_area_m2", "rx_radius_m",
                      "rx_dz_m", "max_sep_m"],
    "grid"        => ["nz", "zstart", "dz", "extendfrac"],
    "mcmc"        => ["nchains", "tmax", "nsamples", "burninfrac", "save_freq",
                      "nensembles", "seed_base", "nmin", "nmax", "lambda_gp",
                      "delta", "fbounds", "smoke"],
    "noise"       => ["sigma_hat_branch", "s_prior_center", "s_prior_tau",
                      "s_sdev", "s_bounds_log10"],
    "kernel"      => ["n_draws"],
    "occam"       => ["regtype", "w0", "lam2_log10", "beta2", "uncertainty"],
    "detect"      => ["theta", "min_thick_m", "occam_chi2_max", "doi_widthfrac",
                      "doi_widthfrac_alt", "metric"],
    "gate"        => ["thresholds_mgL"],
    "convergence" => ["rhat", "min_ess", "splithalf_dw", "splithalf_max",
                      "trend_sigma", "depths_m"],
)

nearest(k, known) = begin
    c = [x for x in known if first(x) == first(k) && abs(length(x) - length(k)) <= 2]
    isempty(c) ? "" : "  did you mean $(join(sort(c), " or "))?"
end

"Refuse any section or key the schema does not know."
function check_keys(c::Dict, path)
    bad = String[]
    for sec in keys(c)
        haskey(KNOWN, sec) || (push!(bad, "[$sec] is not a known section$(nearest(sec, collect(keys(KNOWN))))"); continue)
        c[sec] isa AbstractDict || (push!(bad, "[$sec] is not a table"); continue)
        for k in keys(c[sec])
            k in KNOWN[sec] || push!(bad, "[$sec] $k$(nearest(k, KNOWN[sec]))")
        end
    end
    isempty(bad) || error("$path: unknown section(s) or key(s):\n  " * join(bad, "\n  ") *
        "\n\nAn ignored key changes nothing and warns about nothing, so this is refused.")
end

"Refuse values that are the right type and the wrong number."
function check_values(c::Dict, path)
    bad = String[]
    chk(cond, msg) = cond || push!(bad, msg)
    # Sections that are not in REQUIRED are optional by the schema, but
    # check_values indexed c["aem"] and c["snmr"] unguarded, so deleting either
    # gave a raw KeyError instead of a message. A review found both. Sections
    # the code genuinely needs are now named as required, and everything else
    # is probed with haskey.
    for sec in ("aem", "snmr", "grid", "mcmc", "gate", "detect")
        haskey(c, sec) || push!(bad, "[$sec] is missing, and $(sec == "aem" ?
            "stages/udf_system.jl reads it" : "the stages read it")")
    end
    isempty(bad) || error("$path: missing section(s):\n  " * join(bad, "\n  "))
    pos(sec, k) = haskey(c[sec], k) && chk(c[sec][k] > 0, "[$sec] $k must be > 0, got $(c[sec][k])")
    frac(sec, k) = haskey(c[sec], k) && chk(0 < c[sec][k] <= 1, "[$sec] $k must be in (0,1], got $(c[sec][k])")

    # TYPES FIRST. A review got nz = true, nz = 2.5, npar = 1.5, dz = inf and
    # seed_base = "abc" past every numeric check, and nz = "60" produced a raw
    # MethodError rather than a message. Bool is checked before Integer on
    # purpose: in Julia `true isa Integer` is true.
    ints   = (("run","npar"), ("run","minfree_gb"), ("grid","nz"), ("mcmc","nchains"),
              ("mcmc","nsamples"), ("mcmc","save_freq"), ("mcmc","nensembles"),
              ("mcmc","seed_base"), ("mcmc","nmin"), ("mcmc","nmax"),
              ("kernel","n_draws"), ("snmr","min_pulse_moments"),
              ("convergence","min_ess"))
    reals  = (("grid","dz"), ("grid","zstart"), ("grid","extendfrac"), ("mcmc","tmax"),
              ("mcmc","burninfrac"), ("aem","nominal_unc"), ("aem","tx_area_m2"),
              ("aem","rx_radius_m"), ("aem","rx_dz_m"), ("aem","max_sep_m"),
              ("aem","lm_ramp_hold"), ("aem","hm_ramp_hold"), ("detect","theta"),
              ("detect","min_thick_m"), ("detect","occam_chi2_max"),
              ("detect","doi_widthfrac"), ("convergence","rhat"),
              ("convergence","splithalf_dw"), ("convergence","splithalf_max"),
              ("convergence","trend_sigma"), ("noise","s_prior_center"),
              ("noise","s_prior_tau"), ("occam","w0"), ("occam","beta2"))
    for (sec, k) in ints
        haskey(c, sec) && haskey(c[sec], k) || continue
        v = c[sec][k]
        v isa Bool && (chk(false, "[$sec] $k must be an integer, got the boolean $v"); continue)
        v isa Integer || chk(false, "[$sec] $k must be an integer, got $(repr(v))")
    end
    for (sec, k) in reals
        haskey(c, sec) && haskey(c[sec], k) || continue
        v = c[sec][k]
        v isa Bool && (chk(false, "[$sec] $k must be a number, got the boolean $v"); continue)
        (v isa Real && isfinite(v)) ||
            chk(false, "[$sec] $k must be a finite number, got $(repr(v))")
    end
    isempty(bad) || error("$path: bad type(s):\n  " * join(bad, "\n  "))

    for (sec, k) in (("run","npar"), ("run","minfree_gb"), ("grid","nz"), ("grid","dz"),
                     ("grid","zstart"), ("mcmc","nchains"), ("mcmc","nsamples"),
                     ("mcmc","save_freq"), ("mcmc","nensembles"), ("kernel","n_draws"),
                     ("snmr","min_pulse_moments"), ("aem","tx_area_m2"),
                     ("aem","max_sep_m"), ("convergence","min_ess"))
        haskey(c, sec) && pos(sec, k)
    end
    for (sec, k) in (("aem","nominal_unc"), ("detect","theta"), ("detect","doi_widthfrac"),
                     ("mcmc","burninfrac"))
        haskey(c, sec) && frac(sec, k)
    end
    haskey(c["grid"], "extendfrac") && chk(c["grid"]["extendfrac"] >= 1,
        "[grid] extendfrac must be >= 1, got $(c["grid"]["extendfrac"])")
    haskey(c["mcmc"], "tmax") && chk(c["mcmc"]["tmax"] > 1,
        "[mcmc] tmax must be > 1 for a temperature ladder, got $(c["mcmc"]["tmax"])")
    if haskey(c["mcmc"], "nmin") && haskey(c["mcmc"], "nmax")
        chk(c["mcmc"]["nmin"] < c["mcmc"]["nmax"],
            "[mcmc] nmin must be < nmax, got $(c["mcmc"]["nmin"]) and $(c["mcmc"]["nmax"])")
    end
    for (sec, k) in (("mcmc","fbounds"), ("occam","lam2_log10"), ("noise","s_bounds_log10"))
        if haskey(c, sec) && haskey(c[sec], k)
            v = c[sec][k]
            chk(length(v) == 2 && v[1] < v[2], "[$sec] $k must be [lo, hi] with lo < hi, got $v")
        end
    end
    for (sec, k) in (("aem","lm_gates"), ("aem","hm_gates"))
        if haskey(c[sec], k)
            v = c[sec][k]
            chk(length(v) == 2 && v[1] >= 1 && v[1] <= v[2],
                "[$sec] $k must be [first, last] gate numbers with 1 <= first <= last, got $v")
        end
    end
    haskey(c["aem"], "lm_ramp_hold") && haskey(c["aem"], "hm_ramp_hold") &&
        chk(0 < c["aem"]["lm_ramp_hold"] < c["aem"]["hm_ramp_hold"],
            "[aem] ramp holds must satisfy 0 < lm < hm, got $(c["aem"]["lm_ramp_hold"]) and $(c["aem"]["hm_ramp_hold"])")
    if haskey(c["gate"], "thresholds_mgL")
        v = c["gate"]["thresholds_mgL"]
        chk(!isempty(v) && all(>(0), v) && issorted(v),
            "[gate] thresholds_mgL must be positive and ascending, got $v")
    end
    haskey(c, "detect") && haskey(c["detect"], "metric") &&
        chk(c["detect"]["metric"] in ("integral", "threshold"),
            "[detect] metric must be \"integral\" or \"threshold\", got \"$(c["detect"]["metric"])\"")
    haskey(c["snmr"], "loop_shape") &&
        chk(c["snmr"]["loop_shape"] in ("square", "circle"),
            "[snmr] loop_shape must be \"square\" or \"circle\", got \"$(c["snmr"]["loop_shape"])\"")

    # Checks a review found entirely absent.
    haskey(c["run"], "name") && chk(!isempty(strip(c["run"]["name"])),
        "[run] name must not be empty")
    chk(length(unique(c["run"]["sites"])) == length(c["run"]["sites"]),
        "[run] sites contains duplicates: $(c["run"]["sites"])")
    chk(!isempty(c["run"]["stages"]), "[run] stages must not be empty")
    chk(length(unique(c["run"]["stages"])) == length(c["run"]["stages"]),
        "[run] stages contains duplicates: $(c["run"]["stages"])")
    for k in ("outputs", "logs")
        haskey(c["paths"], k) && chk(!isempty(strip(string(c["paths"][k]))),
            "[paths] $k must not be empty")
    end
    if haskey(c["gate"], "thresholds_mgL")
        v = c["gate"]["thresholds_mgL"]
        chk(length(unique(v)) == length(v), "[gate] thresholds_mgL contains duplicates: $v")
    end
    if haskey(c["aem"], "lowpass_hz")
        chk(all(x -> x isa Real && isfinite(x) && x > 0, c["aem"]["lowpass_hz"]),
            "[aem] lowpass_hz must all be positive and finite, got $(c["aem"]["lowpass_hz"])")
    end
    for (sec, k) in (("occam","regtype"), ("occam","uncertainty"),
                     ("noise","sigma_hat_branch"))
        haskey(c, sec) && haskey(c[sec], k) &&
            chk(c[sec][k] isa AbstractString && !isempty(c[sec][k]),
                "[$sec] $k must be a non-empty string, got $(repr(c[sec][k]))")
    end
    haskey(c["snmr"], "mat_preference") &&
        chk(c["snmr"]["mat_preference"] isa AbstractVector &&
            !isempty(c["snmr"]["mat_preference"]),
            "[snmr] mat_preference must be a non-empty list of patterns")
    # Keys this module WRITES, which a config must not supply.
    for k in ("mode_overrides", "config_file")
        haskey(c["run"], k) && chk(false,
            "[run] $k is written by the loader and must not be set in a config")
    end

    isempty(bad) || error("$path: bad value(s):\n  " * join(bad, "\n  "))
end

"""
    load_config(tomlfile) -> Dict

Read, validate and normalise a control file. Errors name the offending key.
"""
function load_config(tomlfile::AbstractString)
    isfile(tomlfile) || error("no control file at $tomlfile")
    c = TOML.parsefile(tomlfile)

    check_keys(c, tomlfile)
    for ks in REQUIRED
        getin(c, ks) === nothing && error("control file $tomlfile is missing [$(ks[1])] $(ks[2])")
    end
    check_values(c, tomlfile)

    # --- paths ------------------------------------------------------------
    for k in PATH_KEYS
        haskey(c["paths"], k) || continue
        p = c["paths"][k]
        c["paths"][k] = isabspath(p) ? normpath(p) : normpath(joinpath(HERE, p))
    end
    for (k, what) in (("sites", "site table"), ("snmr_meta", "SNMR metadata"))
        isfile(c["paths"][k]) || error("[paths] $k: no $what at $(c["paths"][k])")
    end
    for (k, what) in (("snmr_dir", "SNMR directory"), ("aem_sounds", "AEM sounding directory"))
        isdir(c["paths"][k]) || error("[paths] $k: no $what at $(c["paths"][k])")
    end
    des = c["paths"]["aem_system"] * ".des"
    isfile(des) || error("[paths] aem_system: no .des at $des (the stem takes no extension)")

    # --- mode -------------------------------------------------------------
    mode = c["run"]["mode"]
    mode in ("smoke", "full") || error("[run] mode must be \"smoke\" or \"full\", got \"$mode\"")
    applied = String[]
    if mode == "smoke"
        ov = get(c["mcmc"], "smoke", Dict{String,Any}())
        isempty(ov) && error("[run] mode = \"smoke\" but [mcmc.smoke] is empty")
        for (k, v) in ov
            if k == "kernel_draws"
                push!(applied, @sprintf("kernel.n_draws %s -> %s", c["kernel"]["n_draws"], v))
                c["kernel"]["n_draws"] = v
            else
                haskey(c["mcmc"], k) || error("[mcmc.smoke] $k does not override anything in [mcmc]")
                push!(applied, @sprintf("mcmc.%s %s -> %s", k, c["mcmc"][k], v))
                c["mcmc"][k] = v
            end
        end
    end
    delete!(c["mcmc"], "smoke")
    # RE-CHECK AFTER THE OVERLAY, and BEFORE the loader writes its own keys.
    # check_values runs before the overlay, so every value the overlay
    # introduced was unvalidated: a review got nsamples = -5, nchains = 0,
    # kernel_draws = -3 and nsamples = "many" through. The ordering matters --
    # re-checking after setting mode_overrides made the "written by the loader,
    # do not supply it" rule fire on the loader's own output.
    isempty(applied) || check_values(c, tomlfile * " [after the mcmc.smoke overlay]")
    c["run"]["mode_overrides"] = applied

    # --- sites ------------------------------------------------------------
    raw, hdr = readdlm(c["paths"]["sites"], ',', String; header = true)
    icol = findfirst(==("Field_ID"), strip.(vec(hdr)))
    icol === nothing && error("[paths] sites: $(c["paths"]["sites"]) has no Field_ID column")
    have = Set(strip.(raw[:, icol]))
    missing_sites = [s for s in c["run"]["sites"] if !(s in have)]
    isempty(missing_sites) ||
        error("[run] sites not in $(basename(c["paths"]["sites"])): $(join(missing_sites, ", "))")
    isempty(c["run"]["sites"]) && error("[run] sites is empty")

    # --- the stages, in the only order that works -------------------------
    # detect must precede profile: run_aem_profile reads the DETECTED interval,
    # and with these reversed the profile silently skipped with "no detection
    # file" rather than failing.
    known = ["aem_occam", "aem_transd", "snmr_transd", "detect", "profile", "gate"]
    for s in c["run"]["stages"]
        s in known || error("[run] unknown stage \"$s\"; known: $(join(known, ", "))")
    end
    ix = [findfirst(==(s), known) for s in c["run"]["stages"]]
    issorted(ix) || error("[run] stages are out of order; required order: $(join(known, " -> "))")
    if "profile" in c["run"]["stages"] && !("detect" in c["run"]["stages"])
        error("[run] stage \"profile\" needs \"detect\" before it: it reads the detected interval")
    end

    c["run"]["config_file"] = normpath(abspath(tomlfile))
    CFG[] = c
    return c
end

"A resolved path from [paths]; errors rather than returning a wrong default."
function cfgpath(c::Dict, key::AbstractString)
    haskey(c["paths"], key) || error("[paths] has no $key")
    c["paths"][key]
end

is_smoke(c::Dict) = c["run"]["mode"] == "smoke"

"One block naming every parameter that could change a number. Print it per run."
function describe(c::Dict; io = stdout)
    @printf(io, "control file : %s\n", c["run"]["config_file"])
    @printf(io, "run          : %s   mode %s\n", c["run"]["name"], c["run"]["mode"])
    for o in c["run"]["mode_overrides"]
        @printf(io, "               override %s\n", o)
    end
    @printf(io, "sites        : %s\n", join(c["run"]["sites"], ", "))
    @printf(io, "stages       : %s\n", join(c["run"]["stages"], " -> "))
    @printf(io, "grid         : nz %d  zstart %.3g  dz %.3g  extendfrac %.4g\n",
            c["grid"]["nz"], c["grid"]["zstart"], c["grid"]["dz"], c["grid"]["extendfrac"])
    @printf(io, "mcmc         : %d samples, %d chains, tmax %.3g, seed_base %d\n",
            c["mcmc"]["nsamples"], c["mcmc"]["nchains"],
            get(c["mcmc"], "tmax", NaN), c["mcmc"]["seed_base"])
    @printf(io, "kernel draws : %d\n", c["kernel"]["n_draws"])
    @printf(io, "noise branch : %s\n", c["noise"]["sigma_hat_branch"])
    @printf(io, "detect       : theta %.3g  min_thick %.3g  metric %s\n",
            c["detect"]["theta"], get(c["detect"], "min_thick_m", NaN),
            get(c["detect"], "metric", "-"))
    @printf(io, "gate         : %s mg/L\n", join(c["gate"]["thresholds_mgL"], ", "))
    return nothing
end

end # module
