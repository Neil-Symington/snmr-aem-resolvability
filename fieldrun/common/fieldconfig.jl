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

"""
    load_config(tomlfile) -> Dict

Read, validate and normalise a control file. Errors name the offending key.
"""
function load_config(tomlfile::AbstractString)
    isfile(tomlfile) || error("no control file at $tomlfile")
    c = TOML.parsefile(tomlfile)

    for ks in REQUIRED
        getin(c, ks) === nothing && error("control file $tomlfile is missing [$(ks[1])] $(ks[2])")
    end

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
