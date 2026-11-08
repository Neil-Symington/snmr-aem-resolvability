# stage_prelude.jl -- what every field stage must do before anything else.
#
# The stages were written as standalone scripts with their parameters in four
# places. Ported into the repository they all need the same two things first,
# and in this order:
#
#   1. the control file loaded, because stages/udf_system.jl parses the SkyTEM
#      .des AT LOAD TIME and reads its path from FieldConfig.CFG[]. Including a
#      support module before the config is loaded gives
#      `UndefVarError: FieldConfig not defined` from inside udf_system.jl, which
#      names neither the cause nor the fix.
#   2. the resolved output and log directories, so a stage writes where the
#      control file says rather than beside itself.
#
# Each stage includes this first, then its own support modules.
const COMMON = normpath(@__DIR__)
const FIELDRUN = normpath(joinpath(@__DIR__, ".."))

include(joinpath(COMMON, "fieldconfig.jl"))
using .FieldConfig
using Printf

"The --config argument, defaulting to the shipped two-site example."
const CFGPATH = let i = findfirst(==("--config"), ARGS)
    i === nothing ? joinpath(FIELDRUN, "config", "udf_example.toml") : ARGS[i+1]
end

const FCFG = FieldConfig.load_config(CFGPATH)

"--site NAME restricts a per-site stage to one site, so several can run as
independent processes; they share nothing but the manifest."
const ONLY = let i = findfirst(==("--site"), ARGS)
    i === nothing ? "" : ARGS[i+1]
end

const OUT  = FieldConfig.cfgpath(FCFG, "outputs")
const LOGS = haskey(FCFG["paths"], "logs") ? FieldConfig.cfgpath(FCFG, "logs") : OUT
mkpath(OUT)

# THE MCMC SETTINGS, from [mcmc]. The stages were written with their own argv
# defaults -- 200,001 samples and 4 chains for the AEM trans-D, 150,001 for the
# SNMR -- and ignored the control file entirely, so `mode = "smoke"` did nothing
# for the two stages that take hours. The first two-site run was launched at
# production sample counts without meaning to be.
#
# That is the exact failure a control file exists to prevent: the same parameter
# in two places, and the one a reader edits being the one that is ignored. An
# explicit --nsamples or --nchains still wins, for a one-off, but the default is
# now what the config says.
const MC = FCFG["mcmc"]
const ARGN = let d = Dict{String,String}(), i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) &&
            (d[ARGS[i][3:end]] = ARGS[i+1]; i += 1)
        i += 1
    end
    d
end
const NSAMPLES_CFG = haskey(ARGN, "nsamples") ? parse(Int, ARGN["nsamples"]) : Int(MC["nsamples"])
const NCHAINS_CFG  = haskey(ARGN, "nchains")  ? parse(Int, ARGN["nchains"])  : Int(MC["nchains"])
const TMAX_CFG     = Float64(get(MC, "tmax", 2.5))

@printf("[%s] %d samples, %d chains, Tmax %.3g  (mode %s)\n",
        basename(CFGPATH), NSAMPLES_CFG, NCHAINS_CFG, TMAX_CFG, FCFG["run"]["mode"])

# The sites this run covers, honouring --site.
const RUN_SITES = isempty(ONLY) ? FCFG["run"]["sites"] :
                  (ONLY in FCFG["run"]["sites"] ? [ONLY] :
                   error("--site $ONLY is not in [run] sites: " *
                         join(FCFG["run"]["sites"], ", ")))
