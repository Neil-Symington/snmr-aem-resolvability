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

# The sites this run covers, honouring --site.
const RUN_SITES = isempty(ONLY) ? FCFG["run"]["sites"] :
                  (ONLY in FCFG["run"]["sites"] ? [ONLY] :
                   error("--site $ONLY is not in [run] sites: " *
                         join(FCFG["run"]["sites"], ", ")))
