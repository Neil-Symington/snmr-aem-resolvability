#!/usr/bin/env julia
# run_site.jl -- run the Section 3 workflow for one site.
#
#   julia --project=. run_site.jl configs/synth_2P.toml
#   julia --project=. run_site.jl configs/synth_2P.toml --from 3 --to 6
#   julia --project=. run_site.jl configs/synth_2P.toml --dry-run
#
# Stages, in order:
#   0  inputs      read the truth (synthetic) or the field data, seal them
#   1  AEM         trans-dimensional conductivity posterior          [hours]
#   2  SNMR        water-content posterior on the Stage 1 kernel     [hours]
#   3  hypothesis  the pre-registered interval and water criterion
#   4  gate        the admissible conductance band, and the floor
#   5  quality     pore-water salinity under the Stage 0 priors
#   6  verdict     the table
#
# STAGES 1 AND 2 ARE LONG. Each is a trans-dimensional MCMC of 100k-200k
# samples over four chains: hours on a workstation, not minutes. They are run
# through tools/, not through this driver, because they want a machine you are
# willing to leave alone and because their outputs are reused across configs.
# `run_site.jl` checks for their products and tells you what to run if they are
# missing, rather than silently starting a job that will outlive your terminal.
#
# EVERYTHING FROM STAGE 3 ON IS MINUTES, and that is the part most people want
# to re-run: it is where the hypothesis, the priors and the thresholds live.
using Printf
const HERE = @__DIR__
include(joinpath(HERE, "src", "config.jl")); using .SiteConfig

function parseargs(args)
    isempty(args) && error("usage: run_site.jl <config.toml> [--from N] [--to N] [--dry-run]")
    cfgpath = args[1]
    opt = Dict{String,Any}("from" => 0, "to" => 6, "dry" => false)
    i = 2
    while i <= length(args)
        a = args[i]
        if a == "--dry-run"; opt["dry"] = true; i += 1
        elseif a in ("--from","--to") && i < length(args)
            opt[a[3:end]] = parse(Int, args[i+1]); i += 2
        else; error("unknown argument: $a"); end
    end
    cfgpath, opt
end

function main()
    cfgpath, opt = parseargs(ARGS)
    cfg = SiteConfig.resolve(cfgpath)
    # A dry run validates the config but does not insist the generated inputs
    # are present -- that is the question a dry run is usually asked.
    info = SiteConfig.assertvalid(cfg; root = HERE, require_files = !opt["dry"])
    name = cfg["name"]; odir = joinpath(HERE, cfg["output"]["dir"])
    opt["dry"] || mkpath(odir)
    @printf("site %s  (%s)\n", name, cfg["label"])
    @printf("  column %s, role %s%s\n", cfg["data"]["column"], info.role,
            isempty(info.parent) ? "" : ", from parent " * info.parent)
    @printf("  output -> %s\n\n", cfg["output"]["dir"])

    # Stages 1 and 2 are not driven from here; report whether their products exist.
    if opt["from"] <= 3
        s1 = joinpath(HERE, cfg["data"]["stage1_draws"]); s2 = cfg["data"]["snmr_ensemble"]
        isfile(s1) || @printf("""
  STAGE 1 INPUT MISSING: %s
    run  julia --project=. tools/aem_HYP.jl 200001 <tag>
    then julia --project=. tools/aem_draws.jl <tag>

""", cfg["data"]["stage1_draws"])
        (isempty(s2) || isfile(joinpath(HERE, s2))) || @printf("""
  STAGE 2 INPUT MISSING: %s
    run  julia --project=. stages/stage2_snmr.jl --tier typical --conds C2 --tag <tag>

""", s2)
    end

    steps = [(0, "stages/stage0_inputs.jl"),
             (3, "stages/stage3_hypothesis.jl"),
             (4, "stages/stage4_gate_v2.jl"),
             (4, "stages/stage4_post.jl"),
             (5, "stages/stage5_quality.jl"),
             (6, "stages/stage6_verdict.jl")]
    for (st, script) in steps
        (st < opt["from"] || st > opt["to"]) && continue
        @printf("-- stage %d  %s\n", st, basename(script))
        if opt["dry"]; println("     (dry run)"); continue; end
        run(`$(Base.julia_cmd()) --project=$HERE $(joinpath(HERE, script)) --config $cfgpath`)
    end
    if opt["to"] >= 4 && !opt["dry"]
        @printf("-- stage 4  floor under stated ceilings\n")
        run(`$(Base.julia_cmd()) --project=$HERE $(joinpath(HERE,"stages","stage4_floor.jl")) $(cfg["data"]["column"])`)
    end
    opt["dry"] ? println("\ndry run: config valid, nothing executed") :
                 println("\ndone. verdict table: ", joinpath(cfg["output"]["dir"], "verdict_$name.md"))
end
main()
