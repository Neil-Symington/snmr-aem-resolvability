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
        # These used to print partial commands -- aem_HYP.jl without its chain
        # count or label, Stage 2 without --label or the sample counts, and no
        # mention of the cp or of smr_data.jl at all -- so a reader following
        # them ran Stage 2 at its 100,001-sample default for hours. A review
        # had to reconstruct the chain twice. The whole chain, once:
        chain = """
  THE CHAIN, for a tag <T> and an optional run label <L>:

    julia --project=. synthetics/scripts/build_truth.jl synthetics/config/confined_aquifer.toml
    mkdir -p data
    cp synthetics/output/<name>_truth.jld2 data/truth_murray<T>.jld2
    julia --project=. tools/aem_HYP.jl 200001 <T> 4 <L>
    julia --project=. tools/aem_draws.jl <T> 600 <L>
    julia --project=. tools/smr_data.jl transd --tag <T> --label <L>
    julia --project=. stages/stage2_snmr.jl --tier typical --conds C2 --tag <T>T --label <L>

  For a SMOKE run use configs/smoke_2P.toml, which ships for this purpose. Note
  that aem_HYP.jl and aem_draws.jl take POSITIONAL arguments, not flags -- an
  earlier version of this hint said to add --nsamples to them, and a reviewer
  found that it makes parse(Int, "--nsamples") fail. Only stage2_snmr.jl takes
  flags. The smoke forms are:

    julia --project=. tools/aem_HYP.jl 2001 <T> 2 <L>
    julia --project=. tools/aem_draws.jl <T> 50 <L>
    julia --project=. stages/stage2_snmr.jl --tier typical --conds C2 --tag <T>T --label <L> --nsamples 2001

  Stage 2's default is 100,001 samples and takes hours. It runs four chains and
  four ensembles whatever you ask: both are constants in that stage and it reads
  no flag for either. This hint printed --nchains 2 until 12 Nov 2026, which the
  stage accepted and ignored.

  Or run the whole chain with ./run_smoke.sh or ./run_production.sh, each of
  which edits at the top.
"""
        isfile(s1) || @printf("\n  STAGE 1 INPUT MISSING: %s\n%s",
                              cfg["data"]["stage1_draws"], chain)
        (isempty(s2) || isfile(joinpath(HERE, s2))) ||
            @printf("\n  STAGE 2 INPUT MISSING: %s\n%s", s2, chain)
    end

    # Stage 4 reads Stage 0's sealed artefact, so --from 3 and --from 4 only work
    # once Stage 0 has run for this site. The README called --from 3 "the usual
    # entry point" without saying so, and on a clean tree it fails with an
    # ArgumentError from inside Stage 4 that names a path and not a fix.
    if opt["from"] > 0 && !opt["dry"]
        s0 = joinpath(HERE, cfg["output"]["dir"], "stage0_$(name).jld2")
        isfile(s0) || error(
            "Stage 0 has not been run for this site: " * s0 * " is missing.\n\n" *
            "  Stage 4 reads that artefact, so --from " * string(opt["from"]) *
            " cannot work yet. Run Stage 0\n" *
            "  once first -- it takes seconds and re-seals the inputs:\n\n" *
            "    julia --project=. run_site.jl " * cfgpath * " --to 0\n\n" *
            "  then --from " * string(opt["from"]) * " works for every later re-run.")
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
    vt = joinpath(cfg["output"]["dir"], "verdict_$name.md")
    if opt["dry"]
        println("\ndry run: config valid, nothing executed")
    elseif opt["to"] >= 6 && isfile(joinpath(HERE, vt))
        println("\ndone. verdict table: ", vt)
    else
        # --to 0 used to announce a verdict table it had not written, which a
        # review caught. Say what actually happened.
        @printf("\ndone, stages %d-%d. the verdict table needs stage 6.\n",
                opt["from"], opt["to"])
    end
end
main()
