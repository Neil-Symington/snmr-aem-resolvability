#!/usr/bin/env julia
# Build a synthetic column from a TOML scenario and write it to JLD2.
#
#   julia --project=. synthetics/scripts/build_truth.jl synthetics/config/confined_aquifer.toml
#
# Writes output/<name>_truth.jld2 and prints the section as built.
using Pkg
isfile(joinpath(@__DIR__, "..", "..", "Project.toml")) || error("run from the repo root with --project=.")
include(joinpath(@__DIR__, "..", "src", "SNMRAEMSynthetics.jl"))
using .SNMRAEMSynthetics, JLD2

length(ARGS) >= 1 || error("usage: build_truth.jl <config.toml> [outdir]")
cfgpath = ARGS[1]
outdir  = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "..", "output")
mkpath(outdir)

units, cfg = load_config(cfgpath)
col = build_column(units, cfg)
show(stdout, col)
check_column(col)

name = get(cfg, "name", splitext(basename(cfgpath))[1])
out  = joinpath(outdir, "$(name)_truth.jld2")

# The SNMR grid comes from the inversion's own grid definition, so a column
# built here is a drop-in for the stages. See smr_acquisition's docstring for
# why it is this grid and not [smr]'s uniform one.
include(joinpath(@__DIR__, "..", "..", "src", "prerun_config.jl"))
zsmr = PrerunConfig.grid()[3]        # cell BOUNDARIES, 1.0 -> 147.438 m

save_column(col, out; zsmr = zsmr)
println("\nwrote ", out)
