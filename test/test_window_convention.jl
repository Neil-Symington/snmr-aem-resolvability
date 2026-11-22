# test_window_convention.jl -- the one window convention, on known values.
#
# The SNMR model is piecewise-constant over 60 LAYERS whose bottoms are zsmr,
# the first running 0 to 1.0 m. Reading zsmr as 60 cell BOUNDARIES instead --
# 59 layers, theta interpolated between them, the 0-1.0 m layer dropped -- moves
# the 3P window mean from 0.2811 to 0.2786 and the 2P one from 0.2672 to 0.2747.
# That is small enough to pass for run-to-run spread and large enough to move a
# salinity ceiling, so the convention is pinned here to published numbers.
using Test, JLD2, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "window.jl")); using .Window

const RUN3P = joinpath(ROOT, "data", "stage2", "3PT_pooled_20261118-145749",
                       "cond_C2_typical3PT.jld2")
const TRU3P = joinpath(ROOT, "data", "truth_murray3P.jld2")

@testset "window convention" begin
    @testset "arithmetic" begin
        z = [1.0, 2.0, 3.0]              # layers 0-1, 1-2, 2-3
        θ = [1.0, 2.0, 4.0]
        @test window_integral(θ, z, 0.0, 3.0) == (7.0, 3.0)      # 1+2+4
        @test window_mean(θ, z, 0.0, 3.0) ≈ 7/3
        @test window_mean(θ, z, 1.0, 3.0) ≈ 3.0                  # (2+4)/2
        @test window_mean(θ, z, 1.5, 2.5) ≈ 3.0                  # half of each
        # the FIRST layer runs from 0, not from z[1]
        @test window_integral(θ, z, 0.0, 1.0) == (1.0, 1.0)
        # a mismatched length is an error, not a silent truncation
        @test_throws DimensionMismatch window_integral([1.0, 2.0], z, 0.0, 3.0)
    end

    if isfile(RUN3P) && isfile(TRU3P)
        @testset "published 3P values" begin
            d = load(RUN3P); z = collect(Float64, d["zsmr"])
            W = vcat([window_mean_draws(Float64.(M), z, 7.5, 23.5) for M in d["ens"]]...)
            @test length(W) == 40016
            @test isapprox(quantile(W, 0.05), 0.2760; atol = 5e-5)
            @test isapprox(quantile(W, 0.50), 0.2811; atol = 5e-5)
            @test isapprox(quantile(W, 0.95), 0.2893; atol = 5e-5)

            T = load(TRU3P)
            wt = window_mean(collect(Float64, T["θ_smrgrid"]),
                             collect(Float64, T["zsmr"]), 7.5, 23.5)
            @test isapprox(wt, 0.2704; atol = 5e-5)
            # the support of a window inside the grid is its own length
            @test window_integral(collect(Float64, T["θ_smrgrid"]),
                                  collect(Float64, T["zsmr"]), 7.5, 23.5)[2] ≈ 16.0
        end
    else
        @info "skipped: the 3P Stage 2 run of record is not on disk" run = RUN3P
    end
end

@testset "every bare flag is in FLAGS" begin
    # A bare flag missing from FLAGS falls through to the "--key value" branch.
    # If it is the LAST argument that branch is skipped and the flag vanishes
    # without a word: --truekernel did exactly that on 22 Nov 2026 and both
    # control runs came back as duplicate pooled runs.
    src = read(joinpath(ROOT, "stages", "stage2_snmr.jl"), String)
    m = match(r"const FLAGS = \(([^)]*)\)", src)
    @test m !== nothing
    listed = Set(strip.(replace.(split(m[1], ","), '"' => ""), Ref(' ')))
    # every --flag the file tests for with haskey(ARG, ...) and no value
    for f in ("truekernel", "force")
        @test "--$f" in listed
    end
    # and the file really does branch on them
    @test occursin("haskey(ARG, \"truekernel\")", src)
end
