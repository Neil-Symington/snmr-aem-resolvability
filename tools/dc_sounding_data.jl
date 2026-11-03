# dc_sounding_data.jl -- Stage 0. Generate the synthetic Schlumberger soundings.
#
#   julia --project=. dc_sounding_data.jl
#
# STAGE 0 WORK: this is the only place a truth column is read for the sounding.
# Writes dc_sounding_<tag>.jld2 holding AB/2, MN/2, the clean and noisy
# apparent resistivities and the per-point sigma. Nothing downstream reads the
# truth.
#
# ARRAY. AB/2 at 20 log-spaced values from 2 to 300 m, MN/2 = AB/10. The finite
# MN matters: at (l/L) = 0.2 the dipole-limit kernel is wrong by up to 8% on a
# 10:1 contrast, which is larger than the 3% noise, so the finite-MN form is
# used for both the data and the inversion.
#
# NOISE. sigma_i = rho_a_i * sqrt(0.03^2 + 0.01^2) = 3.162% of the datum --
# 3% Gaussian in quadrature with a 1% RELATIVE floor. Being relative, the floor
# scales with the datum and so never binds; it raises the effective error from
# 3.00% to 3.16% and does nothing else. It replaces an absolute 0.5 Ohm.m floor
# that was 5-25% of the datum on the conductive columns. Same seed as that run,
# so the standard normal draws are identical and only sigma has changed.
using JLD2, Random, Printf, Statistics
include(joinpath(@__DIR__, "..", "src", "dc_forward.jl")); using .DCForward
const OUT = joinpath(@__DIR__, "..", "data")
const SEED = 20261103          # recorded in config.toml as [data.dc].seed
const NAB  = 20
const ABLO, ABHI = 2.0, 300.0
const PCT, FLOOR_REL = 0.03, 0.01

"""Layer resistivities and thicknesses from a truth column."""
function truth_layers(T)
    z = T["z"]; dz = T["dz"]
    rho = 1.0 ./ T["σ_bulk"]
    h = fill(dz, length(rho)-1)
    rho, h
end

function main()
    L = 10 .^ range(log10(ABLO), log10(ABHI), length = NAB)
    l = L ./ 5                                  # MN/2 = AB/10 => l = L/5
    for tag in ("2P","2S","3P","3S")
        f = joinpath(OUT, "truth_murray$tag.jld2")
        isfile(f) || (@printf("  skip %s (no truth)\n", tag); continue)
        T = load(f); rho, h = truth_layers(T)
        clean = rho_schlumberger_finite(rho, h, L, l)
        sd = clean .* sqrt(PCT^2 + FLOOR_REL^2)
        rng = MersenneTwister(SEED + parse(Int, tag[1]) + (tag[2] == 'S' ? 100 : 0))
        obs = clean .+ sd .* randn(rng, NAB)
        jldsave(joinpath(OUT, "dc_sounding_$tag.jld2");
                AB_2 = collect(L), MN_2 = collect(l), rho_a_clean = clean,
                rho_a = obs, sigma = sd, seed = SEED, pct = PCT, floor_rel = FLOOR_REL,
                column = tag, array = "Schlumberger, finite MN",
                forward = "dc_forward.jl, Ghosh J0 filter, verified to 3.3e-5 on two layers")
        @printf("%s  rho_a %.2f - %.2f Ohm.m   sigma %.3f - %.3f   chi2/N of the truth %.4f\n",
                tag, minimum(clean), maximum(clean), minimum(sd), maximum(sd),
                mean(((obs .- clean) ./ sd).^2))
    end
    @printf("\nAB/2 (m): %s\n", join((@sprintf("%.1f", x) for x in L), " "))
    println("wrote dc_sounding_<tag>.jld2 for each column")
end
main()
