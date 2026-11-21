# window.jl -- the ONE window-integral convention, shared by every stage.
#
# THE SNMR MODEL IS PIECEWISE-CONSTANT OVER 60 LAYERS whose BOTTOMS are `zsmr`.
# The first layer runs 0 to zsmr[1] = 1.0 m. `SNMRForward.MRSForward_square` is
# handed that same vector as the layer specification, and a posterior draw has
# one theta per layer, so a window integral must weight each layer by its
# OVERLAP with the window -- not interpolate between the zsmr values.
#
# This existed only inside tools/stage2_report.jl until 22 Nov 2026. Reading
# `zsmr` as 60 cell BOUNDARIES instead -- 59 intervals, theta linear between
# them, the 0-1.0 m layer dropped -- gives 0.2786 on 3P where the right answer
# is 0.2811, and 0.2747 on 2P where it is 0.2672. The error is small enough to
# look like a run-to-run difference and large enough to move a ceiling, which
# is why the convention now lives in one place with a test on known values.
module Window
export window_integral, window_mean, window_mean_draws

"""
    window_integral(θ, zsmr, a, b) -> (∫θ dz, support)

`θ` is one value per layer; layer i spans (`zsmr[i-1]`, `zsmr[i]`], with
`zsmr[0]` taken as 0. Both returns are in metres: the second is the window's
support, which equals `b - a` whenever the window lies inside the grid.
"""
function window_integral(θ::AbstractVector, zsmr::AbstractVector, a::Real, b::Real)
    length(θ) == length(zsmr) ||
        throw(DimensionMismatch("theta has $(length(θ)) layers, zsmr has $(length(zsmr)) " *
                                "bottoms; they must match -- zsmr is LAYER BOTTOMS, not boundaries"))
    top = vcat(0.0, zsmr[1:end-1])
    s = 0.0; w = 0.0
    for i in eachindex(zsmr)
        ov = max(0.0, min(zsmr[i], b) - max(top[i], a))
        s += θ[i] * ov; w += ov
    end
    (s, w)
end

"Mean theta over the window: the integral divided by the support."
function window_mean(θ::AbstractVector, zsmr::AbstractVector, a::Real, b::Real)
    s, w = window_integral(θ, zsmr, a, b)
    w > 0 || throw(ArgumentError("window $a-$b has no support on this grid"))
    s / w
end

"Mean theta over the window for every column of `Θ` (layers x draws)."
window_mean_draws(Θ::AbstractMatrix, zsmr::AbstractVector, a::Real, b::Real) =
    [window_mean(view(Θ, :, k), zsmr, a, b) for k in axes(Θ, 2)]
end
