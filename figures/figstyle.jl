# figstyle.jl -- shared style for every paper figure.
#
# The brief's standing requirements, in one place so no figure drifts from them:
#   PNG at 300 dpi; 170 mm wide for two-column, 85 mm for one-column;
#   panel labels (a), (b), ... top-left of each panel; fonts >= 8 pt at print
#   size; NO suptitles, no editorial text, no truth values in any title --
#   truths are drawn as a line and named in the legend.
module FigStyle
using PyPlot, Printf
export MM170, MM85, DPI, C_FRESH, C_SALINE, C_TRUTH, C_TRUTH_ONBLACK, C_MGR, newfig, panel,
       aquifer_band, finish, OUTDIR

const OUTDIR = @__DIR__
const DPI   = 300
const MM170 = 170 / 25.4          # 6.693 in, two-column
const MM85  = 85 / 25.4           # 3.346 in, one-column

# Colours fixed by the brief.
const C_FRESH  = "#1f5fa9"        # fresh column, blue
const C_SALINE = "#e07b22"        # saline control, orange
const C_TRUTH  = "black"
# On the black-background density figures a black truth line is only its halo,
# and C_FRESH (#1f5fa9) is too dark against inferno's low end. This is the truth
# colour for those panels only; every white-background figure keeps black.
const C_TRUTH_ONBLACK = "#36c9ff"
const C_MGR    = "#8a8a8a"        # manager lines, grey dashed

function __init__()
    # 8 pt floor at print size. Everything scales from `base`.
    base = 8
    rc("font", size = base, family = "sans-serif")
    rc("axes", labelsize = base, titlesize = base, linewidth = 0.6)
    rc("xtick", labelsize = base - 0.5); rc("ytick", labelsize = base - 0.5)
    rc("xtick.major", width = 0.6, size = 2.5)
    rc("ytick.major", width = 0.6, size = 2.5)
    # frameon with a semi-opaque white face: a frameless legend drawn over
    # curves is unreadable, which it was on the first pass of Figure 2.
    rc("legend", fontsize = base - 1, frameon = true, handlelength = 1.6,
       borderaxespad = 0.3, labelspacing = 0.25, framealpha = 0.85,
       edgecolor = "none", facecolor = "white")
    rc("lines", linewidth = 1.0)
    rc("savefig", dpi = DPI, bbox = "tight", pad_inches = 0.02)
    rc("figure", dpi = 110)
    rc("mathtext", default = "regular")
end

"""A figure of the brief's width. `h` is height in inches."""
newfig(; width = MM170, h = 4.0, kw...) = subplots(; figsize = (width, h), kw...)

"""Panel label in the top-left, outside the data but inside the axes box."""
function panel(ax, letter; dx = 0.02, dy = 0.975, weight = "bold")
    # a white backing box: without it a label can land on a curve, which it did
    ax.text(dx, dy, "($letter)", transform = ax.transAxes, va = "top",
            ha = "left", fontweight = weight, fontsize = 8.5,
            zorder = 20,   # above density, wash and every overlay line
            bbox = Dict("facecolor" => "white", "edgecolor" => "none",
                        "alpha" => 0.85, "pad" => 0.8))
end

"""Shade an aquifer interval on a depth-on-y axis (or x if horiz=false)."""
function aquifer_band(ax, z1, z2; horiz = false, c = "#000000", a = 0.055, label = nothing)
    if horiz
        ax.axvspan(z1, z2; color = c, alpha = a, lw = 0,
                   label = label === nothing ? "" : label)
    else
        ax.axhspan(z1, z2; color = c, alpha = a, lw = 0,
                   label = label === nothing ? "" : label)
    end
end

"""Save as PNG at 300 dpi and report the pixel size, so the width can be checked."""
function finish(fig, name; width = MM170)
    p = joinpath(OUTDIR, name * ".png")
    fig.savefig(p, dpi = DPI)
    close(fig)
    mm = width * 25.4
    @printf("  wrote %-34s %.0f mm wide at %d dpi\n", name * ".png", mm, DPI)
    p
end
end # module
