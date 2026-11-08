# figstyle.jl -- shared style for every paper figure.
#
# The brief's standing requirements, in one place so no figure drifts from them:
#   PNG at 300 dpi; 170 mm wide for two-column, 85 mm for one-column;
#   panel labels (a), (b), ... top-left of each panel; fonts >= 8 pt at print
#   size; NO suptitles, no editorial text, no truth values in any title --
#   truths are drawn as a line and named in the legend.
module FigStyle
using PyPlot, Printf
# JLD2 directly, not `@eval Main using JLD2` + invokelatest. That dance cannot
# see a binding created in the same world age, so assert_column threw
# UndefVarError: `load` whenever it was called from anywhere that had not
# already loaded JLD2 at top level. The same trap has now cost time three times
# in this repository; a plain import has none of it, and JLD2 is a direct
# project dependency.
using JLD2
export MM170, MM85, DPI, C_FRESH, C_SALINE, C_TRUTH, C_TRUTH_ONBLACK, C_MGR, newfig, panel,
       aquifer_band, finish, OUTDIR, assert_column, provenance, archie_m

const OUTDIR = @__DIR__
const ROOT   = normpath(joinpath(@__DIR__, ".."))

"""
    assert_column(tag)

Refuse to draw a panel labelled `tag` from another column's data.

NOTHING CHECKED THIS. A reviewer copied the 2P column's files under the names
3P, 2S and 3S, ran fig05, fig07 and figS8, and got PNGs and number CSVs with
rows labelled "3P" and "2S" side by side, identical to the digit, with no
warning. Every figure took the tag from its own hard-coded list and never
compared it with the `column` the Stage 0 artefact records. A figure in the
paper could be built from the wrong column and look entirely normal.

Called for every tag a figure draws. If Stage 0 has not run for that column
there is nothing to check against, and that is said rather than passed over.
"""
function assert_column(tag::AbstractString)
    a = joinpath(ROOT, "outputs", "synth_$tag", "stage0_synth_$tag.jld2")
    if !isfile(a)
        # ERROR, not a warning. The first version warned and drew anyway, and a
        # reviewer pointed out that six scripts then produce a figure from
        # unverified data -- which for a paper figure is the wrong default. The
        # escape hatch is explicit so that choosing it is a decision.
        haskey(ENV, "FIGS_ALLOW_UNVERIFIED") || error("""
    NO STAGE 0 ARTEFACT FOR $tag, so its column label cannot be verified.

      expected: $a

      Run Stage 0 for this column, or set FIGS_ALLOW_UNVERIFIED=1 to draw
      anyway. Do not use an unverified figure in the paper: a figure drawn from
      the wrong column looks entirely normal.
    """)
        @warn "FIGS_ALLOW_UNVERIFIED is set: $tag is NOT verified" path=a
        return nothing
    end
    A = JLD2.load(a)
    got = get(A, "column", nothing)
    got === nothing && error("$a records no `column`; cannot verify the label $tag")
    String(got) == tag || error("""
    FIGURE LABEL DOES NOT MATCH THE DATA.

      This panel is labelled $tag, but $a
      was sealed from column $(got).

      A figure drawn from the wrong column looks entirely normal, which is why
      this is an error and not a warning.
    """)

    # AND THE TRUTH FILE BEHIND IT. `column` is copied from the config, so a
    # config declaring column = "3P" while pointing at truth_murray2P.jld2 is
    # sealed by the real Stage 0 and passes the check above -- a reviewer built
    # exactly that and drew fig04 and fig05 from it.
    #
    # Stage 0 records truth_file_sha, so compare it against the hash of the
    # truth this column is SUPPOSED to come from. That is the naming convention
    # every committed config follows, and it makes the check an identity rather
    # than a filename match. Recording the path instead was the first attempt
    # and Stage 0's key whitelist refused it, correctly.
    want = joinpath(ROOT, "data", "truth_murray$tag.jld2")
    sha  = get(A, "truth_file_sha", "")
    if isfile(want) && !isempty(string(sha))
        mine = bytes2hex(reinterpret(UInt8, [hash(read(want))]))
        mine == string(sha) || error("""
    THIS COLUMN WAS NOT SEALED FROM ITS OWN TRUTH FILE.

      panel label   : $tag
      expected truth: $(relpath(want, ROOT))
      sealed sha    : $(sha)
      that file's   : $(mine)

      $a records column $tag but was sealed from a different truth file. A
      config may declare any column it likes; this catches one that declares
      $tag while pointing somewhere else.
    """)
    end
    nothing
end

"""
    archie_m(T)

Archie's cementation exponent from a truth file, under either spelling.

The paper's truth files write `M_ARCHIE`; the shipped generator writes
`m_archie`. Three figures read `M_ARCHIE` bare, so on a GENERATOR-BUILT column
they died with `KeyError: key "M_ARCHIE" not found` -- a reviewer hit it on
fig02, fig06 and export_numbers_v2 after following the documented chain. The
two spellings are the same quantity.
"""
function archie_m(T)
    haskey(T, "M_ARCHIE") && return T["M_ARCHIE"]
    haskey(T, "m_archie") && return T["m_archie"]
    b = get(T, "build", nothing)
    b !== nothing && haskey(b, "m_archie") && return b["m_archie"]
    error("""
    the truth file carries no Archie exponent.

      Looked for "M_ARCHIE" (the paper's truth files), "m_archie" (what
      synthetics/scripts/build_truth.jl writes) and build["m_archie"].
    """)
end

"""
    provenance(tags)

One line naming the columns a figure was drawn from and the config that sealed
each, printed so a figure's origin is in the build log. A reviewer noted that a
figure built from a SMOKE run carries no marker distinguishing it from a
production one; this at least records which config it came from.
"""
function provenance(tags)
    for t in tags
        a = joinpath(ROOT, "outputs", "synth_$t", "stage0_synth_$t.jld2")
        isfile(a) || continue
        A = JLD2.load(a)
        @printf("  %-4s column %-4s  config %s\n", t, get(A, "column", "?"),
                basename(string(get(A, "config_source", "?"))))
    end
end
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
