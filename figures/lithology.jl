# lithology.jl -- the Synthetic Section Builder's symbology, exactly.
#
# Taken from tools/section_builder.html: the LITH table (line 214), LITH_ORDER,
# guessLith, and the dark-theme rule `mix(tint, "#121922", 0.52)`. Nothing here
# is invented; if the builder changes, this is what has to change with it.
#
# A block is drawn as the builder draws one: TINT, then the PATTERN over it in
# ink at 70% opacity, then ONE thin outline. Contacts are solid lines.
module Lithology

using PyPlot

export LITH_TINT, LITH_ORDER, guesslith, lithof, drawcolumn, lithkey

# tints, light-theme values, from the builder's LITH table
const LITH_TINT = Dict(
    "gravel"            => "#df8a2e",
    "clean sand"        => "#f3d152",
    "sand"              => "#e6bd3f",
    "silty sand"        => "#c6ae56",
    "sandy clay"        => "#8da667",
    "clay"              => "#6a93a6",
    "weathered bedrock" => "#c4738a",
    "saprock"           => "#a66a90",
    "fresh bedrock"     => "#8a6d9e",
)

const LITH_ORDER = ["gravel", "clean sand", "sand", "silty sand", "sandy clay",
                    "clay", "weathered bedrock", "saprock", "fresh bedrock"]

# Pattern ink, from the builder's --ink: #18212c light, #e4eaf1 dark. Figure 2
# is a LIGHT figure, so it takes the light ink and the UNMIXED tints. The 52%
# mix towards #121922 is the builder's DARK-theme rule and is applied only when
# a figure is drawn on the dark ground -- applying it on white made the tints
# muddy and the near-white pattern invisible.
# Pattern opacity. The builder uses 0.70; at figure scale the sand stipple was
# invisible, so the paper figures take 0.95 and a slightly heavier stroke.
const PATALPHA = 0.95
const INK_LIGHT = "#18212c"
const INK_DARK  = "#e4eaf1"
const DARKMIX = "#121922"      # what the builder mixes tints towards, dark only
const MIXT = 0.52

hex2rgb(h) = (parse(Int, h[2:3], base=16)/255,
              parse(Int, h[4:5], base=16)/255,
              parse(Int, h[6:7], base=16)/255)

"Mix `a` towards `b` by `t`, as the builder's mix() does."
function mix(a, b, t)
    ra, ga, ba = hex2rgb(a); rb, gb, bb = hex2rgb(b)
    (ra + (rb-ra)*t, ga + (gb-ga)*t, ba + (bb-ba)*t)
end

"The tint a block is filled with. `dark` applies the builder's 52% mix."
tintof(l; dark = false) = dark ? mix(get(LITH_TINT, l, LITH_TINT["sand"]), DARKMIX, MIXT) :
                                 hex2rgb(get(LITH_TINT, l, LITH_TINT["sand"]))

"""
The builder's guessLith, rule for rule (section_builder.html).
Used when a [[unit]] carries no `lithology` key -- which, on both paper columns,
is every unit.
"""
function guesslith(phi, clay)
    phi  < 0.05 && return "fresh bedrock"
    phi  < 0.09 && return "saprock"
    clay > 0.50 && return "clay"
    clay > 0.30 && return phi < 0.20 ? "weathered bedrock" : "sandy clay"
    clay > 0.13 && return "silty sand"
    clay > 0.05 && return "sand"
    return "clean sand"
end

lithof(u) = haskey(u, "lithology") && !isempty(String(u["lithology"])) ?
            String(u["lithology"]) :
            guesslith(Float64(get(u, "phi", 0.0)), Float64(get(u, "clay", 0.0)))

# ---- patterns, tile sizes as in the builder -----------------------------
# Each returns a list of primitives in TILE coordinates; the caller tiles them
# across the block. Sizes are the builder's: sand 12x14, clean sand 16x16,
# gravel 18x17, silty sand 16x16, sandy clay 16x16, clay 16x14, weathered 14x14,
# saprock 14x14, fresh bedrock 12x12.
const TILE = Dict(
    "sand" => (12.0, 14.0), "clean sand" => (16.0, 16.0), "gravel" => (18.0, 17.0),
    "silty sand" => (16.0, 16.0), "sandy clay" => (16.0, 16.0), "clay" => (16.0, 14.0),
    "weathered bedrock" => (14.0, 14.0), "saprock" => (14.0, 14.0),
    "fresh bedrock" => (12.0, 12.0))

"Primitives for one tile: (:dot,x,y,r), (:dash,x,y,len), (:ell,x,y,rx,ry), (:line,x1,y1,x2,y2)."
function tileprims(l)
    l == "sand"        && return [(:dot,3.0,3.5,0.95),(:dot,8.5,7.0,0.95),(:dot,5.0,11.0,0.95)]
    l == "clean sand"  && return [(:dot,4.5,5.0,0.8),(:dot,11.0,11.5,0.8)]
    l == "gravel"      && return [(:ell,5.0,5.5,3.2,2.1),(:ell,12.5,12.0,3.6,2.3),(:dot,13.5,4.0,1.0)]
    l == "silty sand"  && return [(:dot,4.0,4.5,0.9),(:dot,11.5,11.0,0.9),(:dash,3.0,11.5,5.0)]
    l == "sandy clay"  && return [(:dash,2.5,4.0,6.0),(:dash,8.0,11.0,6.0),(:dot,12.5,4.5,0.9)]
    l == "clay"        && return [(:dash,2.0,4.0,6.5),(:dash,8.5,10.0,6.5)]
    l == "weathered bedrock" && return [(:line,3.0,10.0,6.5,4.5),(:line,8.0,11.0,11.5,5.5)]
    l == "saprock"     && return [(:line,0.0,7.0,14.0,7.0),(:line,7.0,0.0,7.0,7.0),(:line,0.0,14.0,14.0,14.0),(:line,0.0,7.0,0.0,14.0)]
    l == "fresh bedrock" && return [(:line,0.0,0.0,12.0,12.0),(:line,12.0,0.0,0.0,12.0)]
    return []
end

"""
    drawcolumn(ax, units, zmax; width=1.0, dark=true)

One stratigraphic column. `units` is a vector of (top, base, lithology).
Tint, then pattern, then one thin outline; contacts as solid lines.
"""
function drawcolumn(ax, units, zmax; dark = false, px_per_m = 6.0)
    ink = dark ? INK_DARK : INK_LIGHT
    for (ztop, zbot, l) in units
        ztop >= zmax && continue
        zb = min(zbot, zmax)
        ax.add_patch(matplotlib.patches.Rectangle((0.0, ztop), 1.0, zb - ztop,
            facecolor = tintof(l; dark = dark), edgecolor = "none", zorder = 2))
        # pattern, tiled over the block in data coordinates
        tw, th = get(TILE, l, (14.0, 14.0))
        dx, dy = tw / 100.0, th / px_per_m
        prims = tileprims(l)
        # the tile grid overruns the block's right and lower edges by up to one
        # tile, so every primitive is clipped to the block itself. Without this
        # the pattern bleeds into whatever sits beside the column.
        clipr = matplotlib.patches.Rectangle((0.0, ztop), 1.0, zb - ztop,
                                             transform = ax.transData)
        clip!(a) = (a.set_clip_path(clipr); a)
        yy = ztop
        while yy < zb
            xx = 0.0
            while xx < 1.0
                for p in prims
                    if p[1] === :dot
                        clip!(ax.add_patch(matplotlib.patches.Circle(
                            (xx + p[2]/100, yy + p[3]/px_per_m), p[4]/100,
                            facecolor = ink, edgecolor = "none", alpha = PATALPHA, zorder = 3)))
                    elseif p[1] === :dash
                        clip!(ax.plot([xx + p[2]/100, xx + (p[2]+p[4])/100],
                                [yy + p[3]/px_per_m, yy + p[3]/px_per_m],
                                color = ink, lw = 0.5, alpha = PATALPHA, zorder = 3,
                                solid_capstyle = "butt")[1])
                    elseif p[1] === :ell
                        clip!(ax.add_patch(matplotlib.patches.Ellipse(
                            (xx + p[2]/100, yy + p[3]/px_per_m), 2p[4]/100, 2p[5]/px_per_m,
                            facecolor = "none", edgecolor = ink, lw = 0.5, alpha = PATALPHA, zorder = 3)))
                    else
                        clip!(ax.plot([xx + p[2]/100, xx + p[4]/100],
                                [yy + p[3]/px_per_m, yy + p[5]/px_per_m],
                                color = ink, lw = 0.6, alpha = PATALPHA, zorder = 3)[1])
                    end
                end
                xx += dx
            end
            yy += dy
        end
        # one thin outline, and the contact
        ax.add_patch(matplotlib.patches.Rectangle((0.0, ztop), 1.0, zb - ztop,
            facecolor = "none", edgecolor = ink, lw = 0.4, alpha = 0.55, zorder = 4))
    end
    ax.set_xlim(0, 1); ax.set_xticks([])
end

"The lithologies present, in the builder's order."
lithkey(units) = [l for l in LITH_ORDER if l in unique(x[3] for x in units)]

end # module
