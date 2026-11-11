# builder.jl — open the section builder in your browser.
#
#   julia --project=. tools/builder.jl
#
# The builder is one self-contained HTML file next to this script. It needs no
# server and no network: the only external reference is the Google Fonts
# stylesheet, which falls back to a system stack when you are offline.
#
# Build a section by dragging contacts and setting porosity, saturation, clay
# fraction and pore-water conductivity, then save the scenario TOML and feed it
# to build_truth.jl:
#
#   julia --project=. synthetics/scripts/build_truth.jl <saved>.toml
#
# The forward models it previews are the repository's own, from
# synthetics/src/SNMRAEMSynthetics.jl, and the lithology presets and clay
# surface conductivities are from synthetics/src/Lithology.jl. What it draws is
# piecewise constant with sharp contacts, while build_truth.jl adds
# depth-correlated variability and gradational boundaries -- so use it to judge
# contrast, not detail. The page says so too.
const PAGE = joinpath(@__DIR__, "section_builder.html")
isfile(PAGE) || error("missing $PAGE")

# A file:// URL, so the page's own "Save .toml" button works. That button is
# hidden when the same file is served over http, because a hosted copy runs in a
# sandbox that blocks page-initiated downloads and the button would do nothing.
url = "file://" * PAGE

opener = if Sys.islinux(); "xdg-open"
         elseif Sys.isapple(); "open"
         elseif Sys.iswindows(); "cmd"
         else ""; end

println("section builder: ", PAGE)
if isempty(opener) || Sys.which(opener) === nothing
    println("\nCould not find a browser opener on this system. Open this yourself:\n\n  ", url, "\n")
else
    try
        Sys.iswindows() ? run(`cmd /c start "" $url`) : run(`$opener $url`)
        println("opened in your default browser")
    catch e
        println("\n", opener, " failed (", e, "). Open this yourself:\n\n  ", url, "\n")
    end
end
