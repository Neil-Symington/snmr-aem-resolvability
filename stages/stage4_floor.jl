# stage4_floor.jl -- Stage 4 addition. A salinity FLOOR under stated ceilings.
#
#   julia --project=. stage4_floor.jl [tag ...]
#
# POST-PROCESSING ONLY. Reads finished stage outputs; runs no inversion and
# changes none.
#
# The gate converts the band's UPPER edge clay-free, which overstates salinity
# and so can only PASS or leave the question open. This is the mirror: cap
# everything that could be explaining the measured conductance instead of
# salt, and whatever conductance is left must be water.
#
#   sigma_bulk_min = C_lo / L                       the data's lowest conductance
#   clay ceiling   = sigma_clay_max * C_max(theta)  C_max = c0 (1 - theta/phi_ref)
#   sigma_w floor  = (sigma_bulk_min - clay ceiling) / porosity_max^m, >= 0
#
# AN OPEN LOWER EDGE GIVES A FLOOR OF ZERO, by construction and not by
# arithmetic: an open edge is where the solver's sigma_min sits, not where the
# data stop, so there is no measured minimum conductance to subtract from.
# The pieces are still reported, because they are what the field sites will
# need and because a reader should see why the answer is zero.
#
# m: sigma_w = X / phi^m with phi < 1, so phi^m falls as m rises and sigma_w
# RISES with m. The floor is therefore at the m prior's LOWER bound. The
# brief's parenthetical said upper; both are printed and the lower one is the
# floor. [ceilings].m_for_floor records the choice.
using JLD2, Statistics, Printf, TOML
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE,".."))
include(joinpath(ROOT,"src","config.jl")); using .SiteConfig
const SYN = joinpath(HERE, "..", "data")
const OUT = joinpath(ROOT,"outputs")
const MGL = 6400.0
const ENS = Dict("3P"=>"3PT","2P"=>"2PT","3S"=>"3ST","2S"=>"2ST",
                 "2PJ"=>"2PJT","2SJ"=>"2SJT")

cellsof(z)=(mid=(z[1:end-1].+z[2:end])./2; (top=vcat(0.0,mid), bot=vcat(mid,z[end]+(z[end]-z[end-1])/2)))
function winw(z,a,b); g=cellsof(z); w=[max(0.0,min(g.bot[i],b)-max(g.top[i],a)) for i in eachindex(z)]
    ix=findall(>(1e-12),w); (ix=ix,w=w[ix]); end

"""Stage 2 median water content over the window: the window-mean theta of each
draw, then the median over draws."""
function theta_snmr(ensfile, a, b)
    E = load(ensfile); Θ = Matrix{Float64}(reduce(hcat, E["ens"]))
    ww = winw(E["zsmr"], a, b); L = b - a
    median([sum(Θ[ww.ix,j] .* ww.w)/L for j in 1:size(Θ,2)])
end

function floorfor(tag; cfgpath = joinpath(ROOT,"configs","synth_$tag.toml"))
    cfg = SiteConfig.resolve(cfgpath)
    ce  = cfg["ceilings"]; pr = cfg["priors"]
    mlo, mhi = SiteConfig.assertvalid(cfg; root=ROOT, require_files=false).m_prior
    P4 = load(joinpath(OUT,"synth_$tag","stage4post_synth_$tag.jld2"))
    w  = (P4["window"].top, P4["window"].base); L = w[2]-w[1]
    Clo = P4["band_C"].lo; openlo = P4["band_open"].lo
    θ   = theta_snmr(joinpath(SYN,"cond_C2_typical$(ENS[tag]).jld2"), w...)
    Cmax = ce["clay_c0"] * (1 - θ/ce["clay_phi_ref"])
    clayceil = ce["sigma_clay_max"] * Cmax
    φmax = ce["porosity_max"]
    lead = Clo/L - clayceil
    σw(m) = max(0.0, lead) / φmax^m
    σw_lo, σw_hi = σw(mlo), σw(mhi)              # mlo is the floor
    floor_tds = openlo ? 0.0 : σw_lo*MGL
    Xlim = Float64(cfg["hypothesis"]["salinity"]["value"])
    (tag=tag, window=w, L=L, C_lo=Clo, edge_lo_status = openlo ? "OPEN" : "CLOSED",
     theta=θ, Cmax=Cmax, clay_ceiling=clayceil, sigma_bulk_min=Clo/L,
     lead=lead, floor_tds=floor_tds, floor_tds_at_mhi = openlo ? 0.0 : σw_hi*MGL,
     Xlim=Xlim, phi_max=φmax, m_used=mlo, m_other=mhi,
     verdict = (!openlo && σw_lo*MGL > Xlim) ? "FAIL under stated ceilings" :
               (openlo ? "no floor (lower edge OPEN)" : "floor below the line"))
end

function main(tags)
    rows = [floorfor(t) for t in tags]
    c = SiteConfig.resolve(joinpath(ROOT,"configs","synth_$(tags[1]).toml"))["ceilings"]
    @printf("STAGE 4 SALINITY FLOOR under stated ceilings\n")
    @printf("  sigma_clay_max %.2f S/m   C_max = %.2f (1 - theta/%.2f)   porosity_max %.2f   m = %s bound\n\n",
            c["sigma_clay_max"], c["clay_c0"], c["clay_phi_ref"], c["porosity_max"], c["m_for_floor"])
    @printf("  %-5s %9s %7s %8s %9s %11s %10s  %s\n", "col", "C_lo (S)", "edge",
            "theta", "clay S/m", "C_lo/L S/m", "floor mg/L", "verdict")
    println("  " * "-"^88)
    for r in rows
        @printf("  %-5s %9.5f %7s %8.4f %9.5f %11.5f %10.0f  %s\n",
                r.tag, r.C_lo, r.edge_lo_status, r.theta, r.clay_ceiling,
                r.sigma_bulk_min, r.floor_tds, r.verdict)
    end
    println("  " * "-"^88)
    # WHY the answer is zero, which is not only the open edges. Compared against
    # the conductivity the survey actually sees (Stage 5's sigma_bulk p50 -- no
    # truth is read here), the clay ceiling alone exceeds the whole measurement.
    # sigma_bulk_floor is the conductivity below which no floor can reach the
    # salinity line even with zero clay: X_lim/MGL * porosity_max^m.
    println("  WHY ZERO -- the clay ceiling against what the survey measures")
    for r in rows
        S5 = load(joinpath(OUT,"synth_$(r.tag)","stage5_synth_$(r.tag).jld2"))
        sb = S5["sigma_bulk_p50"]
        sbfloor = r.Xlim/MGL * r.phi_max^r.m_used
        breakeven = sb / r.Cmax
        @printf("    %-5s AEM sigma_bulk p50 %.5f | clay ceiling %.5f = %5.1fx it | sigma_clay_max would need < %.4f S/m\n",
                r.tag, sb, r.clay_ceiling, r.clay_ceiling/sb, breakeven)
        @printf("    %-5s a %.0f mg/L floor needs sigma_bulk >= %.5f S/m at phi=%.2f, m=%.2f -- %s\n",
                "", r.Xlim, sbfloor, r.phi_max, r.m_used,
                sb >= sbfloor ? "this column clears that" : "this column CANNOT, at any clay ceiling")
    end
    mkpath(joinpath(ROOT,"outputs","stage4_floor"))
    jldopen(joinpath(ROOT,"outputs","stage4_floor","floor.jld2"),"w") do f
        for r in rows; f[r.tag] = r; end
    end
    open(joinpath(ROOT,"outputs","stage4_floor","floor.csv"),"w") do f
        println(f, "column,C_lo_S,edge_lo_status,theta_snmr,clay_ceiling_S_per_m,floor_tds_mg_per_L,verdict")
        for r in rows
            @printf(f, "\"%s\",%.6g,\"%s\",%.6g,%.6g,%.6g,\"%s\"\n",
                    r.tag, r.C_lo, r.edge_lo_status, r.theta, r.clay_ceiling,
                    r.floor_tds, r.verdict)
        end
    end
    println("\nwrote outputs/stage4_floor/floor.{jld2,csv}")
end
main(isempty(ARGS) ? ["3P","2P","3S","2S"] : ARGS)
