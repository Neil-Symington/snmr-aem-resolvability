# check_open_edge.jl -- is the band's upper edge set by the DATA or by sigma_min?
#
#   julia --project=. check_open_edge.jl
#
# The tier-P gates run with the solver clamped to sigma in [1e-3, 10] S/m, and on
# every column chi2 never reaches N + 3.84 before the nodes stop being
# admissible. Two explanations: the data genuinely do not bound the window's
# conductance above, or the sigma_min FLOOR does -- raising the window's
# conductance forces the surrounding cells down to compensate, and once they sit
# on 1e-3 the constraint cannot be met at all.
#
# The package-budget check does not separate them: sigma_min over the whole grid
# is 0.27 S against a Stage 1 package conductance of 26-141 S, so the budget is
# nowhere near binding. The compensation must be local, which the budget
# arithmetic cannot see. Moving the floor is the only test that can.
#
# So: re-run at lambda2 = 1 with sigma_min = 1e-4, everything else fixed, walk
# the upper edge up from the known admissible one and bisect any crossing.
using JLD2, Statistics, Printf, LinearAlgebra, Random, TOML, SparseArrays
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "config.jl")); using .SiteConfig
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "occam_core.jl")); using .OccamCore
const SK = HiQGA.transD_GP.SkyTEM1DInversion
const MGL = 6400.0
cellsof(z)=(mid=(z[1:end-1].+z[2:end])./2; (top=vcat(0.0,mid), bot=vcat(mid,z[end]+(z[end]-z[end-1])/2)))
function winw(z,a,b); g=cellsof(z); w=[max(0.0,min(g.bot[i],b)-max(g.top[i],a)) for i in eachindex(z)]
    ix=findall(>(1e-12),w); (ix=ix,w=w[ix]); end

function operator(A)
    ρ0 = fill(100.0, length(A["z"])); ρ0[A["z"] .< 0] .= 1e12
    aem = SK.dBzdt(; timeslow=A["times_low"], ramplow=A["ramp_low"], zRxlow=A["z_rx"], zTxlow=A["z_tx"],
        timeshigh=A["times_high"], ramphigh=A["ramp_high"], zRxhigh=A["z_rx"], zTxhigh=A["z_tx"],
        rRx=A["r_rx"], rTx=A["r_tx"], z=copy(A["z"]), ρ=ρ0, lowpassfcs=A["lowpass_fc"], calcjacobian=true)
    aem.dlow=copy(A["d_low"]); aem.dhigh=copy(A["d_high"])
    aem.σlow=copy(A["sigma_low"]); aem.σhigh=copy(A["sigma_high"])
    aem.selectlow=copy(A["select_low"]); aem.selecthigh=copy(A["select_high"])
    aem.ndatalow=Int(A["ndata_low"]); aem.ndatahigh=Int(A["ndata_high"])
    sg=vcat(aem.σlow[aem.selectlow], aem.σhigh[aem.selecthigh])
    aem.W=sparse(Diagonal(1 ./ sg)); aem.res=zeros(length(sg))
    aem.J=zeros(length(sg), length(aem.ρ)-1); aem
end

function main()
    setbounds!(-4.0, 1.0)                       # sigma_min = 1e-4 S/m
    @printf("solver bounds now %s  (sigma %.0e to %.0f S/m)\n\n", string(getbounds()), 1e-4, 10)
    for tag in ("3P","2P","3S","2S")
        cfg = SiteConfig.resolve(joinpath(ROOT,"configs","synth_$tag.toml"))
        info = SiteConfig.assertvalid(cfg; root=ROOT)
        odir = joinpath(ROOT, cfg["output"]["dir"])
        A  = load(joinpath(odir,"stage0_synth_$tag.jld2"))
        G  = load(joinpath(odir,"stage4_synth_$tag.jld2"))
        P4 = load(joinpath(odir,"stage4post_synth_$tag.jld2"))
        E  = load(joinpath(HERE, "..", "data","cond_C2_typical$(tag)T.jld2"))
        aem = operator(A); zall=A["zall"]; n=length(aem.ρ)-1; ndata=A["ndata"]
        dzc = diff(vcat(0.0, collect(Float64,zall)))
        R = makeR1(n); m0 = fill(-2.0, n)
        w = (G["window"].top, G["window"].base); thick = w[2]-w[1]
        celltop = vcat(0.0, collect(Float64,zall)[1:end-1]); cellbot = collect(Float64,zall)
        ov = [max(0.0, min(cellbot[i],w[2]) - max(celltop[i],w[1])) for i in 1:n]
        fixidx = findall(>(1e-9), ov); dzeff = copy(dzc); dzeff[fixidx] .= ov[fixidx]
        mref, _ = occam(aem, m0, m0, R, ndata; mode=:minimise)
        target = ndata + 3.84
        probe(C) = begin
            m, chi, _H, nf = occam(aem, mref, m0, R, ndata; fixidx=fixidx, fixval=C,
                mode=:map, ctype=:cond, dz=dzeff, lambda2=1.0, beta2=0.1, withH=false, info=true)
            (chi, nf.nbound==0 && nf.conres<1e-8 && nf.settled, nf.nbound)
        end
        # walk up from the old upper edge until chi2 crosses or the node fails
        lt = log10(P4["band_C"].hi / thick)
        @printf("=== %s  window %.2f-%.2f m, old upper edge sigma_eff %.5g S/m (C %.4g S)\n",
                tag, w..., 10^lt, P4["band_C"].hi)
        prev = (lt, probe(10^lt * thick))
        @printf("    %-10s %-11s %-6s %s\n","sigma_eff","chi2","bound","status")
        @printf("    %-10.5g %-11.3f %-6d %s\n", 10^lt, prev[2][1], prev[2][3],
                prev[2][2] ? "admissible" : "INADMISSIBLE")
        crossed = false; cur = lt
        for _ in 1:14
            cur += 0.2; cur > 1.0 && break
            chi, ok, nb = probe(10^cur * thick)
            @printf("    %-10.5g %-11.3f %-6d %s\n", 10^cur, chi, nb,
                    !ok ? "INADMISSIBLE" : chi > target ? "chi2 > target  <-- CROSSED" : "admissible")
            if ok && chi > target; crossed = true; break; end
            if !ok; break; end
            prev = (cur, (chi, ok, nb))
        end
        if crossed
            a, b = prev[1], cur
            for _ in 1:12
                mid=(a+b)/2; chi, ok, _ = probe(10^mid*thick)
                (ok && chi <= target) ? (a=mid) : (b=mid)
                abs(b-a) < 5e-3 && break
            end
            edge = (a+b)/2
            ww = winw(E["zsmr"], w...); Θ=Matrix{Float64}(reduce(hcat,E["ens"]))
            Random.seed!(7); NS=20000
            ms = info.m_prior[1] .+ (info.m_prior[2]-info.m_prior[1]).*rand(NS)
            js = rand(1:size(Θ,2), NS)
            K = [sum((@view Θ[ww.ix,js[s]]) .^ ms[s] .* ww.w) for s in 1:NS]
            bnd = quantile((10^edge*thick) ./ K .* MGL, 0.95)
            @printf("    NEW UPPER EDGE sigma_eff %.5g S/m (C %.4g S); chi2 reaches N+3.84 -> CLOSED\n",
                    10^edge, 10^edge*thick)
            @printf("    bound at that edge: TDS p95 = %.0f mg/L   (was %.0f at the 1e-3 floor)\n\n",
                    bnd, P4["tds_data_only_p95"])
        else
            @printf("    NO CROSSING up to sigma_eff %.3g S/m: chi2 never reaches N+3.84 -> still OPEN\n\n",
                    10^cur)
        end
        flush(stdout)
    end
end
main()
