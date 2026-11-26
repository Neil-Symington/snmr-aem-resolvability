# report_html.jl -- the interactive viewer for the field run: one page per site
# (AEM ensemble, SNMR ensemble and data fit, the profile's gate models, the
# misfit profile and the rock-physics sliders) plus an index of every site.
#
#   julia --project=. fieldrun/stages/report_html.jl --config <toml> [--site S] [--index] [--zmax 100]
#
#   no flags   every site in [run] sites, then the index
#   --site S   that site's page only (one process per site in a batch)
#   --index    the index only, from whatever the stages have written
#
# Writes <outputs>/viewer/<site>.html (self-contained, opened from disk), the
# same payload as <site>.json, and index.html. Runs no inversion.
#
# THE PAGE IS NEIL'S. viewer/page.html and viewer/script.js are tools/make_report.jl's
# PAGE and SCRIPT (branch m1-production-rerun, 24 Nov 2026), vendored so the field
# workflow does not depend on the paper branch. The changes are all guarded by
# D.field or by a missing value, and the synthetic page is unchanged by them:
#   - no truth at a field site: the white line, the "true" column and key entry,
#     the clay-truth row and the truth conductance on the misfit plot go;
#   - the profile's nodes are the window's GEOMETRIC-MEAN CONDUCTIVITY (S/m),
#     not a conductance (S), so the labels say so and the thickness is 1 --
#     the page's sigma_bulk = C / thick is then the node value itself;
#   - ticks to 1e-3 S/m on the misfit plot; a summary line under the title;
#   - a site with no detected interval has no profile or gate, and the page
#     draws the two ensembles and the data fit only.
# A diff against make_report.jl's blocks shows exactly that.
#
# WHAT EACH PANEL IS, at a field site.
#   AEM       aem_td_<site>.jld2 sigma_draws: the trans-D draws the stage kept
#             (a thinned subset of the T = 1 chain), on the AEM grid. zbnd is
#             the layer tops, so the cells are [zbnd[i], zbnd[i+1]].
#   SNMR      smr_transd_<site>.jld2 theta: every kept draw pooled over kernel
#             draws and replicate ladders, on the kernel nodes zbnd -- the
#             grid detect_aquifer.jl uses, with cell edges midway between.
#   data fit  detect_<site>.jld2: V0 at each pulse moment, the noise
#             sigma_eff = sqrt(s p50) sigma-hat, and the p05/p50/p95 of the
#             forward-modelled draws detect computed.
#   gate      aem_profile_<site>.jld2: the min-chi2 model at each forced value,
#             admissible where chi2 <= rlevel (ndata + 3.84). The gate's
#             interval, theta and m ranges come from gate_<site>.jld2.
include(joinpath(@__DIR__, "..", "common", "stage_prelude.jl"))
using JLD2, Printf, Statistics, Dates

const VIEW = joinpath(OUT, "viewer")
const INDEX_ONLY = any(==("--index"), ARGS)
const ZMAX = let i = findfirst(==("--zmax"), ARGS)
    i === nothing ? 100.0 : parse(Float64, ARGS[i+1])
end
const PAGE   = read(joinpath(@__DIR__, "viewer", "page.html"), String)
const SCRIPT = read(joinpath(@__DIR__, "viewer", "script.js"), String)
const COMMIT = try
    strip(read(`git -C $(joinpath(FIELDRUN, "..")) rev-parse --short HEAD`, String))
catch
    "unknown"
end

# ------------------------------------------------------------------ JSON
# Hand-written, as make_report.jl's is: there is no JSON package in the
# project. "</" is escaped so no string can close the <script> it sits in.
jesc(s) = replace(string(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n",
                  "\r" => "", "\t" => " ", "</" => "<\\/")
js(::Nothing) = "null"
js(x::Bool) = x ? "true" : "false"
js(x::Integer) = string(x)
js(x::Real) = isfinite(x) ? @sprintf("%.6g", x) : "null"
js(s::AbstractString) = "\"" * jesc(s) * "\""
js(s::Symbol) = js(string(s))
js(v::AbstractVector) = "[" * join((js(x) for x in v), ",") * "]"
js(t::NamedTuple) = "{" * join((js(string(k)) * ":" * js(v) for (k, v) in pairs(t)), ",") * "}"
js(d::AbstractDict) = "{" * join((js(string(k)) * ":" * js(d[k]) for k in sort(collect(keys(d)))), ",") * "}"

function binrow(vals, lo, hi, nb)
    c = zeros(Int, nb); w = (hi - lo) / nb
    @inbounds for x in vals
        isfinite(x) || continue
        c[clamp(floor(Int, (x - lo) / w) + 1, 1, nb)] += 1
    end
    c
end
rowq(M, p) = [quantile(view(M, i, :), p) for i in axes(M, 1)]
h(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
fmt(x, d = 2) = (x isa Real && isfinite(x)) ? @sprintf("%.*f", d, x) : "–"
fmtg(x) = (x isa Real && isfinite(x)) ? @sprintf("%.3g", x) : "–"

const AEM_LO, AEM_HI, AEM_NB = -4.0, 0.6, 92        # log10 S/m; the trans-D prior is -4 to 0.5
const SNMR_LO, SNMR_HI, SNMR_NB = 0.0, 0.45, 60     # v/v; [mcmc] fbounds tops out at 0.45

ofile(kind, fid) = joinpath(OUT, "$(kind)_$(fid).jld2")
lazy(f, k, d = nothing) = jldopen(f, "r") do F; haskey(F, k) ? F[k] : d; end

"Everything the index needs about one site, read lazily from the stage outputs."
function site_record(fid)
    r = Dict{String,Any}("site" => fid)
    have = Dict(k => isfile(ofile(k, fid)) for k in
                ("aem_occam", "aem_td", "smr_transd", "detect", "aem_profile", "gate"))
    r["have"] = have
    if have["aem_occam"]
        f = ofile("aem_occam", fid)
        r["line"] = lazy(f, "line", missing); r["sep_m"] = lazy(f, "sep_m", NaN)
        r["occam_chi2n"] = lazy(f, "chi2n", NaN)
        r["noise_cal"] = lazy(f, "noise_calibration", "inflate")
        r["noise_x"] = r["noise_cal"] == "two_term" ? lazy(f, "noise_conv_factor", NaN) : lazy(f, "noise_scale", NaN)
    end
    if have["aem_td"]
        r["ml_sigma"] = lazy(ofile("aem_td", fid), "ml_sigma_factor_p50", NaN)
    end
    if have["smr_transd"]
        f = ofile("smr_transd", fid)
        ls = lazy(f, "log_s", Float64[])
        r["s_p50"] = isempty(ls) ? NaN : exp(median(ls))
        r["s_treatment"] = lazy(f, "s_treatment", "sampled")
        r["gate_basis"] = lazy(f, "gate_basis", "")
        r["rhat"] = lazy(f, "rhat_replicate_max", lazy(f, "rhat_max", NaN))
        r["ess"] = lazy(f, "ess_replicate", lazy(f, "ess_bulk_min", NaN))
        r["conv_pass"] = lazy(f, "gates_pass", false)
        r["nens"] = lazy(f, "nensembles", 0); r["nrep"] = lazy(f, "nreplicates", 1)
    end
    if have["detect"]
        f = ofile("detect", fid)
        W = lazy(f, "W", Float64[])
        r["verdict"] = lazy(f, "verdict", ""); r["win"] = lazy(f, "win", [NaN, NaN])
        r["win_source"] = lazy(f, "win_source", ""); r["zlim"] = lazy(f, "zlim", NaN)
        r["W50"] = isempty(W) ? NaN : median(W)
        r["W05"] = isempty(W) ? NaN : quantile(W, 0.05)
        r["W95"] = isempty(W) ? NaN : quantile(W, 0.95)
        r["PW2"] = isempty(W) ? NaN : mean(W .>= 2)
        r["detectability"] = lazy(f, "detect_rule", "not recorded")
    end
    if have["aem_profile"]
        f = ofile("aem_profile", fid)
        r["sig_lo"] = lazy(f, "edge_lo", NaN); r["sig_hi"] = lazy(f, "edge_hi", NaN)
        r["lo_open"] = lazy(f, "edge_lo_open", false); r["hi_open"] = lazy(f, "edge_hi_open", false)
    end
    if have["gate"]
        f = ofile("gate", fid)
        r["gate_verdict"] = lazy(f, "verdict", "")
        r["Xc95_mgl"] = lazy(f, "Xcert_p95", NaN) * lazy(f, "mgl", NaN)
        pf = lazy(f, "pfresh", Dict{String,Float64}())
        r["pfresh"] = pf
    end
    r
end

"A gate prior tag without its numbers: \"clay prior 0.10-2.0 S/m\" -> \"clay prior\"."
shorttag(k) = strip(replace(String(k), r"\s*[\d(].*$" => ""))

"The one-line HTML summary under a site page's title."
function summary_html(r)
    b = String[]
    push!(b, "<a href=\"index.html\">&larr; all sites</a>")
    if haskey(r, "verdict")
        w = r["win"]
        push!(b, "<b>$(h(r["verdict"]))</b>")
        isfinite(w[1]) && push!(b, @sprintf("interval %.1f–%.1f m (%s)", w[1], w[2], h(r["win_source"])))
        push!(b, "W p50 $(fmt(r["W50"])) m (p05–p95 $(fmt(r["W05"]))–$(fmt(r["W95"])))")
        push!(b, "depth limit $(fmt(r["zlim"], 0)) m <em>($(h(r["detectability"])) rule)</em>")
    end
    if haskey(r, "s_p50")
        cls = r["conv_pass"] ? "ok" : "bad"
        push!(b, "SNMR s p50 $(fmt(r["s_p50"])) ($(h(r["s_treatment"])))")
        push!(b, "<span class=\"$cls\">R&#770; $(fmt(r["rhat"])) ESS $(fmt(r["ess"], 0)) $(r["conv_pass"] ? "converged" : "not converged")</span>")
    end
    haskey(r, "occam_chi2n") &&
        push!(b, "AEM Occam &chi;&sup2;/N $(fmt(r["occam_chi2n"], 3)), noise $(h(r["noise_cal"])) &times;$(fmt(r["noise_x"], 3))")
    if haskey(r, "sig_hi")
        lo = (r["lo_open"] ? "&le;" : "") * fmtg(r["sig_lo"])
        hi = fmtg(r["sig_hi"]) * (r["hi_open"] ? "+" : "")
        push!(b, "AEM admits $(lo)–$(hi) S/m")
    end
    if haskey(r, "gate_verdict")
        push!(b, "<b>$(h(r["gate_verdict"]))</b>, X<sub>cert</sub> p95 $(fmt(r["Xc95_mgl"], 0)) mg/L")
        for (k, v) in sort(collect(r["pfresh"]))
            push!(b, "P(fresh | $(h(shorttag(k)))) $(fmt(v, 2))")
        end
    end
    join(("<span>" * x * "</span>" for x in b), "")
end

"The page's data object for one site."
function site_payload(fid)
    r = site_record(fid)
    hv = r["have"]
    (hv["aem_td"] && hv["smr_transd"]) ||
        error("$fid: needs aem_td and smr_transd outputs (have: $(join([k for (k, v) in hv if v], ", ")))")

    # ---- AEM ensemble
    A = load(ofile("aem_td", fid))
    σ = Float64.(A["sigma_draws"]); zall = Float64.(A["zall"]); zb = Float64.(A["zbnd"])
    top = zb; bot = vcat(zb[2:end], 2zall[end] - zb[end])
    ah = [binrow(log10.(view(σ, i, :)), AEM_LO, AEM_HI, AEM_NB) for i in axes(σ, 1)]

    # ---- SNMR ensemble
    S = load(ofile("smr_transd", fid))
    θ = Float64.(S["theta"]); zs = Float64.(S["zbnd"])
    sh = [binrow(view(θ, i, :), SNMR_LO, SNMR_HI, SNMR_NB) for i in axes(θ, 1)]

    # ---- detect: interval, depth limit, the data fit
    win = [NaN, NaN]; dlim = nothing; fit = nothing
    if hv["detect"]
        Dt = load(ofile("detect", fid))
        win = Float64.(Dt["win"]); dlim = Float64(Dt["zlim"])
        if haskey(Dt, "pred_V0")
            P = Dt["pred_V0"]
            fit = (q = Float64.(Dt["q"]), V = Float64.(Dt["V0"]), sig = Float64.(Dt["sigma_eff"]),
                   p05 = P[:, 1], p50 = P[:, 2], p95 = P[:, 3],
                   chi2n = haskey(Dt, "chi2n_draws") ? median(Dt["chi2n_draws"]) : NaN,
                   tier = "", src = "detect_$(fid).jld2", haspred = true)
        end
    end

    # ---- profile + gate
    gate = nothing; pet = nothing; tabs = NaN; trel = NaN; chiref = NaN; mb = [AEM_LO, AEM_HI]; ndata = 0
    if hv["aem_profile"] && hv["gate"]
        P = load(ofile("aem_profile", fid)); G = load(ofile("gate", fid))
        C = Float64.(P["targets"]); χ = Float64.(P["chi2_prof"]); M = P["models"]
        mb = [Float64(P["bound_lo"]), Float64(P["bound_hi"])]
        tabs = Float64(P["rlevel"]); chiref = Float64(P["chi2_ref"]); trel = chiref + 3.84
        ndata = Int(P["ndata"])
        ok = [all(isfinite, M[:, k]) && isfinite(χ[k]) for k in eachindex(C)]
        ix = findall(ok)
        adm = [χ[k] <= tabs for k in ix]
        if any(adm)
            ob = [findall(x -> x <= mb[1] + 1e-9 || x >= mb[2] - 1e-9, M[:, k]) .- 1 for k in ix]
            best = C[ix][argmin(χ[ix])]
            gate = (C = C[ix], chi2 = χ[ix], adm = Int.(adm), best = best,
                    m = [Float64.(M[:, k]) for k in ix], ob = ob,
                    ms = nothing, mschi2 = nothing, fb = nothing, msn = 0,
                    src = "aem_profile_$(fid).jld2")
            th = Float64.(G["theta_interval"]); md = Float64.(G["m_draws"])
            pet = (mgl = Float64(G["mgl"]), phi = [quantile(th, 0.05), quantile(th, 0.95)],
                   m = [minimum(md), maximum(md)], scmax = nothing, c0 = nothing,
                   phiref = nothing, cp = [0.0, 1e9], src = "gate_$(fid).jld2")
        else
            @warn "$fid: no profile node is admissible; the page shows the ensembles only"
        end
    end

    ndraw_aem = size(σ, 2)
    prov = @sprintf("AEM trans-D %s draws kept of %s &middot; SNMR %d kernel draws &times; %d ladders, %s draws, s %s &middot; ",
                    string(ndraw_aem), string(get(A, "ndraws_total", ndraw_aem)),
                    get(S, "nensembles", 0), get(S, "nreplicates", 1), string(size(θ, 2)),
                    h(get(S, "s_treatment", "sampled"))) *
           (gate === nothing ? "no gate &middot; " :
            @sprintf("profile %d nodes, %.0e–%.1f S/m &middot; ", length(gate.C), 10^mb[1], 10^mb[2])) *
           "run $(h(FCFG["run"]["name"])) &middot; commit $(h(COMMIT)) &middot; " *
           Dates.format(now(), "yyyy-mm-dd HH:MM")

    (field = true, gen = Dates.format(now(), "yyyy-mm-dd HH:MM"),
     win = round.(win, digits = 2), zmax = ZMAX, zb = isfinite(win[2]) ? min(ZMAX, something(dlim, ZMAX)) : nothing,
     thick = 1.0, ndata = ndata, tabs = tabs, trel = trel, chiref = chiref,
     chibest = gate === nothing ? NaN : minimum(gate.chi2), chirefstore = chiref, mb = mb,
     aem = (z = zall, top = top, bot = bot, lo = AEM_LO, hi = AEM_HI, nb = AEM_NB, h = ah,
            rmax = [max(1, maximum(x)) for x in ah],
            p05 = rowq(σ, 0.05), p50 = rowq(σ, 0.5), p95 = rowq(σ, 0.95), truth = Float64[],
            ndraw = ndraw_aem, nens = 1, src = "aem_td_$(fid).jld2"),
     snmr = (z = zs, lo = SNMR_LO, hi = SNMR_HI, nb = SNMR_NB, h = sh,
             rmax = [max(1, maximum(x)) for x in sh],
             p05 = rowq(θ, 0.05), p50 = rowq(θ, 0.5), p95 = rowq(θ, 0.95), truth = Float64[],
             ndraw = size(θ, 2), nens = get(S, "nensembles", 0), src = "smr_transd_$(fid).jld2",
             note = r["conv_pass"] ? "converged" : @sprintf("not converged: R-hat %.2f", r["rhat"]),
             dlim = dlim),
     gate = gate, truthC = nothing, truthC_fine = nothing, snmr_run = basename(OUT),
     colname = fid, aem_dirs = "", sens = nothing, pet = pet, claytrue = nothing, fit = fit,
     prov = prov, summ = summary_html(r), csvname = "viewer_$(fid).csv")
end

function write_page(fid)
    mkpath(VIEW)
    D = site_payload(fid)
    J = js(D)
    jf = joinpath(VIEW, "$fid.json")
    write(jf, J)
    # validate before the page is written, as make_report.jl does
    ok = try strip(read(`python3 -c "import json,sys; json.load(open(sys.argv[1])); print('ok')" $jf`, String)) == "ok"
         catch; nothing end
    ok === false && error("$fid: payload is not valid JSON ($jf)")
    i = findfirst("</style>", PAGE)
    head, body = PAGE[1:last(i)], PAGE[last(i)+1:end]
    title = "$fid — AEM, SNMR and gate"
    head = replace(head, "@TITLE@" => h(title)); body = replace(body, "@H1@" => h(title))
    f = joinpath(VIEW, "$fid.html")
    open(f, "w") do io
        print(io, "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
              "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
              head, "\n</head>\n<body>\n", body,
              "\n<script>\nconst D=", J, ";\n</script>\n", SCRIPT, "\n</body>\n</html>\n")
    end
    @printf("  %-14s %s  (%.0f kB)%s\n", fid, relpath(f, OUT), filesize(f) / 1024,
            ok === nothing ? "  [JSON not checked: no python3]" : "")
end

function write_index(sites)
    mkpath(VIEW)
    rows = String[]
    for s in sites
        r = try site_record(s) catch e
            @warn "index: $s unreadable" exception = e; Dict{String,Any}("site" => s, "have" => Dict())
        end
        hv = get(r, "have", Dict())
        linked = isfile(joinpath(VIEW, "$s.html"))
        name = linked ? "<a href=\"$(h(s)).html\">$(h(s))</a>" : h(s)
        missing_st = [k for k in ("aem_occam", "aem_td", "smr_transd", "detect", "aem_profile", "gate")
                      if !get(hv, k, false)]
        w = get(r, "win", [NaN, NaN])
        conv = haskey(r, "rhat") ?
            "<span class=\"$(r["conv_pass"] ? "ok" : "bad")\">$(fmt(r["rhat"])) / $(fmt(r["ess"], 0))</span>" : "–"
        pf = get(r, "pfresh", Dict())
        pfs = isempty(pf) ? "–" : join(["$(fmt(v, 2))" for (k, v) in sort(collect(pf))], " / ")
        push!(rows, string("<tr><td class=\"s\">", name, "</td>",
            "<td>", h(string(get(r, "line", "–"))), "</td>",
            "<td>", fmt(get(r, "sep_m", NaN), 0), "</td>",
            "<td>", fmt(get(r, "occam_chi2n", NaN), 3), "</td>",
            "<td>", fmt(get(r, "s_p50", NaN)), "</td>",
            "<td>", conv, "</td>",
            "<td>", isfinite(w[1]) ? @sprintf("%.1f–%.1f", w[1], w[2]) : "–", "</td>",
            "<td>", fmt(get(r, "W50", NaN)), "</td>",
            "<td>", fmt(get(r, "PW2", NaN), 3), "</td>",
            "<td>", fmt(get(r, "zlim", NaN), 0), "</td>",
            "<td class=\"v\">", h(get(r, "verdict", "–")), "</td>",
            "<td>", fmtg(get(r, "sig_hi", NaN)), get(r, "hi_open", false) ? "+" : "", "</td>",
            "<td>", fmt(get(r, "Xc95_mgl", NaN), 0), "</td>",
            "<td>", pfs, "</td>",
            "<td class=\"v\">", h(get(r, "gate_verdict", "–")), "</td>",
            "<td class=\"m\">", isempty(missing_st) ? "" : "missing: " * join(missing_st, ", "), "</td></tr>"))
    end
    pfkeys = let k = String[]
        for s in sites
            f = ofile("gate", s); isfile(f) || continue
            k = sort(collect(keys(lazy(f, "pfresh", Dict())))); break
        end
        k
    end
    css = """
:root{--bg:#f1f4f4;--panel:#fff;--ink:#16262b;--soft:#4a5f65;--mute:#7b9095;--rule:#ccd8d8;
  --acc:#0d6a68;--hot:#b4502f;color-scheme:light}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#0d1517;--panel:#162023;
  --ink:#e3ebea;--soft:#aec0c0;--mute:#8098a0;--rule:#2a393d;--acc:#51b5ad;--hot:#dd7d60;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#0d1517;--panel:#162023;--ink:#e3ebea;--soft:#aec0c0;--mute:#8098a0;
  --rule:#2a393d;--acc:#51b5ad;--hot:#dd7d60;color-scheme:dark}
body{background:var(--bg);color:var(--ink);margin:0;
  font:13px/1.45 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.app{padding:14px 16px 22px;max-width:1600px;margin:0 auto}
h1{font-size:16px;margin:0 0 2px}
.prov{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:10.5px;color:var(--mute);margin-bottom:12px}
.tw{overflow-x:auto;background:var(--panel);border:1px solid var(--rule);border-radius:4px}
table{border-collapse:collapse;width:100%;font-size:12px}
th{font-size:10px;font-weight:600;letter-spacing:.05em;text-transform:uppercase;color:var(--mute);
  text-align:right;padding:7px 8px;border-bottom:1px solid var(--rule);white-space:nowrap;vertical-align:bottom}
td{padding:5px 8px;border-bottom:1px solid var(--rule);text-align:right;white-space:nowrap;
  font-variant-numeric:tabular-nums;font-family:ui-monospace,Menlo,Consolas,monospace;font-size:11.5px}
td.s,th.s{text-align:left;font-family:inherit;font-weight:600}
td.v{text-align:left;font-family:inherit;white-space:normal;min-width:180px}
td.m{text-align:left;color:var(--hot);font-family:inherit;font-size:11px}
a{color:var(--acc)}
.ok{color:var(--acc);font-weight:600}.bad{color:var(--hot);font-weight:600}
.note{font-size:11.5px;color:var(--soft);margin-top:10px;max-width:900px}
"""
    f = joinpath(VIEW, "index.html")
    open(f, "w") do io
        print(io, "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
              "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
              "<title>$(h(FCFG["run"]["name"])) sites</title>\n<style>", css, "</style>\n</head>\n<body>\n<div class=\"app\">\n",
              "<h1>$(h(FCFG["run"]["name"])) — $(length(sites)) sites</h1>\n",
              "<div class=\"prov\">config $(h(basename(CFGPATH))) &middot; commit $(h(COMMIT)) &middot; ",
              Dates.format(now(), "yyyy-mm-dd HH:MM"), "</div>\n<div class=\"tw\"><table>\n<thead><tr>",
              "<th class=\"s\">site</th><th>AEM line</th><th>sep m</th><th>Occam &chi;&sup2;/N</th>",
              "<th>SNMR s p50</th><th>R&#770; / ESS</th><th>interval m</th><th>W p50 m</th><th>P(W&ge;2)</th>",
              "<th>depth limit m</th><th class=\"s\">detection</th><th>&sigma;<sub>hi</sub> S/m</th>",
              "<th>X<sub>cert</sub> p95 mg/L</th><th>P(fresh)</th>",
              "<th class=\"s\">gate</th><th class=\"s\"></th></tr></thead>\n<tbody>\n",
              join(rows, "\n"), "\n</tbody></table></div>\n",
              "<div class=\"note\">R&#770; / ESS is the SNMR convergence gate as the stage recorded it (red: not converged; ",
              "the products may still be stable across ladders). &sigma;<sub>hi</sub> is the upper edge of the conductivity ",
              "the AEM data admit over the interval (+ = open). X<sub>cert</sub> p95 is the clay-free certification. ",
              "P(fresh) is under each gate prior, in order: ", h(join(shorttag.(pfkeys), " / ")), ".</div>\n",
              "</div>\n</body>\n</html>\n")
    end
    println("  index: ", relpath(f, OUT))
end

function main()
    sites = RUN_SITES
    if !INDEX_ONLY
        nok = 0
        for s in sites
            try write_page(s); nok += 1 catch e
                @warn "viewer: $s failed" exception = (e, catch_backtrace())
            end
        end
        nok == 0 && (println("FAILED: no site page written"); exit(1))
    end
    # the index lists every site in the run, not just this process's --site
    (isempty(ONLY) || INDEX_ONLY) && write_index(FCFG["run"]["sites"])
end

main()
