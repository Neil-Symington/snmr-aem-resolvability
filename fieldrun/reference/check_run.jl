# check_run.jl -- compare a finished batch run against stored reference values.
#
#   julia --project=. fieldrun/reference/check_run.jl --config <control file>
#
# WHY THIS EXISTS. A batch run used to end by printing "complete" and a table,
# and the only way to judge it was to look at the table and nod. The first full
# two-site run printed "snmr_transd ok" for both sites after both had thrown,
# and all three sweep stages then printed an empty table and exited 0. A run
# has to be checkable by something that cannot be satisfied by plausible-looking
# output, so this script exits non-zero when a run disagrees with the reference.
#
# WHAT IS CHECKED, and what deliberately is not:
#
#   [deterministic]  The Occam inversion has no Monte Carlo in it. chi2/N, the
#                    data count and sigma(z) at fixed depths reproduce to the
#                    printed digits on any machine, so they are compared with a
#                    tight relative tolerance and a disagreement is a defect.
#
#   [structural]     Everything downstream of trans-D is a posterior. In smoke
#                    mode it is 2001 samples and NOT CONVERGED -- the stage says
#                    so itself ("gates FAIL") -- so comparing its quantiles to
#                    a stored number would be checking noise. What is checked is
#                    what must hold for ANY valid run: the files exist, the
#                    conductivity band is ordered and positive, the target
#                    interval lies inside the grid, and every verdict is one of
#                    the strings the stages can actually emit.
#
#   [full]           Converged reference values. Absent until a full-mode run
#                    has been done; the checker says so rather than inventing
#                    them, and smoke numbers are NOT promoted into this section.
#
# So: a passing smoke check means the workflow is wired correctly end to end.
# It does not mean the posterior is quotable, and nothing here should be read
# as saying it is.

const COMMON = normpath(joinpath(@__DIR__, "..", "common"))
include(joinpath(COMMON, "fieldconfig.jl"))
using .FieldConfig
using JLD2, TOML, Printf, Statistics

const CFGPATH = let i = findfirst(==("--config"), ARGS)
    i === nothing && error("check_run.jl needs --config <control file>")
    ARGS[i+1]
end
const FCFG = FieldConfig.load_config(CFGPATH)
const OUT  = FieldConfig.cfgpath(FCFG, "outputs")
const MODE = String(get(FCFG["run"], "mode", "full"))
const SITES = String.(FCFG["run"]["sites"])
const STAGES = String.(FCFG["run"]["stages"])

const REFPATH = let i = findfirst(==("--reference"), ARGS)
    i === nothing ? joinpath(@__DIR__, "udf_example.toml") : ARGS[i+1]
end
isfile(REFPATH) || error("no reference file at $REFPATH")
const REF = TOML.parsefile(REFPATH)

# Verdicts the stages can emit. The gate's are exact strings; detection's are
# built as a base clause plus " above NN m" plus an optional pinning clause, so
# only the base can be matched and it is matched as a PREFIX. A verdict that
# matches nothing means the stage's vocabulary changed without the reference,
# which is worth failing on -- the verdict is the sentence that gets quoted.
const DETECT_PREFIXES = String.(get(REF, "detect_verdict_prefixes", String[]))
const GATE_VERDICTS   = Set(String.(get(REF, "gate_verdicts", String[])))

const FAILS = String[]
fail(msg) = (push!(FAILS, msg); @printf("  FAIL  %s\n", msg))
pass(msg) = @printf("  ok    %s\n", msg)

"Compare x against a reference value with a relative tolerance."
function cmp_rel(label, x, want, tol)
    if !isfinite(x)
        fail("$label is $x, expected $want")
        return
    end
    rel = abs(x - want) / max(abs(want), eps())
    rel <= tol ? pass(@sprintf("%s = %.6g (ref %.6g, %.1e rel)", label, x, want, rel)) :
                 fail(@sprintf("%s = %.6g, reference %.6g, %.2e rel > %.1e",
                               label, x, want, rel, tol))
end

load_or(path) = isfile(path) ? load(path) : nothing

function check_occam(site)
    f = joinpath(OUT, "aem_occam_$(site).jld2")
    A = load_or(f)
    A === nothing && return fail("missing $(basename(f))")
    r = get(get(REF, "deterministic", Dict()), site, nothing)
    if r === nothing
        @printf("  --    %s: no deterministic reference, structure only\n", site)
    else
        tol = Float64(get(r, "tol", 1e-4))
        cmp_rel("$site chi2/N", Float64(A["chi2n"]), Float64(r["chi2n"]), tol)
        Int(A["ndata"]) == Int(r["ndata"]) ?
            pass("$site ndata = $(A["ndata"])") :
            fail("$site ndata = $(A["ndata"]), reference $(r["ndata"])")
        # Geometric means over depth windows, exactly as the stage reports
        # them: 10^mean(log10 sigma) over the cells in [za, zb]. NOT the
        # arithmetic mean -- the two differ by about 1.6x in this ground, so
        # the wrong one would look like a real disagreement.
        if haskey(r, "sigma_windows_m")
            z = Float64.(A["zall"]); m = Float64.(A["m_log10sigma"])
            wins = r["sigma_windows_m"]      # [[0, 20], [20, 60], [60, 120]]
            for (w, want) in zip(wins, Float64.(r["sigma_geomean"]))
                za, zb = Float64(w[1]), Float64(w[2])
                k = findall(zz -> za <= zz <= zb, z)
                if isempty(k)
                    fail(@sprintf("%s window %.0f-%.0f m contains no cells", site, za, zb))
                    continue
                end
                cmp_rel(@sprintf("%s sigma %.0f-%.0f m", site, za, zb),
                        10^mean(m[k]), want, tol)
            end
        end
    end
    # True of any Occam solve, reference or not.
    all(isfinite, A["sigma"]) || fail("$site Occam sigma(z) has non-finite values")
    all(>(0), A["sigma"])     || fail("$site Occam sigma(z) is not positive")
end

function check_transd(site)
    for (tag, f) in (("AEM trans-D", "aem_td_$(site).jld2"),
                     ("SNMR trans-D", "smr_transd_$(site).jld2"))
        p = joinpath(OUT, f)
        isfile(p) ? pass("$site $tag present ($(round(filesize(p)/1024, digits=1)) kB)") :
                    fail("$site $tag missing: $f")
    end
end

function check_detect(site)
    D = load_or(joinpath(OUT, "detect_$(site).jld2"))
    D === nothing && return fail("$site detect output missing")
    za, zb = Float64.(D["win"])
    zgrid = Float64.(D["z"])
    if !(isfinite(za) && isfinite(zb))
        fail(@sprintf("%s target interval is not finite (%.3g-%.3g)", site, za, zb))
    elseif !(za < zb)
        fail(@sprintf("%s target interval is not ordered (%.2f-%.2f m)", site, za, zb))
    elseif za < minimum(zgrid) - 1e-9 || zb > maximum(zgrid) + 1e-9
        fail(@sprintf("%s target interval %.2f-%.2f m is outside the grid %.2f-%.2f m",
                      site, za, zb, minimum(zgrid), maximum(zgrid)))
    else
        pass(@sprintf("%s target interval %.2f-%.2f m, inside the grid", site, za, zb))
    end
    v = String(D["verdict"])
    if isempty(DETECT_PREFIXES) || any(p -> startswith(v, p), DETECT_PREFIXES)
        pass("$site detect verdict \"$v\"")
    else
        fail("$site detect verdict \"$v\" starts with none of the known clauses: " *
             join(DETECT_PREFIXES, " | "))
    end
    isfinite(Float64(D["zlim"])) || fail("$site zlim is not finite")
end

function check_profile(site)
    P = load_or(joinpath(OUT, "aem_profile_$(site).jld2"))
    P === nothing && return fail("$site profile output missing")
    lo, hi = Float64(P["edge_lo"]), Float64(P["edge_hi"])
    if !(isfinite(lo) && isfinite(hi) && lo > 0)
        fail(@sprintf("%s conductivity band is not positive and finite (%.4g, %.4g)", site, lo, hi))
    elseif !(lo < hi)
        fail(@sprintf("%s conductivity band is not ordered (lo %.5f >= hi %.5f S/m)", site, lo, hi))
    else
        pass(@sprintf("%s band %.5f-%.5f S/m (%.2f decades)", site, lo, hi, log10(hi/lo)))
    end
end

function check_gate(site)
    G = load_or(joinpath(OUT, "gate_$(site).jld2"))
    G === nothing && return fail("$site gate output missing")
    σlo, σhi = Float64(G["sigma_lo"]), Float64(G["sigma_hi"])
    σlo < σhi || fail(@sprintf("%s gate band is not ordered (%.5f, %.5f)", site, σlo, σhi))
    x = Float64(G["Xcert_p95"])
    isfinite(x) && x > 0 ? pass(@sprintf("%s X_cert p95 = %.4f (%.0f mg/L)", site, x, x * Float64(G["mgl"]))) :
                           fail("$site X_cert p95 is $x")
    v = String(G["verdict"])
    if isempty(GATE_VERDICTS) || v in GATE_VERDICTS
        pass("$site gate verdict \"$v\"")
    else
        fail("$site gate verdict \"$v\" is not one of " *
             join(sort(collect(GATE_VERDICTS)), " | "))
    end
end

function main()
    @printf("reference check: %s\n", basename(CFGPATH))
    @printf("  outputs   %s\n", OUT)
    @printf("  reference %s\n", REFPATH)
    @printf("  mode      %s\n\n", MODE)

    want_mode = String(get(REF, "mode", MODE))
    want_mode == MODE ||
        @printf("NOTE: reference was recorded in mode \"%s\", this run is \"%s\".\n\n",
                want_mode, MODE)

    if MODE != "smoke" && !haskey(REF, "full")
        println("NOTE: this is a $MODE run but the reference has no [full] section.")
        println("      Only the Occam numbers and the structural invariants are checked.")
        println("      Converged reference values have not been recorded yet.\n")
    end

    for s in SITES
        @printf("%s\n", s)
        "aem_occam"   in STAGES && check_occam(s)
        ("aem_transd" in STAGES || "snmr_transd" in STAGES) && check_transd(s)
        "detect"      in STAGES && check_detect(s)
        "profile"     in STAGES && check_profile(s)
        "gate"        in STAGES && check_gate(s)
        println()
    end

    if isempty(FAILS)
        @printf("PASS: %d sites agree with %s\n", length(SITES), basename(REFPATH))
        exit(0)
    end
    @printf("FAIL: %d check(s) disagree with the reference\n", length(FAILS))
    for f in FAILS
        println("  - ", f)
    end
    exit(1)
end

main()
