# stage6_verdict.jl -- Stage 6. The verdict table.
#
#   julia --project=. stage6_verdict.jl --config ../configs/synth_3P.toml
#
# Exactly the rows of the paper's illustrative table, in its order, machine-
# readable (CSV + JSON) and rendered (Markdown). There is deliberately NO
# "next measurement" row.
#
# VERDICTS ARE THREE-STATE. PASS, FAIL and UNDETERMINED are different
# statements and the third is not evidence against the hypothesis. The field
# work made this concrete: on 19 soundings the lower edge of the admissible
# band ran to the sweep floor at every one, and that was being read as physics
# when it was the edge of the search. An UNDETERMINED verdict says the data do
# not decide; a FAIL says they decide against.
#
# A SATURATED CRITERION IS NOT A PASS EARNED. Where Stage 3 reports that every
# window the interval spec allows clears the water criterion, the water row is
# marked PASS (uninformative) -- the hypothesis was satisfied by the spec, not
# by the ground.
#
# THIS IS THE ONLY STAGE THAT MAY SCORE AGAINST TRUTH, and only after every
# verdict above is fixed. The scoring block reads the truth file, compares, and
# cannot feed back into any verdict because the verdicts are already written.
using JLD2, Statistics, Printf, Dates
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "config.jl")); using .SiteConfig

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) ? (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (i += 1)
    end
end
const CFGPATH = get(ARG, "config", "")
isempty(CFGPATH) && error("usage: stage6_verdict.jl --config <path>")
f3(x) = isfinite(x) ? @sprintf("%.3f", x) : "n/a"
f4(x) = isfinite(x) ? @sprintf("%.4g", x) : "n/a"
f0(x) = isfinite(x) ? @sprintf("%.0f", x) : "n/a"
jsonesc(s) = replace(string(s), "\\"=>"\\\\", "\""=>"\\\"", "\n"=>"\\n")

function main()
    cfg  = SiteConfig.resolve(CFGPATH)
    info = SiteConfig.assertvalid(cfg; root = ROOT)
    d, hy, out, pr = cfg["data"], cfg["hypothesis"], cfg["output"], cfg["priors"]
    odir = joinpath(ROOT, out["dir"])
    S3 = load(joinpath(odir, "stage3_$(cfg["name"]).jld2"))
    P4 = load(joinpath(odir, "stage4post_$(cfg["name"]).jld2"))
    S5 = load(joinpath(odir, "stage5_$(cfg["name"]).jld2"))

    iv = hy["interval"]; wv = hy["water"]; sv = hy["salinity"]
    Xlim = Float64(sv["value"]); fresh_is_less = sv["direction"] == "less_than"
    q(v,p) = quantile(v, p)

    # ---------- the rows ----------------------------------------------
    rows = Vector{Tuple{String,String,String}}()   # (row, value, note)
    push!(rows, ("Hypothesis",
        @sprintf("interval L %.0f-%.0f m within %.0f-%.0f m; W >= %.1f m; TDS %s %.0f mg/L",
                 iv["L_min"], iv["L_max"], iv["depth_min"], iv["depth_max"],
                 Float64(wv["value"]), fresh_is_less ? "<" : ">", Xlim),
        "pre-registered"))

    tops, bases, lens = S3["tops"], S3["bases"], S3["lengths"]
    multi = haskey(S3, "P_perdraw")  # always true; modes recomputed below
    srt = sort(tops); gaps = diff(srt); bigg = findall(>(3.0), gaps)
    intnote = isempty(bigg) ? "unimodal" :
        "MULTIMODAL: " * join([begin
            edges = vcat(1, bigg .+ 1, length(srt)+1)
            @sprintf("%.1f-%.1f m (%.0f%%)", srt[edges[b]], srt[edges[b+1]-1],
                     100*(edges[b+1]-edges[b])/length(srt))
        end for b in 1:length(bigg)+1], " / ")
    push!(rows, ("Interval: top",
        @sprintf("%.2f m  (p5-p95 %.2f - %.2f)", q(tops,.5), q(tops,.05), q(tops,.95)), intnote))
    push!(rows, ("Interval: base",
        @sprintf("%.2f m  (p5-p95 %.2f - %.2f)", q(bases,.5), q(bases,.05), q(bases,.95)), ""))
    push!(rows, ("Interval: length",
        @sprintf("%.2f m  (p5-p95 %.2f - %.2f)", q(lens,.5), q(lens,.05), q(lens,.95)),
        S3["length_at_Lmax"] > 0.5 ?
          @sprintf("PINNED at L_max in %.0f%% of draws: spec-determined, not inferred",
                   100*S3["length_at_Lmax"]) : ""))

    sat = get(S3, "saturated", false)
    fwp = get(S3, "frac_windows_pass", NaN)
    # THE DISCRIMINATING NUMBER, as its own row. P(W >= w) for the BEST window is
    # 1.0000 on every column including the thin confined one, so it tells a
    # reader nothing. The fraction of allowed windows that pass does: 54% on the
    # unconfined column against 8% on the confined one.
    push!(rows, ("Water criterion: windows passing",
        isfinite(fwp) ? @sprintf("%.1f%% of %d allowed windows", 100fwp, get(S3,"n_windows",0)) : "n/a",
        "the discriminating statistic; P for the best window saturates at 1.0 on every column"))
    push!(rows, ("Water criterion",
        isfinite(S3["P_perdraw"]) ?
            @sprintf("P(W >= %.1f m) = %.4f", Float64(wv["value"]), S3["P_perdraw"]) :
            "not determined",
        isfinite(fwp) ? @sprintf("%.1f%% of %d allowed windows satisfy it%s",
                                 100fwp, get(S3,"n_windows",0),
                                 sat ? "; SATURATED -- not a pass earned" : "") : ""))

    bl, bh = P4["band_C"].lo, P4["band_C"].hi
    bo = P4["band_open"]
    push!(rows, ("Admissible conductance band",
        @sprintf("%s - %s S, width %s dex", f4(bl), f4(bh), f3(P4["band_dex"])),
        @sprintf("lower edge open = %s; upper edge open = %s", bo.lo, bo.hi)))

    push!(rows, ("Salinity criterion",
        isfinite(P4["P_wrong"]) ?
            @sprintf("P(TDS %s %.0f mg/L) = %.4f,  Delta = %s",
                     fresh_is_less ? ">" : "<", Xlim, P4["P_wrong"], f4(P4["delta"])) :
            "not determined (no admissible band)",
        @sprintf("%s (tolerance %.2f)", P4["verdict_label"], P4["delta_tol"])))

    push!(rows, ("Data-only bound",
        @sprintf("TDS at 95%% = %s mg/L", f0(P4["tds_data_only_p95"])),
        "clay-free, from the band's upper edge; an UPPER bound on salinity"))

    pp = S5["priors"]
    # ACCEPTANCE RATE as its own row. Draws with clay_product >= C/thickness
    # give a negative sigma_w: the clay alone would explain everything
    # measured. A low acceptance rate means the Stage 5 estimate is carried by
    # the prior's tail, so it is stated rather than buried.
    acc = get(S5, "acceptance", NaN)
    push!(rows, ("Stage 5 acceptance rate",
        isfinite(acc) ? @sprintf("%.1f%% of draws give sigma_w > 0", 100acc) : "n/a",
        !isfinite(acc) ? "no draws (undetermined)" :
        acc < 0.10 ? @sprintf("BELOW 10%%: Stage 5 estimate is FRAGILE -- it rests on the clay prior's tail") :
                     "clay prior is consistent with the measured conductance"))
    push!(rows, ("Stage 5 estimate",
        isfinite(S5["TDS_p50"]) ?
            @sprintf("TDS p50 %s mg/L  (p5-p95 %s - %s)", f0(S5["TDS_p50"]),
                     f0(S5["TDS_p5"]), f0(S5["TDS_p95"])) :
            "not determined",
        @sprintf("m~U(%.2f,%.2f), phi~U(%.2f,%.2f), clay product ~ log-U(%.0e,%.0e) S/m%s",
                 pp.m[1], pp.m[2], pp.porosity[1], pp.porosity[2],
                 pp.clay_product[1], pp.clay_product[2],
                 get(S5, "fragile", false) ? "; FRAGILE" : "")))

    # ---------- the verdict -------------------------------------------
    # Three-state, per criterion. UNDETERMINED wherever the quantity the
    # criterion needs was never determined -- an open band edge, a
    # prior-carried Delta, or a saturated criterion whose pass is the spec's.
    Pw3 = S3["P_perdraw"]
    v_water = !isfinite(Pw3) ? "UNDETERMINED" :
              Pw3 >= 0.95 ? (sat ? "PASS (uninformative)" : "PASS") :
              Pw3 <= 0.05 ? "FAIL" : "UNDETERMINED"
    # THE SALINITY VERDICT IS ONE-SIDED, AND THIS IS NOT A CHOICE.
    #
    # Stage 4's P comes from a CLAY-FREE conversion: all the measured
    # conductance is attributed to pore water, so the salinity it returns is an
    # UPPER BOUND on the true salinity, biased in one direction by construction.
    # A one-sided bound cannot support a two-sided verdict:
    #
    #   clay-free says FRESH  ->  the water IS fresh. Even the most pessimistic
    #                             attribution clears the limit. PASS is sound.
    #   clay-free says SALTY  ->  it may be salt, or it may be clay, and the AEM
    #                             cannot tell. That is UNDETERMINED, not FAIL.
    #
    # Reading a high clay-free P(too salty) as FAIL is how this table first
    # rejected column 2P, whose true water is 561 mg/L against a 1,000 mg/L
    # limit -- a confident FAIL on a hypothesis that is true. The clay-free
    # number was correct as a bound (59,926 mg/L does bound 561 mg/L); using it
    # as a verdict was the error.
    #
    # A FAIL therefore requires the CLAY-SUBTRACTED Stage 5 estimate to put at
    # least 95% of its mass on the wrong side, and that is flagged as
    # prior-dependent wherever it fires.
    # STAGE 4 FLOOR UNDER STATED CEILINGS -- the only route to a FAIL that does
    # not go through Stage 5. Cap every term that could be explaining the
    # measured conductance instead of salt (clay at [ceilings].sigma_clay_max
    # times C_max, porosity at [ceilings].porosity_max, m at its prior's lower
    # bound) and whatever conductance is left must be water. If that leftover
    # alone puts the water past the line, the hypothesis fails -- but only
    # under those ceilings, which is why the label says so and why they live
    # in the config. A floor needs a CLOSED lower edge: an open one sits at the
    # solver's sigma_min, not at the data, so there is nothing to subtract from
    # and stage4_floor.jl returns zero.
    floor_tds = 0.0
    let ff = joinpath(ROOT, "outputs", "stage4_floor", "floor.jld2"),
        tag = replace(cfg["name"], "synth_" => "")
        if isfile(ff)
            D = load(ff)
            haskey(D, tag) && (floor_tds = D[tag].floor_tds)
        end
    end
    floor_fails = isfinite(floor_tds) && floor_tds > Float64(sv["value"])

    P_fresh = isfinite(P4["P_wrong"]) ? 1 - P4["P_wrong"] : NaN
    P5 = get(S5, "P_criterion", NaN)          # Stage 5: P(criterion MET), clay subtracted
    v_sal_note = ""
    v_sal = if !isfinite(P_fresh) || !isfinite(P4["delta"])
                v_sal_note = "no admissible band: the gate did not determine a conductance"
                "UNDETERMINED"
            elseif floor_fails
                v_sal_note = @sprintf("the floor under the stated ceilings is %.0f mg/L, past the %.0f mg/L line; ceilings are in config.toml [ceilings]",
                                      floor_tds, Float64(sv["value"]))
                "FAIL under stated ceilings"
            elseif isfinite(P5) && P5 <= 0.05 && !get(S5, "fragile", false)
                v_sal_note = "FAIL rests on the clay-subtracted Stage 5 estimate, so it is PRIOR-DEPENDENT"
                "FAIL"
            elseif isfinite(P5) && P5 <= 0.05
                # Stage 5 is the only route to a FAIL, and a FRAGILE Stage 5 is
                # not a route at all: below 10% acceptance the estimate is the
                # clay prior's tail, not the ground.
                v_sal_note = @sprintf("Stage 5 says salty but its acceptance rate is %.1f%%: FRAGILE, so no FAIL", 100acc)
                "UNDETERMINED"
            elseif P_fresh >= 0.95 && P4["delta"] <= P4["delta_tol"]
                "PASS"
            elseif P_fresh >= 0.95
                v_sal_note = @sprintf("clay-free bound clears the limit but Delta = %.4f exceeds %.2f: not robust to lambda2 or to dropping an end node",
                                      P4["delta"], P4["delta_tol"])
                "UNDETERMINED"
            else
                v_sal_note = "clay-free bound cannot clear the limit; clay and salt are not separable here, so this is NOT a FAIL"
                "UNDETERMINED"
            end
    v_int = isempty(bigg) ? "PASS" : "PASS (multimodal)"
    overall = any(v -> startswith(v, "FAIL"), (v_water, v_sal)) ? v_sal :
              any(v -> startswith(v, "UNDETERMINED"), (v_water, v_sal, v_int)) ? "UNDETERMINED" :
              "PASS"
    fl = get(S5, "flips", String[])
    push!(rows, ("Inputs that flip the verdict alone",
        isempty(fl) ? "none" : join(fl, ", "),
        "swing table: each input at p10/p90, others at their medians"))
    push!(rows, ("Verdict: interval", v_int, ""))
    push!(rows, ("Verdict: water", v_water, ""))
    push!(rows, ("Stage 4 floor under stated ceilings",
        @sprintf("TDS >= %.0f mg/L", floor_tds),
        floor_tds > 0 ? "clay, porosity and m all at their config ceilings" :
                        "zero: the band's lower edge is OPEN, so there is no measured minimum conductance"))
    push!(rows, ("Verdict: salinity", v_sal, v_sal_note))
    push!(rows, ("Verdict: OVERALL", overall, ""))

    # ---------- render -------------------------------------------------
    println()
    println("="^96)
    @printf("VERDICT TABLE -- %s   (%s)\n", cfg["name"], cfg["label"])
    @printf("lineage: column %s, parent %s, role %s   |   config %s\n",
            d["column"], isempty(info.parent) ? "-" : info.parent, info.role,
            basename(cfg["_source"]))
    println("="^96)
    @printf("%-28s  %-44s  %s\n", "ROW", "VALUE", "NOTE")
    println("-"^96)
    for (r, v, nt) in rows; @printf("%-28s  %-44s  %s\n", r, v, nt); end
    println("="^96)

    mkpath(odir)
    # CSV
    open(joinpath(odir, "verdict_$(cfg["name"]).csv"), "w") do f
        println(f, "site,column,parent,role,row,value,note")
        for (r, v, nt) in rows
            println(f, join(["\"$(cfg["name"])\"", "\"$(d["column"])\"",
                "\"$(get(d,"parent",""))\"", "\"$(info.role)\"",
                "\"$r\"", "\"$(replace(v, "\""=>"\"\""))\"",
                "\"$(replace(nt, "\""=>"\"\""))\""], ","))
        end
    end
    # JSON
    open(joinpath(odir, "verdict_$(cfg["name"]).json"), "w") do f
        println(f, "{")
        println(f, "  \"site\": \"", jsonesc(cfg["name"]), "\",")
        println(f, "  \"label\": \"", jsonesc(cfg["label"]), "\",")
        println(f, "  \"column\": \"", jsonesc(d["column"]), "\",")
        println(f, "  \"parent\": \"", jsonesc(get(d,"parent","")), "\",")
        println(f, "  \"role\": \"", jsonesc(info.role), "\",")
        println(f, "  \"config\": \"", jsonesc(cfg["_source"]), "\",")
        println(f, "  \"generated\": \"", Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"), "\",")
        println(f, "  \"verdict_overall\": \"", overall, "\",")
        println(f, "  \"rows\": [")
        for (i,(r,v,nt)) in enumerate(rows)
            println(f, "    {\"row\": \"", jsonesc(r), "\", \"value\": \"", jsonesc(v),
                    "\", \"note\": \"", jsonesc(nt), "\"}", i < length(rows) ? "," : "")
        end
        println(f, "  ]")
        println(f, "}")
    end
    # Markdown
    open(joinpath(odir, "verdict_$(cfg["name"]).md"), "w") do f
        println(f, "# Verdict table — ", cfg["name"], "\n")
        println(f, "**", cfg["label"], "**  \n")
        println(f, "Column `", d["column"], "`, parent `", get(d,"parent","—"),
                "`, role `", info.role, "`. Config `", basename(cfg["_source"]), "`.\n")
        println(f, "| Row | Value | Note |")
        println(f, "|---|---|---|")
        for (r,v,nt) in rows; println(f, "| ", r, " | ", v, " | ", nt, " |"); end
    end
    println("wrote verdict_$(cfg["name"]).{csv,json,md} in ", out["dir"])

    # ---------- scoring against truth, AFTER the verdicts are fixed ----
    if get(out, "score_against_truth", false) && haskey(d, "truth")
        tp = isfile(joinpath(ROOT, d["truth"])) ? joinpath(ROOT, d["truth"]) : d["truth"]
        if isfile(tp)
            println("\nSCORE AGAINST TRUTH (read only now; every verdict above is already fixed)")
            T = load(tp)
            zt, sbt, swt = T["z"], T["σ_bulk"], T["σw"]
            dzt = diff(vcat(0.0, collect(Float64, zt)))
            wt, wb = P4["window"].top, P4["window"].base
            it = findall(i -> zt[i] > wt && zt[i] <= wb, eachindex(zt))
            Ctrue = sum(sbt[it] .* dzt[it])
            swtrue = mean(swt[it]) * P4["mg_per_sm"]
            hasband = isfinite(bl) && isfinite(bh)
            inband = hasband ? (bl <= Ctrue <= bh) : false
            @printf("  true conductance over the reported window = %.4f S  (sigma_eff %.5f S/m)\n",
                    Ctrue, Ctrue / P4["thick"])
            @printf("  reported band %s - %s S  ->  %s\n", f4(bl), f4(bh),
                    hasband ? "truth INSIDE: $(inband)" : "no band to score against")
            @printf("  true water TDS in the window = %.0f mg/L   (criterion %s %.0f)\n",
                    swtrue, fresh_is_less ? "<" : ">", Xlim)
            truth_fresh = fresh_is_less ? (swtrue < Xlim) : (swtrue > Xlim)
            @printf("  truth satisfies the salinity criterion: %s   | verdict said: %s\n",
                    truth_fresh, v_sal)
            if !isfinite(P4["tds_data_only_p95"])
                println("  data-only bound: not defined (no admissible band) -- nothing to score")
            else
                @printf("  data-only bound %s mg/L %s the truth (%s)\n",
                        f0(P4["tds_data_only_p95"]),
                        P4["tds_data_only_p95"] >= swtrue ? "BOUNDS" : "FAILS TO BOUND",
                        P4["tds_data_only_p95"] >= swtrue ? "correct: it is an upper bound" :
                            "WRONG: a clay-free bound must not fall below the truth")
            end
        end
    end
end
main()
