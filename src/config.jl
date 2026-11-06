# config.jl -- load a site config, resolve `extends`, and assert the
# invariants the schema promises. Every stage receives the resolved Dict
# and reads nothing else.
module SiteConfig

using TOML

export loadconfig, resolve, assertvalid, LITHOLOGY_M

const SCHEMA_SUPPORTED = (1,)

# The cementation-exponent prior, keyed to lithology. Ranges and their
# references are tabulated in docs/choosing_priors.md.
const LITHOLOGY_M = Dict(
    "clean_sand"             => (1.3, 1.6),
    "clayey_sand"            => (1.5, 1.9),
    "consolidated_sandstone" => (1.8, 2.1),
)

"""Recursive per-key merge. Arrays REPLACE, never append, so an override
cannot silently extend a ladder or a threshold list."""
function mergecfg(base::AbstractDict, over::AbstractDict)
    out = Dict{String,Any}(base)
    for (k, v) in over
        out[k] = (v isa AbstractDict && get(base, k, nothing) isa AbstractDict) ?
                 mergecfg(base[k], v) : v
    end
    out
end

"""Load `path`, following `extends` chains. Paths inside a file are
relative to THAT file. Cycles abort rather than recurse."""
function resolve(path::AbstractString; seen = String[])
    ap = abspath(path)
    ap in seen && error("config `extends` cycle: " * join([seen; ap], " -> "))
    isfile(ap) || error("config not found: $ap")
    cfg = TOML.parsefile(ap)
    cfg["_source"] = ap
    if haskey(cfg, "extends")
        parent = abspath(joinpath(dirname(ap), cfg["extends"]))
        base = resolve(parent; seen = [seen; ap])
        delete!(cfg, "extends")
        src = cfg["_source"]
        cfg = mergecfg(base, cfg)
        cfg["_source"] = src          # the leaf file, not the base
        cfg["_base"] = parent
    end
    cfg
end

"""Every check that can be made before a single forward model is run.
Fails loudly: a config that cannot be trusted must not reach Stage 1."""
function assertvalid(c::AbstractDict; root = pwd(), require_files = true)
    probs = String[]
    push_if(cond, msg) = cond || push!(probs, msg)

    sv = get(c, "schema_version", nothing)
    push_if(sv in SCHEMA_SUPPORTED,
            "schema_version $(repr(sv)) not supported; this code knows $(SCHEMA_SUPPORTED)")

    for s in ("data", "priors", "hypothesis", "gate", "inversion", "output")
        push_if(haskey(c, s), "missing section [$s]")
    end
    isempty(probs) || error("config invalid:\n  " * join(probs, "\n  "))

    d, pr, hy, g, inv = c["data"], c["priors"], c["hypothesis"], c["gate"], c["inversion"]

    # -- Stage 0 is the only stage that may name a truth file ----------
    synthetic = haskey(d, "truth") && !isempty(get(d, "truth", ""))
    push_if(synthetic ⊻ (haskey(d, "aem_file") && !isempty(get(d, "aem_file", ""))),
            "exactly one of [data].truth (synthetic) or [data].aem_file (field) must be set")

    # -- loop geometry: the square/circular flag has been inverted here
    #    before, so shape+side must reproduce the recorded area --------
    if haskey(d, "loop")
        L = d["loop"]
        shape, side = get(L, "shape", ""), get(L, "side", NaN)
        area = shape == "square"   ? side^2 :
               shape == "circular" ? pi * side^2 : NaN
        want = get(L, "area_check", NaN)
        push_if(shape in ("square", "circular"), "[data.loop].shape must be square|circular")
        push_if(isfinite(area) && isfinite(want) && abs(area - want) <= 1e-6 * want,
                "[data.loop]: shape=$shape side=$side gives area " *
                "$(repr(area)) m^2, but area_check = $(repr(want)) m^2 " *
                "(relative difference $(abs(area-want)/want)) -- one of them is wrong")
    end

    # -- the m prior: keyword or explicit range, not both, not neither -
    haskw, hasrange = haskey(pr, "m_lithology"), haskey(pr, "m_range")
    push_if(haskw ⊻ hasrange, "[priors]: set exactly one of m_lithology or m_range")
    if haskw
        push_if(haskey(LITHOLOGY_M, pr["m_lithology"]),
                "[priors].m_lithology = $(repr(pr["m_lithology"])) unknown; " *
                "known: $(join(sort(collect(keys(LITHOLOGY_M))), ", "))")
    end
    mlo, mhi = haskw && haskey(LITHOLOGY_M, pr["m_lithology"]) ?
               LITHOLOGY_M[pr["m_lithology"]] : (hasrange ? Tuple(pr["m_range"]) : (NaN, NaN))
    push_if(isfinite(mlo) && mlo < mhi, "[priors]: m prior must have lo < hi, got ($mlo, $mhi)")

    # -- porosity ------------------------------------------------------
    if haskey(pr, "porosity_range")
        p1, p2 = pr["porosity_range"]
        push_if(0 < p1 < p2 < 1, "[priors].porosity_range must satisfy 0 < lo < hi < 1, got $([p1,p2])")
    end
    if haskey(pr, "clay")
        cl = pr["clay"]
        push_if(haskey(cl, "product"),
                "[priors.clay] needs `product` = [lo, hi] S/m, a log-uniform prior on " *
                "sigma_clay*C. The separate fraction/sigma_clay form is superseded: it " *
                "concentrated mass above the measured conductance (98.1% rejected on 3P)")
        if haskey(cl, "product")
            c1, c2 = cl["product"]
            push_if(0 < c1 < c2, "[priors.clay].product must satisfy 0 < lo < hi")
        end
        for dead in ("fraction_dist", "fraction_median", "logit_sd", "correlation_m", "sigma_clay")
            push_if(!haskey(cl, dead),
                    "[priors.clay].$dead is from the superseded prior; remove it " *
                    "(use `product` instead) so no run silently mixes the two")
        end
    end

    # -- the hypothesis must be BOUNDED --------------------------------
    if haskey(hy, "interval")
        iv = hy["interval"]
        if get(iv, "mode", "") == "search"
            for k in ("L_min", "L_max", "depth_min", "depth_max")
                push_if(haskey(iv, k), "[hypothesis.interval]: mode=search needs $k")
            end
            if all(haskey(iv, k) for k in ("L_min", "L_max", "depth_min", "depth_max"))
                push_if(0 < iv["L_min"] <= iv["L_max"],
                        "[hypothesis.interval]: need 0 < L_min <= L_max")
                push_if(0 <= iv["depth_min"] < iv["depth_max"],
                        "[hypothesis.interval]: need depth_min < depth_max")
                # An unbounded screen pins to the base of the grid, where
                # the depth limit has no teeth.
                push_if(iv["L_max"] <= iv["depth_max"] - iv["depth_min"],
                        "[hypothesis.interval]: L_max $(iv["L_max"]) m does not fit inside " *
                        "the depth range $(iv["depth_min"])-$(iv["depth_max"]) m")
            end
        end
    end

    # -- the gate grid must span the solver's own bounds ---------------
    if haskey(g, "sigma_bounds")
        b1, b2 = g["sigma_bounds"]
        push_if(0 < b1 < b2, "[gate].sigma_bounds must satisfy 0 < lo < hi")
        # Nodes near the bounds are unreachable; the schema's answer is to
        # exclude them per node, which only works if that is switched on.
        push_if(get(g, "exclude_on_bound", false),
                "[gate].exclude_on_bound = false with a grid that runs to the model " *
                "bounds: the end nodes cannot satisfy the constraint and must not be " *
                "normalised over (set it true, or narrow the grid deliberately)")
    end

    # -- Stage 2's per-draw kernel requires within-kernel statistics ----
    if haskey(inv, "snmr")
        sn = inv["snmr"]
        if get(sn, "kernel", "") == "per_draw"
            push_if(get(sn, "statistics_scope", "") == "within_kernel",
                    "[inversion.snmr]: kernel=per_draw requires statistics_scope=within_kernel " *
                    "-- pooling draw-to-data statistics across kernels is what broke three " *
                    "field diagnostics")
        end
        push_if(get(sn, "hierarchical_noise", false),
                "[inversion.snmr].hierarchical_noise = false: the records are effectively " *
                "noiseless under a fixed noise model and the deep posterior collapses to " *
                "the GP mean")
    end

    # -- lineage: a control is meaningless without its parent ----------
    role = get(d, "role", "primary")
    push_if(role in ("primary", "saline_control"),
            "[data].role must be primary|saline_control")
    if role == "saline_control"
        push_if(!isempty(get(d, "parent", "")),
                "[data]: a saline_control must name the [data].parent it was built from " *
                "-- comparing a control against the wrong parent has already produced one " *
                "retracted finding in this project")
    end

    # -- files exist ---------------------------------------------------
    # Two kinds, and conflating them hid a defect: `generator` and `system` are
    # files the repository ships, so they must be there even on a dry run, while
    # `truth` and `snmr_ensemble` are built by earlier stages and legitimately
    # absent. A dry run used to skip both and report "config valid" for a
    # config naming a generator script that had been deleted.
    for k in ("generator", "system")
        v = get(d, k, "")
        (v isa String && !isempty(v)) || continue
        push_if(isfile(joinpath(root, v)) || isfile(v),
                "[data].$k not found: $v (a file this repository should ship)")
    end
    if require_files
        # These are products of earlier stages, so "not found" on a clean clone
        # is the normal case and the message has to say what produces them. A
        # review hit the README's own quick-start line and got two paths and no
        # remedy.
        made_by = Dict(
            "truth" => "synthetics/scripts/build_truth.jl, then copied to data/truth_murray<T>.jld2",
            "snmr_ensemble" => "stages/stage2_snmr.jl (which needs tools/smr_data.jl first)")
        for k in ("truth", "snmr_ensemble")
            v = get(d, k, "")
            (v isa String && !isempty(v)) || continue
            push_if(isfile(joinpath(root, v)) || isfile(v),
                    "[data].$k not found: $v\n      produced by $(made_by[k])" *
                    "\n      run run_site.jl with --dry-run to see the whole chain")
        end
    end

    isempty(probs) || error("config invalid ($(get(c, "_source", "?"))):\n  " *
                            join(probs, "\n  "))
    (; m_prior = (mlo, mhi), synthetic, role, parent = get(d, "parent", ""))
end

loadconfig(path; kw...) = (c = resolve(path); (c, assertvalid(c; kw...)))

end # module
