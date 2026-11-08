# stage0_inputs.jl -- Stage 0 of the Section 3 workflow.
#
#   julia --project=. stage0_inputs.jl --config ../configs/synth_3P.toml
#
# THE ONLY STAGE PERMITTED TO READ A TRUTH FILE.
#
# It reads the truth column, forward-models the AEM data, adds noise, and then
# writes a SEALED ARTEFACT holding only what a field crew would come home with:
# the observed data, the noise estimates, the system geometry, the depth grid,
# and the settings. Stages 1 to 6 open that artefact and the config, and
# nothing else.
#
# WHY THIS EXISTS. On synthetic data the easiest way to flatter a workflow is
# to let a later stage peek at the answer, and it has happened repeatedly in
# this project: screens placed at the truth's contacts, the cementation
# exponent set to the truth's 1.7, a conductance integrated over the interval
# the truth calls the aquifer. None of those were deliberate; each one was a
# convenient `T["..."]` three functions deep. A firewall that is a rule gets
# broken. A firewall that is a FILE BOUNDARY does not, because the later
# stages have nothing to read.
#
# The artefact is written against a WHITELIST, not a blacklist: only the keys
# named in SEALED_KEYS are allowed out, so adding a field to the truth file
# cannot silently widen what Stages 1-6 can see. The whitelist is asserted
# before the file is written.
using JLD2, Statistics, Printf, Random, TOML, LinearAlgebra
const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "config.jl")); using .SiteConfig
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "aem_jac.jl"))

const ARG = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        startswith(ARGS[i], "--") && i < length(ARGS) ? (ARG[ARGS[i][3:end]] = ARGS[i+1]; i += 2) : (i += 1)
    end
end
const CFGPATH = get(ARG, "config", "")
isempty(CFGPATH) && error("usage: stage0_inputs.jl --config <path/to/site.toml>")

# Everything Stages 1-6 are allowed to see. Nothing else leaves this stage.
const SEALED_KEYS = Set([
    # --- the observation -------------------------------------------------
    "d_low", "d_high", "sigma_low", "sigma_high", "select_low", "select_high",
    "ndata_low", "ndata_high", "ndata",
    # --- the system, as acquired -----------------------------------------
    "times_low", "times_high", "ramp_low", "ramp_high",
    "r_tx", "r_rx", "z_tx", "z_rx", "lowpass_fc",
    # --- the inversion grid (an acquisition/processing choice) -----------
    "z", "zall", "znall", "zbnd", "nfixed",
    # --- SNMR acquisition + where its ensemble is ------------------------
    "snmr_ensemble", "snmr_dead_time", "snmr_larmor_hz", "snmr_inclination_deg",
    "snmr_noise_tier", "loop_shape", "loop_side_m", "loop_turns",
    # --- provenance (never values) ---------------------------------------
    "site_name", "site_label", "column", "parent", "role",
    "config_source", "config_hash", "seed", "m_prior", "stage0_written",
    "truth_file_sha",       # the IDENTITY of the truth, not its contents
])

function main()
    cfg  = SiteConfig.resolve(CFGPATH)
    info = SiteConfig.assertvalid(cfg; root = ROOT)
    d, out = cfg["data"], cfg["output"]

    @printf("Stage 0 -- inputs.  site %s  (%s)\n", cfg["name"], cfg["label"])
    @printf("  lineage: column %s, parent %s, role %s\n",
            d["column"], isempty(info.parent) ? "-" : info.parent, info.role)
    @printf("  m prior: U(%.2f, %.2f)  from lithology %s\n",
            info.m_prior[1], info.m_prior[2], get(cfg["priors"], "m_lithology", "explicit"))

    haskey(d, "truth") && !isempty(d["truth"]) ||
        error("field branch not implemented in this script; use stage0_field.jl")
    tpath = isfile(joinpath(ROOT, d["truth"])) ? joinpath(ROOT, d["truth"]) : d["truth"]

    # ---- the one truth read ------------------------------------------
    println("\n  READING TRUTH (the only stage that may): ", basename(tpath))
    T = load(tpath)
    seed = Int(d["seed"])

    # the loop geometry the config asserts must match the truth's acquisition
    if haskey(d, "loop") && haskey(T, "LOOP_SIDE")
        cs, ts = Float64(d["loop"]["side"]), Float64(T["LOOP_SIDE"])
        abs(cs - ts) <= 1e-9 || error("""
            [data.loop].side = $cs m but the truth file was built with
            LOOP_SIDE = $ts m. These are different surveys. Fix the config
            (or you are pointing at the wrong lineage: the P columns use
            100.0 m, the bare 3/2/3W/2W lineage uses 88.6 m).""")
    end

    B = build_aem_jac(T, seed)
    aem, grid = B.aem, B.grid
    nd = Int(aem.ndatalow) + Int(aem.ndatahigh)
    @printf("  forward + noise: %d gates (%d low, %d high), seed %d\n",
            nd, aem.ndatalow, aem.ndatahigh, seed)

    # The noise realisation is a property of the DATA and belongs here. The
    # trap: makenoisydata! reseeds from its own rseed, so every run that did
    # not pass it explicitly shares one realisation. aem_jac.jl passes it.
    #
    # Sanity check on the data: the GENERATING model must fit its own noisy
    # data at chi2/N ~ 1. `aem.res` is not populated by makenoisydata!, so it
    # has to be computed -- reading it straight after generation gives 0.0000
    # and looks like a perfect fit.
    # chi2 is the NOISE-WEIGHTED residual, norm(W*res)^2. Summing res alone
    # gives a number in raw data units (1e-12 V/Am^4) that reads as 0.0000.
    SK = HiQGA.transD_GP.SkyTEM1DInversion
    SK.getresidual(aem, log10.(B.σt[2:end]))
    chi2gen = norm(aem.W * aem.res)^2
    @printf("  generating model fits its own noisy data: chi2/N = %.4f  (must be ~1)\n",
            chi2gen / nd)
    (0.4 < chi2gen/nd < 2.5) || @printf("  WARNING: chi2/N = %.3f is not ~1 -- check the noise model\n", chi2gen/nd)

    # ---- the sealed artefact -----------------------------------------
    sys = Dict{String,Any}()
    A = Dict{String,Any}(
        "d_low" => copy(aem.dlow), "d_high" => copy(aem.dhigh),
        "sigma_low" => copy(aem.σlow), "sigma_high" => copy(aem.σhigh),
        "select_low" => copy(aem.selectlow), "select_high" => copy(aem.selecthigh),
        "ndata_low" => Int(aem.ndatalow), "ndata_high" => Int(aem.ndatahigh),
        "ndata" => nd,
        "times_low" => collect(LM_times), "times_high" => collect(HM_times),
        "ramp_low" => Matrix(LM_ramp), "ramp_high" => Matrix(HM_ramp),
        "r_tx" => rTx, "r_rx" => 13.0, "z_tx" => -40.0, "z_rx" => -42.0,
        "lowpass_fc" => collect(lowpassfcs),
        "z" => copy(aem.z), "zall" => collect(grid.zall),
        "znall" => collect(grid.znall), "zbnd" => collect(grid.zbnd),
        "nfixed" => Int(aem.nfixed),
        "snmr_ensemble" => get(d, "snmr_ensemble", ""),
        "snmr_dead_time" => get(get(d, "snmr", Dict()), "dead_time", NaN),
        "snmr_larmor_hz" => get(get(d, "snmr", Dict()), "larmor_hz", NaN),
        "snmr_inclination_deg" => get(get(d, "snmr", Dict()), "inclination_deg", NaN),
        "snmr_noise_tier" => get(get(d, "snmr", Dict()), "noise_tier", ""),
        "loop_shape" => get(get(d, "loop", Dict()), "shape", ""),
        "loop_side_m" => get(get(d, "loop", Dict()), "side", NaN),
        "loop_turns" => get(get(d, "loop", Dict()), "turns", 1),
        "site_name" => cfg["name"], "site_label" => cfg["label"],
        "column" => d["column"], "parent" => get(d, "parent", ""), "role" => info.role,
        "config_source" => cfg["_source"],
        # NOT a SHA, despite the key names: these are Julia`s 64-bit hash(),
        # which is not reproducible outside Julia and is not guaranteed stable
        # across Julia versions. A reviewer flagged both names as misleading.
        # Kept as-is so existing artefacts stay readable; read them as "a
        # fingerprint from this Julia", not as a checksum anyone can recompute.
        "config_hash" => bytes2hex(reinterpret(UInt8, [hash(read(cfg["_source"], String))])),
        "seed" => seed, "m_prior" => collect(info.m_prior),
        "stage0_written" => string(now_str()),
        # The IDENTITY of the truth, not its contents. figures/figstyle.jl
        # compares this against the hash of data/truth_murray<column>.jld2,
        # which catches a config declaring column = "3P" while pointing at
        # another column's truth -- a reviewer sealed exactly that and drew two
        # figures from it. Recording the PATH instead was the first attempt and
        # this file's own key whitelist refused it, correctly: a path is
        # provenance the firewall has not authorised, and the hash already
        # carries what is needed.
        "truth_file_sha" => filesha(tpath),
    )

    # THE SEAL MUST MATCH THE CONFIG THAT ASKED FOR IT. A reviewer showed that
    # this stage records the seed VARIABLE it used, so a Stage 0 that alters the
    # seed records the altered value and test_stage0_roundtrip -- which rebuilds
    # using the artefact's OWN recorded seed -- passes. `seed = Int(d["seed"]) + 1`
    # sealed data drawn at 124 for a config asking 123, printed "seed 124", and
    # the firewall test said OK.
    #
    # Check against the RESOLVED config, not a fresh TOML.parsefile of the
    # source. The first version of this check re-parsed the file on disk to
    # avoid trusting the local variable, and in doing so bypassed `extends`
    # inheritance: all seven committed configs inherit `seed` from
    # ../config.toml and carry none of their own, so Stage 0 ran the forward
    # model and the noise draw and then died with KeyError: key "seed" not
    # found -- on EVERY committed config. It killed every real run for a round,
    # and neither CI nor the suite runs Stage 0, which is how it shipped.
    #
    # The resolved config is still not the local variable, so this catches what
    # it was built to catch: a Stage 0 that reads the seed and then alters it
    # before drawing. test/test_stage0_seed.jl runs this on a COMMITTED config,
    # which is the gap that let the first version through -- the prototype that
    # verified it wrote `seed` into a config it generated.
    let want = Int(cfg["data"]["seed"])
        seed == want || error("""
    Stage 0 sealed seed $seed but $(cfg["_source"]) asks for $want.

      The artefact records the seed it used, so nothing downstream would notice.
    """)
        chk = build_aem_jac(T, want).aem
        (vec(chk.dlow) == vec(aem.dlow) && vec(chk.dhigh) == vec(aem.dhigh)) || error("""
    Stage 0's sealed data are not the draw at seed $want.

      An independent rebuild at the config's seed differs from what was sealed,
      so the noise realisation every later stage inverts is not the one the
      config specifies.
    """)
    end

    # ---- the whitelist is ASSERTED, not assumed -----------------------
    extra = setdiff(keys(A), SEALED_KEYS)
    isempty(extra) || error("Stage 0 tried to seal keys that are not whitelisted: " *
                            join(sort(collect(extra)), ", "))
    # and nothing that smells of truth got in by name
    banned = filter(k -> occursin(r"(?i)truth|sigma_bulk|σ_bulk|clay|theta_total|θ_|phi|φ|sigma_w|σw|σ_s|units|basin|archie",
                                 k) && k != "truth_file_sha", collect(keys(A)))
    isempty(banned) || error("Stage 0 artefact contains truth-like keys: " * join(banned, ", "))

    odir = joinpath(ROOT, out["dir"]); mkpath(odir)
    opath = joinpath(odir, "stage0_$(cfg["name"]).jld2")
    jldopen(opath, "w") do f
        for (k, v) in A; f[k] = v; end
    end

    println("\n  SEALED. Stages 1-6 may read this and the config, nothing else.")
    @printf("    %d keys, all whitelisted; %d truth-like keys rejected: none\n", length(A), 0)
    @printf("    truth identity recorded as a hash (%s...), never its values\n",
            A["truth_file_sha"][1:12])
    println("  wrote ", relpath(opath, ROOT))
end

now_str() = Libc.strftime("%Y-%m-%dT%H:%M:%S", time())
function filesha(p)
    # identity of the truth file, so a verdict can be traced to the column it
    # came from without the column's values being readable downstream
    bytes2hex(reinterpret(UInt8, [hash(read(p))]))
end

main()
