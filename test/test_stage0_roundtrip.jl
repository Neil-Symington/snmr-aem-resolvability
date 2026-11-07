# test_stage0_roundtrip.jl -- the firewall must not change the inverse problem.
#
# Stage 4 rebuilds the AEM operator from the SEALED Stage 0 artefact instead of
# from the truth column. That is only legitimate if the two operators pose the
# SAME inverse problem. This asserts it: for arbitrary models, chi2 computed by
# the truth-built operator and by the artefact-built operator must agree to
# machine precision.
#
# This is the test that would have caught a stale W. W is Diagonal(1/sigma) and
# is not rebuilt by assigning sigma, so a reconstruction that forgets it weights
# the new data by the old noise and every chi2 downstream is wrong by a factor
# nothing else would reveal.
using JLD2, Statistics, Printf, LinearAlgebra, SparseArrays
const HERE = @__DIR__; const ROOT = normpath(joinpath(HERE, ".."))
include(joinpath(ROOT, "src", "aem_setup.jl"))
include(joinpath(ROOT, "src", "aem_jac.jl"))
const SK = HiQGA.transD_GP.SkyTEM1DInversion
chi2(aem, m) = (SK.getresidual(aem, m); norm(aem.W * aem.res)^2)

function rebuild(A)
    ρ0 = fill(100.0, length(A["z"])); ρ0[A["z"] .< 0] .= 1e12
    aem = SK.dBzdt(; timeslow=A["times_low"], ramplow=A["ramp_low"],
        zRxlow=A["z_rx"], zTxlow=A["z_tx"], timeshigh=A["times_high"],
        ramphigh=A["ramp_high"], zRxhigh=A["z_rx"], zTxhigh=A["z_tx"],
        rRx=A["r_rx"], rTx=A["r_tx"], z=copy(A["z"]), ρ=ρ0,
        lowpassfcs=A["lowpass_fc"], calcjacobian=true)
    aem.dlow = copy(A["d_low"]); aem.dhigh = copy(A["d_high"])
    aem.σlow = copy(A["sigma_low"]); aem.σhigh = copy(A["sigma_high"])
    aem.selectlow = copy(A["select_low"]); aem.selecthigh = copy(A["select_high"])
    aem.ndatalow = Int(A["ndata_low"]); aem.ndatahigh = Int(A["ndata_high"])
    sg = vcat(aem.σlow[aem.selectlow], aem.σhigh[aem.selecthigh])
    aem.W = sparse(Diagonal(1 ./ sg)); aem.res = zeros(length(sg))
    aem.J = zeros(length(sg), length(aem.ρ) - 1)
    aem
end

function main()
fails = 0
COLS = ("3P", "2P", "3S", "2S")   # local: const is illegal inside main()
skipped = String[]
for col in COLS
    tf = joinpath(HERE, "..", "data", "truth_murray$col.jld2")
    af = joinpath(ROOT, "outputs", "synth_$col", "stage0_synth_$col.jld2")
    (isfile(tf) && isfile(af)) || (println("SKIP $col"); push!(skipped, col); continue)
    T = load(tf); A = load(af)
    a1 = build_aem_jac(T, Int(A["seed"])).aem      # built from truth
    a2 = rebuild(A)                                # built from the artefact
    n = length(a1.ρ) - 1
    worst = 0.0
    for m in (fill(-2.0, n), fill(-1.0, n), fill(-2.5, n),
              [-2.0 + 0.5sin(i/5) for i in 1:n])
        c1, c2 = chi2(a1, m), chi2(a2, m)
        worst = max(worst, abs(c1 - c2) / max(c1, 1.0))
    end
    ok = worst < 1e-12
    @printf("  %-4s chi2 agreement over 4 models: worst rel diff %.3e  %s\n",
            col, worst, ok ? "PASS" : "FAIL")
    ok || (fails += 1)
    # and the data itself must be bit-identical
    bits = a1.dlow == a2.dlow && a1.dhigh == a2.dhigh &&
           a1.σlow == a2.σlow && a1.σhigh == a2.σhigh
    @printf("       data and noise bit-identical: %s\n", bits ? "PASS" : "FAIL")
    bits || (fails += 1)
end
println()
# "FIREWALL ROUND-TRIP OK" used to print after SKIPping all four columns, which
# a reviewer caught: CI never has artefacts, so the line always appeared and
# always meant nothing. Say SKIPPED when nothing was compared.
if length(skipped) == length(COLS)
    println("FIREWALL ROUND-TRIP SKIPPED -- no artefacts for any column; nothing was compared")
else
println(fails == 0 ? "FIREWALL ROUND-TRIP OK -- the sealed artefact poses the same inverse problem" :
                     "$fails FAILURE(S)")
end
fails == 0 || exit(1)
end
main()
