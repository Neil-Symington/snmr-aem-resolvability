# both_rules.jl -- M1.1 item 5. The accepted band under BOTH rules, with truth
# inside or outside for each. Reads stored Stage 4 output; runs no solves.
#
#   absolute:  chi2 <= N + 3.84        (current; drives the verdict)
#   relative:  chi2 <= chiref + 3.84   (profile likelihood, for comparison)
#
# LINEAGE. The stored 2P and 2S gates were run on the RETIRED 2P column
# (HYP_SYNTH/truth_murray2P.jld2), so they are paired here with that column's
# truth, not the release one. 3P and 3S pair with the column still in use.
using JLD2, Printf

const P    = "/home/symo/Documents/GA/AEM_SMR_Resolvability_Quickrun/outputs"
const HYP  = "/home/symo/Documents/GA/AEM_SMR_Resolvability_Quickrun/synthetics/outputs/HYP_SYNTH"
const HERE = normpath(@__DIR__); const ROOT = normpath(joinpath(HERE,".."))
mkpath(joinpath(HERE,"results"))
filesha(p) = bytes2hex(reinterpret(UInt8, [hash(read(p))]))

# column => (window, truth file used BY THAT STORED GATE RUN)
CASES = [("3P",(7.5,23.5),  joinpath(HYP,"truth_murray3P.jld2")),
         ("3S",(7.5,23.5),  joinpath(HYP,"truth_murray3S.jld2")),
         ("2P",(17.5,33.5), joinpath(HYP,"truth_murray2P.jld2")),
         ("2S",(17.5,33.5), joinpath(HYP,"truth_murray2S.jld2"))]

function nodes(D)
    if haskey(D,"coarse")
        C=D["coarse"]; (C["ltc"], C["chi2"], C["good"])
    else
        lam=string(D["lambda_hat"]); S=D["sweeps"]
        k = haskey(S,lam) ? lam : first(keys(S)); sw=S[k]
        (D["lt"], sw.chi2, sw.good)
    end
end

"Widest contiguous span of nodes that are solver-good and below `thr`."
function band(lt, chi2, good, thr, L)
    ok = [good[i] && chi2[i] <= thr for i in eachindex(lt)]
    any(ok) || return (NaN, NaN)
    i = findfirst(ok); j = findlast(ok)
    (10^lt[i]*L, 10^lt[j]*L)
end

io=open(joinpath(HERE,"results","both_rules.md"),"w")
emit(s="") = (println(s); println(io,s))
emit("# M1.1 item 5 — accepted band under both rules")
emit()
emit("`absolute` is the current rule and drives the verdict. `relative` is shown for")
emit("comparison only. Bands are the widest contiguous run of nodes that are both")
emit("solver-good and below the threshold; no new solves.")
emit()
emit("| column | truth hash | N | chiref | true window C (S) | absolute band (S) | truth | relative band (S) | truth |")
emit("|---|---|---|---|---|---|---|---|---|")
for (tag,w,tf) in CASES
    f=joinpath(P,"synth_$tag","stage4_synth_$tag.jld2")
    (isfile(f) && isfile(tf)) || (emit("| $tag | missing | | | | | | | |"); continue)
    D=load(f); lt,chi2,good = nodes(D); L=D["thick"]; N=D["ndata"]; cr=D["chiref"]
    T=load(tf); z=T["z"]; dz=T["dz"]; sb=T["σ_bulk"]
    tw=sum(sb[i]*dz for i in eachindex(z) if w[1]<=z[i]<w[2])
    aa=band(lt,chi2,good,N+3.84,L); rr=band(lt,chi2,good,cr+3.84,L)
    ins(b)= (isnan(b[1]) ? "none" : (tw>=b[1] && tw<=b[2] ? "**in**" : "**OUT**"))
    emit(@sprintf("| %s | `%s` | %d | %.2f | %.4f | %.4f - %.4f | %s | %.4f - %.4f | %s |",
         tag, filesha(tf)[1:12], N, cr, tw, aa[1], aa[2], ins(aa), rr[1], rr[2], ins(rr)))
end
emit()
emit("Windows: test one 7.5-23.5 m, test two 17.5-33.5 m.")
close(io); println("\nwrote ", joinpath(HERE,"results","both_rules.md"))
