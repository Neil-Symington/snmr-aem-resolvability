# salinity_ceiling.jl -- the clay-free salinity ceiling from a band's upper edge.
# No solves; arithmetic on stored numbers.
#
# The gate's upper edge converts to an UPPER BOUND on salinity by assuming no
# clay at all: every bit of the window's conductance is attributed to pore
# water. stages/stage5_quality.jl:17-18 is the relation,
#
#     sigma_bulk = sigma_w * theta_total^m + clay_term
#     clay-free:  sigma_w = sigma_bulk / theta_total^m
#     TDS        = sigma_w * mg_per_sm
#
# so the bound is largest where theta_total^m is smallest. theta_total comes
# from [priors].porosity_range and m from [priors].m_lithology; both are stated
# with every number below, because the bound is a function of them, not of the
# data alone.
using Printf
const MGL = 6400.0
const PHI_LO, PHI_HI = 0.25, 0.40      # [priors].porosity_range
const M_LO, M_HI = 1.5, 1.9            # clayey_sand
const L = 16.0                         # window thickness, m

tds(C, θ, m) = (C / L) / θ^m * MGL

"The bound at the conservative corner, the central point, and the tight corner."
function report(label, C)
    @printf("%-34s C = %8.4f S  sigma_bulk = %.5f S/m\n", label, C, C/L)
    @printf("   %-30s %10s %6s %12s\n", "theta_total, m", "theta^m", "", "clay-free TDS")
    for (nm, θ, m) in (("widest bound  (0.25, 1.9)", PHI_LO, M_HI),
                       ("central       (0.325, 1.7)", (PHI_LO+PHI_HI)/2, (M_LO+M_HI)/2),
                       ("tightest bound (0.40, 1.5)", PHI_HI, M_LO))
        @printf("   %-30s %10.5f %6s %9.0f mg/L\n", nm, θ^m, "", tds(C, θ, m))
    end
    println()
end

println("3P salinity ceiling, clay-free, from the band's UPPER edge")
println("window 7.5-23.5 m (16.0 m); theta_total ~ U(0.25, 0.40); m ~ U(1.5, 1.9); mg_per_sm = 6400")
println("="^78)
report("OLD gate, reported upper edge",        0.3447)
report("NEW absolute edge, lower bracket",     0.3582)
report("NEW absolute edge, upper bracket",     0.4509)
report("NEW relative edge (chiref 25.089)",    0.2845)
# ── the SNMR lower limit on porosity (M1.6) ───────────────────────────────
# theta_total >= lowest p05 of window water content across the ensembles. On 3P
# that is 0.265, from the 22 Oct SNMR vintage, which runs 25% low on the one
# column where both vintages exist -- so this row is provisional. A lower-biased
# measurement makes a WEAKER lower limit, never a wrong one.
const PHI_SNMR = 0.265
println("With the SNMR lower limit applied: theta_total in [", PHI_SNMR, ", ", PHI_HI, "]")
println("="^78)
function report_snmr(label, C)
    @printf("%-34s C = %8.4f S\n", label, C)
    for (nm, θ, m) in (("widest bound  (0.265, 1.9)", PHI_SNMR, M_HI),
                       ("central       (0.333, 1.7)", (PHI_SNMR+PHI_HI)/2, (M_LO+M_HI)/2),
                       ("tightest bound (0.40, 1.5)", PHI_HI, M_LO))
        @printf("   %-30s %10.5f %6s %9.0f mg/L\n", nm, θ^m, "", tds(C, θ, m))
    end
    @printf("   factor across the box: %.2f  (was %.2f without the SNMR limit)\n\n",
            tds(C,PHI_SNMR,M_HI)/tds(C,PHI_HI,M_LO), tds(C,PHI_LO,M_HI)/tds(C,PHI_HI,M_LO))
end
report_snmr("NEW absolute edge, lower bracket", 0.3582)
report_snmr("NEW relative edge (chiref 25.089)", 0.2845)

println("Truth, for comparison:")
@printf("   true window conductance 0.2063 S -> clay-free TDS at the central point %.0f mg/L\n",
        tds(0.2063, 0.325, 1.7))
@printf("   TRUE pore-water TDS is 187 mg/L; the clay-free bound is an UPPER limit, so it\n")
@printf("   should and does sit above that -- the gap is what the clay is carrying.\n")
@printf("   CLAY SHARE AT TRUTH: clay carries %.1f%% of sigma_bulk in this window\n", 67.0)
@printf("   (0.00864 of 0.01290 S/m). On 2P it is 73.8%% (0.03162 of 0.04285 S/m).\n")
@printf("   So two thirds of what the AEM measures here is clay, not water, which is\n")
@printf("   why the clay-free bound sits a factor of 3 above the true pore water.\n")
