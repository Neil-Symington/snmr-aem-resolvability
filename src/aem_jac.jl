# aem_jac.jl — a SkyTEM operator with calcjacobian = true, carrying the SAME
# noisy data the trans-D chain saw. Included at Main scope AFTER
# gadi/src/aem_setup.jl, whose geometry and noise constants it uses.
#
# The shared build_skytem_operator does not set calcjacobian, so the
# Gauss-Newton steps in occam_core.jl cannot use it.
#
# makenoisydata! reseeds from its OWN rseed kwarg, so it must be passed
# explicitly. 123 is its default and therefore the realisation every earlier run
# in this project used — including the chains — so seed = 123 reproduces the main
# dataset and makes posterior chi2 values directly comparable with these. That
# equality is verified: chi2 of the truth model agrees with the value stored by
# aem_equivalence.jl to 12 significant figures.
function build_aem_jac(truth::AbstractDict, seed::Integer)
    grid = build_aem_grid()
    σt = piecewise_interp(grid.z, truth["z"], truth["σ_bulk"]); σt[grid.z .< 0] .= 1e-12
    ρt = 1.0 ./ σt
    SK = HiQGA.transD_GP.SkyTEM1DInversion
    aem = SK.dBzdt(; timeslow=LM_times, ramplow=LM_ramp, zRxlow=-42.0, zTxlow=-40.0,
        timeshigh=HM_times, ramphigh=HM_ramp, zRxhigh=-42.0, zTxhigh=-40.0,
        rRx=13.0, rTx=rTx, z=grid.z, ρ=copy(ρt), lowpassfcs=lowpassfcs, calcjacobian=true)
    SK.makenoisydata!(aem, log10.(ρt[2:end]); rseed=Int(seed), σ_halt_low=LM_noise,
                      σ_halt_high=HM_noise, units=1e-12, showplot=false)
    (aem=aem, grid=grid, σt=σt, T=truth)
end
