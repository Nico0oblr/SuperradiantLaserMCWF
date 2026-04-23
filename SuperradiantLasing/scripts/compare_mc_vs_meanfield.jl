using Random
using PyPlot

#include("../SuperradiantLasing.jl")
#using .SuperradiantLasing

N = 200
mc_model = PopulationModel(
    N,
    N ÷ 2,
    1.0,
    0.0,
    0.0,
    0.3,
    0.0,
)

out = compare_mc_to_meanfield(
    mc_model;
    N_traj = 2000,
    t_max = 1000.0,
    n_grid = 1000,
    S0 = mc_model.Jmax,
    M0 = -mc_model.Jmax,
    g = 0.01,
    rng = Xoshiro(1234),
)

println("Matched parameters:")
println("  thresholds γp⁻, γp⁺ = ", out.thresholds)
println("  lasing regime? ", out.lasing)
println("  max |Δexc| = ", out.errors[1])
println("  max |ΔM|   = ", out.errors[2])

figure()
plot(out.t, out.mc_M, label = "MC")
plot(out.t, out.mf_sz, "--", label = "MF")
xlabel("t")
ylabel("inversion")
legend()
tight_layout()