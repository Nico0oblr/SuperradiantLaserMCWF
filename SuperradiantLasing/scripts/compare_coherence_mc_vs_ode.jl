using Random
using PyPlot
using OrdinaryDiffEq
using LinearAlgebra

#include("../src/SuperradiantLasing.jl")
#using .SuperradiantLasing

N = 200
N_traj = 10000
t_ss = 80.0
τ_max = 10.0
n_grid = 1200
rng = Xoshiro(1234)

pop_model = PopulationModel(
    N,
    N ÷ 2,
    1.0,
    0.0,
    0.0,
    0.3,
    0.0,
)

S_samples, M_samples = get_stationary_population_samples(
    pop_model,
    N_traj;
    t_ss = t_ss,
    S0 = pop_model.Jmax,
    M0 = -pop_model.Jmax,
    rng = rng,
)

println("Mean stationary sample S = ", sum(S_samples) / length(S_samples))
println("Mean stationary sample M = ", sum(M_samples) / length(M_samples))

coh_model = CoherenceModel(pop_model)

τ_grid, Cτ = simulate_coherence_correlator_from_samples(
    S_samples,
    M_samples,
    coh_model;
    τ_max = τ_max,
    n_grid = n_grid,
    rng = rng,
)

params = Dict(
    "global_decay" => pop_model.global_decay,
    "global_pump" => pop_model.global_pump,
    "local_pump" => pop_model.local_pump,
    "local_decay" => pop_model.local_decay,
    "local_dephasing" => pop_model.local_dephasing,
)

basis_pop, idx_pop, Rpop, Npop = set_up_rate_matrix(pop_model.Jmax, params)
pss = steady_state_from_rate_matrix(Rpop)

basis_coh, idx_coh, Rcoh, Ncoh = set_up_coherence_rate_matrix(pop_model.Jmax, params)
c0 = coherence_initial_from_population_ss(basis_pop, pss, idx_coh)
rvec = coherence_readout_vector(basis_coh)

fcoh!(du, u, p, t) = mul!(du, Rcoh, u)
prob = ODEProblem(fcoh!, c0, (0.0, τ_max))
sol = solve(prob, Tsit5(); saveat = τ_grid, reltol = 1e-9, abstol = 1e-9)

C_det = [dot(rvec, u) for u in sol.u]

figure()
plot(τ_grid, Cτ, label = "MC")
plot(sol.t, C_det, "--", label = "ODE")
xlabel("τ")
ylabel("C(τ)")
legend()
tight_layout()

ω_grid_mc, Sω_mc = spectrum_even_fft(τ_grid, Cτ)
ω_grid_det, Sω_det = spectrum_even_fft(sol.t, C_det)

figure()
plot(ω_grid_mc, Sω_mc, label = "MC")
plot(ω_grid_det, Sω_det, "--", label = "ODE")
xlabel(L"\omega")
ylabel(L"S(\omega)")
xlim(-15.0, 15.0)
legend()
tight_layout()