using Random
using PyPlot
using OrdinaryDiffEq
using LinearAlgebra

#include("../SuperradiantLasing.jl")
#using .SuperradiantLasing

obs_exc(S, M, model) = M + model.N / 2

N = 200
N_traj = 2000

model = PopulationModel(
    N,
    N ÷ 2,
    1.0,   # global_decay
    0.0,   # global_pump
    0.0,   # local_decay
    0.3,   # local_pump
    0.0,   # local_dephasing
)

t_grid, avg_exc = simulate_population_ensemble(
    model,
    obs_exc,
    N_traj;
    t_max = 5.0,
    S0 = model.Jmax,
    M0 = -model.Jmax,
    n_grid = 1000,
    rng = Xoshiro(1234),
)

params = Dict(
    "global_decay" => model.global_decay,
    "global_pump" => model.global_pump,
    "local_pump" => model.local_pump,
    "local_decay" => model.local_decay,
    "local_dephasing" => model.local_dephasing,
)

basis, index_mapping, R, nstates = set_up_rate_matrix(model.Jmax, params)
u0 = initial_state(index_mapping, (model.Jmax, -model.Jmax), nstates)

f!(du, u, p, t) = mul!(du, R, u)
prob = ODEProblem(f!, u0, (0.0, maximum(t_grid)))
sol = solve(prob, Tsit5(); saveat = t_grid, reltol = 1e-9, abstol = 1e-9)

obs_vec_exc = [obs_exc(S, M, model) for (S, M) in basis]
ode_exc = [dot(obs_vec_exc, u) for u in sol.u]

figure()
plot(t_grid, avg_exc, label = "MC")
plot(t_grid, ode_exc, "--", label = "ODE")
xlabel("t")
ylabel(L"\langle n_{\mathrm{exc}} \rangle")
legend()
tight_layout()