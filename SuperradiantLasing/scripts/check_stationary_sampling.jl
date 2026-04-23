using Random
using PyPlot

#include("../src/SuperradiantLasing.jl")
#using .SuperradiantLasing

function integer_hist(x::Vector{Int})
    xmin = minimum(x)
    xmax = maximum(x)
    counts = zeros(Int, xmax - xmin + 1)
    @inbounds for xi in x
        counts[xi - xmin + 1] += 1
    end
    xs = collect(xmin:xmax)
    return xs, counts
end

function plot_integer_step(xs, p; label = "")
    plot(xs, p; drawstyle = "steps-mid", linewidth = 1.5, label = label)
end

N = 2000
N_traj = 10000
t_ss = 20.0

pop_model = PopulationModel(
    N,
    N ÷ 2,
    1.0,
    0.0,
    0.0,
    0.3,
    0.0,
)

Ss_erg, Ms_erg = get_stationary_population_samples_ergodic(
    pop_model,
    N_traj;
    t_ss = t_ss,
    t_sample = 20.0,
    S0 = pop_model.Jmax,
    M0 = pop_model.Jmax,
    rng = Xoshiro(1234),
)

Ss, Ms = get_stationary_population_samples(
    pop_model,
    N_traj;
    t_ss = t_ss,
    S0 = pop_model.Jmax,
    M0 = pop_model.Jmax,
    rng = Xoshiro(1234),
)

xs1, c1 = integer_hist(Ss)
xs2, c2 = integer_hist(Ss_erg)
p1 = c1 ./ sum(c1)
p2 = c2 ./ sum(c2)

xm1, cm1 = integer_hist(Ms)
xm2, cm2 = integer_hist(Ms_erg)
pm1 = cm1 ./ sum(cm1)
pm2 = cm2 ./ sum(cm2)

figure()
plot_integer_step(xs1, p1, label = "independent")
plot_integer_step(xs2, p2, label = "ergodic")
xlabel("S")
ylabel("P(S)")
legend()
tight_layout()

figure()
plot_integer_step(xm1, pm1, label = "independent")
plot_integer_step(xm2, pm2, label = "ergodic")
xlabel("M")
ylabel("P(M)")
legend()
tight_layout()
