function compatible_mf_params(mc_model::PopulationModel; g::Float64 = 1.0, κ::Union{Nothing,Float64} = nothing)
    Γeff = mc_model.global_decay
    κ === nothing && (κ = g^2 / Γeff)
    γ = mc_model.local_decay
    γp = mc_model.local_pump
    return (N = mc_model.N, g = g, κ = κ, γ = γ, γp = γp, Γeff = Γeff)
end

function compare_mc_to_meanfield(
    mc_model::PopulationModel;
    N_traj::Int = 2000,
    t_max::Float64 = 5.0,
    n_grid::Int = 1000,
    S0::Int = mc_model.Jmax,
    M0::Int = -mc_model.Jmax,
    g::Float64 = 1.0,
    κ::Union{Nothing,Float64} = nothing,
    rng::AbstractRNG = Random.default_rng(),
    make_plots::Bool = true,
)
    pars = compatible_mf_params(mc_model; g = g, κ = κ)
    N, g, κ, γ, γp = pars.N, pars.g, pars.κ, pars.γ, pars.γp

    obs_exc(S, M, model) = M + model.N / 2
    obs_M(S, M, model) = M

    t_grid, avg_M_mc = simulate_population_ensemble(
        mc_model,
        obs_M,
        N_traj;
        S0 = S0,
        M0 = M0,
        t_max = t_max,
        n_grid = n_grid,
        rng = rng,
    )

    sol = solve_mf_superradiant_laser(
        N, g, κ, γ, γp;
        α0 = 0.0 + 0im,
        s0 = 1im / N,
        sz0 = M0,
        tspan = (0.0, t_max),
        saveat = t_grid,
    )

    sz_mf = inversion.(sol.u)
    exc_mf = sz_mf .+ N / 2

    γp_minus, γp_plus = pump_thresholds(N, g, κ, γ)
    lasing = has_lasing_steady_state(N, g, κ, γ, γp)

    if make_plots
        figure()
        plot(t_grid, avg_M_mc, label="MC")
        plot(t_grid, sz_mf, "--", label="MF")
        xlabel("t")
        ylabel("inversion")
        legend()
        display(gcf())
        close("all")
    end

    println("Matched parameters:")
    println("  thresholds γp⁻, γp⁺ = ", (γp_minus, γp_plus))
    println("  lasing regime? ", lasing)
    println(g * sqrt(N))

    return (
        t = t_grid,
        mc_M = avg_M_mc,
        mf_sz = sz_mf,
        lasing = lasing,
        thresholds = (γp_minus, γp_plus),
        params = pars,
    )
end