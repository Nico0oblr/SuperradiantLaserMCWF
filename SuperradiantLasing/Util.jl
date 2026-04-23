using Printf

# --------------------------------
# Scan gamma, compute spectra, extract linewidths
# --------------------------------
function scan_gamma_linewidth(
    gammas;
    N = 1000,
    N_traj = 4000,
    τ_max = 10.0,
    t_ss = 20.0,
    n_grid = 2000,
    ω_max = 10.0,
    n_ω = 4000,
    pump = 0.3,
    Γ = 1.0,
    deph = 0.0,
    seed = 1234,
    plot_spectra = true,
)
    widths = Float64[]
    spectra = Vector{Vector{Float64}}()

    if plot_spectra
        figure()
    end

    for gamma in gammas
        pop_model = PopulationModel(
            N,
            N ÷ 2,
            Γ,       # global_decay
            0.0,     # global_pump
            gamma,   # local_decay
            pump,    # local_pump
            deph,    # local_dephasing
        )

        τ_grid, Cτ, _, _ = simulate_correlator_ergodic(
            pop_model;
            N_traj = N_traj,
            t_ss = t_ss,
            τ_max = τ_max,
            n_grid = n_grid,
            S0 = pop_model.Jmax,
            M0 = -pop_model.Jmax,
            rng = Xoshiro(seed),
        )

        ω_grid, Sω = spectrum_even_fft(τ_grid, Cτ)

        push!(spectra, Sω)
        push!(widths, spectrum_fwhm(ω_grid, Sω))

        if plot_spectra
            plot(ω_grid, Sω ./ maximum(Sω), label = L"\gamma=%$gamma")
        end
    end

    xlim(left = -200.0, right = 200.0)

    if plot_spectra
        xlabel(L"\omega-\omega_0")
        ylabel("normalized spectrum")
        title("Spectrum scan vs local decay")
        legend()
        display(gcf())
        close("all")
    end

    fit = fit_linewidth_vs_gamma(collect(gammas), widths)
    a, b = coef(fit)

    figure()
    scatter(gammas, widths, label="data")
    xlabel(L"\gamma")
    ylabel("FWHM")
    title("Linewidth vs local decay")
    legend()
    display(gcf())
    close("all")

    return (
        gammas = collect(gammas),
        ω_grid = ω_grid,
        spectra = spectra,
        widths = widths,
        fit = fit,
        fit_params = coef(fit),
    )
end