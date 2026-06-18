# ============================================
# MixingApproximations.jl
#
# Quasi-stationary tail approximations for coherence-sector
# correlators. This file assumes the following are already defined:
#
#   CoherenceModel
#   CoherenceState
#   CoherenceEventBuffer
#   valid_coherence_state
#   coherence_step!
#   build_coherence_events!
#   initial_coherence_weight
#   coherence_readout
#   CoherenceModel(::PopulationModel)
#   get_stationary_population_samples
#   get_stationary_population_samples_ergodic
#
# Convention:
#   coherence label (S,M) represents |S,M-1><S,M|.
#
# LHS correlator:
#   C_LHS(t) = <S†(t) S(0)>
#   source   = S rho_ss
#   label    = (S_sample, M_sample)
#
# RHS correlator:
#   C_RHS(t) = <S(0) S†(t)>
#   source   = rho_ss S
#   label    = (S_sample, M_sample + 1)
#
# Susceptibility commutator integrand:
#   C_comm(t) = C_LHS(t) - C_RHS(t)
# ============================================

# --------------------------------------------
# Initial coherence-sector label from a sampled population state
# --------------------------------------------

@inline function initial_coherence_label(side::Symbol, S::Int, M::Int)
    if side === :LHS
        return S, M
    elseif side === :RHS
        return S, M + 1
    else
        error("Unknown side = $side. Use :LHS or :RHS.")
    end
end


# --------------------------------------------
# Exact burn-in up to t_mix, accumulated directly into Cτ
#
# Returns:
#   survived     : whether trajectory is still alive at t_mix
#   amp_mix      : alive amplitude at t_mix
#   sink_mix     : sink rate at t_mix
#   S_mix, M_mix : state at t_mix
# --------------------------------------------

function accumulate_coherence_sector_to_mix!(
    Cτ::AbstractVector{Float64},
    model::CoherenceModel,
    S0::Int,
    M0::Int,
    τ_grid::AbstractVector,
    t_mix::Float64;
    rng::AbstractRNG = Random.default_rng(),
)
    st = CoherenceState(S0, M0, 0.0, false)
    buf = CoherenceEventBuffer()

    w_in = initial_coherence_weight(S0, M0)

    if !valid_coherence_state(S0, M0, model.Jmax)
        return false, 0.0, 0.0, S0, M0
    end

    idx = 1
    nτ = length(τ_grid)

    while !st.absorbed && st.t < t_mix
        current_val = w_in * coherence_readout(st.S, st.M)

        build_coherence_events!(buf, model, st)

        r_internal = sum_coherence_internal_rates(buf)
        rtot = r_internal + buf.sink_rate

        if rtot <= 0.0
            @inbounds while idx <= nτ && τ_grid[idx] <= t_mix
                Cτ[idx] += current_val
                idx += 1
            end
            return false, 0.0, 0.0, st.S, st.M
        end

        dt = -log(rand(rng)) / rtot
        t_next = st.t + dt

        fill_until_t = min(t_next, t_mix)

        @inbounds while idx <= nτ && τ_grid[idx] <= fill_until_t
            Cτ[idx] += current_val
            idx += 1
        end

        if t_next >= t_mix
            # The jump lies after t_mix, so the state at t_mix is still the current state.
            amp_mix = current_val
            sink_mix = buf.sink_rate
            return true, amp_mix, sink_mix, st.S, st.M
        end

        # The jump occurs before t_mix, so now actually apply it.
        st.t = t_next

        x = rand(rng) * rtot
        acc = 0.0

        jumped_internal = false
        @inbounds for k in 1:buf.n
            acc += buf.rates[k]
            if x <= acc
                st.S += buf.dS[k]
                st.M += buf.dM[k]
                jumped_internal = true
                break
            end
        end

        if !jumped_internal
            st.absorbed = true
            return false, 0.0, 0.0, st.S, st.M
        end
    end

    return false, 0.0, 0.0, st.S, st.M
end


# --------------------------------------------
# Sampled quasi-stationary exponential tail
# --------------------------------------------

@inline function add_sampled_qs_tail!(
    Cτ::AbstractVector,
    τ_grid::AbstractVector,
    i_mix::Int,
    t_mix::Float64,
    amp_mix::Float64,
    λ_qs::Float64,
    rng::AbstractRNG,
)
    if amp_mix == 0.0
        return nothing
    end

    if λ_qs <= 0.0
        @inbounds for i in i_mix+1:length(τ_grid)
            Cτ[i] += amp_mix
        end
        return nothing
    end

    T_tail = -log(rand(rng)) / λ_qs
    t_abs = t_mix + T_tail

    @inbounds for i in i_mix+1:length(τ_grid)
        if τ_grid[i] <= t_abs
            Cτ[i] += amp_mix
        else
            break
        end
    end

    return nothing
end


# --------------------------------------------
# Main QS approximation from stationary population samples
#
# side = :LHS gives <S†(t) S(0)>
# side = :RHS gives <S(0) S†(t)>
#
# tail_mode = :deterministic
#   Adds amp_mix * exp(-λ_qs * (τ - t_mix)) globally.
#
# tail_mode = :sampled
#   Adds amp_mix * 1_{τ < t_mix + Exp(λ_qs)} trajectory by trajectory.
# --------------------------------------------

function simulate_coherence_correlator_from_samples_qs(
    S_samples::AbstractVector{Int},
    M_samples::AbstractVector{Int},
    model::CoherenceModel;
    side::Symbol,
    τ_max::Float64,
    t_mix::Float64,
    n_grid::Int = 1000,
    tail_mode::Symbol = :deterministic,
    rng::AbstractRNG = Random.default_rng(),
)
    @assert side in (:LHS, :RHS)
    @assert length(S_samples) == length(M_samples)
    @assert 0.0 <= t_mix <= τ_max

    τ_grid = collect(range(0.0, τ_max, length = n_grid))
    Cτ = zeros(Float64, n_grid)

    N_samples = length(S_samples)

    survived = falses(N_samples)
    amp_mix = zeros(Float64, N_samples)
    sink_mix = zeros(Float64, N_samples)
    S_mix = zeros(Int, N_samples)
    M_mix = zeros(Int, N_samples)

    # Exact evolution to physical time t_mix.
    # The sampling grid is only used for recording Cτ, not for deciding when mixing stops.
    for n in eachindex(S_samples)
        S0, M0 = initial_coherence_label(side, S_samples[n], M_samples[n])

        survived_n, amp_n, sink_n, S_n, M_n = accumulate_coherence_sector_to_mix!(
            Cτ,
            model,
            S0,
            M0,
            τ_grid,
            t_mix;
            rng = rng,
        )

        survived[n] = survived_n
        amp_mix[n] = amp_n
        sink_mix[n] = sink_n
        S_mix[n] = S_n
        M_mix[n] = M_n
    end

    n_surv = count(survived)

    λ_qs = n_surv == 0 ? 0.0 : sum(sink_mix[survived]) / n_surv

    # This index is now only a grid bookkeeping object for the tail.
    # It does not control physical mixing.
    i_tail = searchsortedlast(τ_grid, t_mix)

    if tail_mode === :deterministic
        amp_total = sum(amp_mix[survived])

        if amp_total != 0.0
            if λ_qs > 0.0
                @inbounds for i in i_tail+1:n_grid
                    Cτ[i] += amp_total * exp(-λ_qs * (τ_grid[i] - t_mix))
                end
            else
                @inbounds for i in i_tail+1:n_grid
                    Cτ[i] += amp_total
                end
            end
        end

    elseif tail_mode === :sampled
        for n in eachindex(S_samples)
            if survived[n]
                add_sampled_qs_tail!(
                    Cτ,
                    τ_grid,
                    i_tail,
                    t_mix,
                    amp_mix[n],
                    λ_qs,
                    rng,
                )
            end
        end

    else
        error("Unknown tail_mode = $tail_mode. Use :deterministic or :sampled.")
    end

    Cτ ./= N_samples

    stats = (
        side = side,
        λ_qs = λ_qs,
        n_survived = n_surv,
        survival_fraction = n_surv / N_samples,
        mean_amp_mix = n_surv == 0 ? 0.0 : sum(amp_mix[survived]) / n_surv,
        mean_sink_mix = λ_qs,
        S_mix = S_mix,
        M_mix = M_mix,
        amp_mix = amp_mix,
        sink_mix = sink_mix,
        survived = survived,
    )

    return τ_grid, Cτ, stats
end


# --------------------------------------------
# Readability wrappers
# --------------------------------------------

simulate_coherence_correlator_from_samples_LHS_qs(args...; kwargs...) =
    simulate_coherence_correlator_from_samples_qs(args...; side = :LHS, kwargs...)

simulate_coherence_correlator_from_samples_RHS_qs(args...; kwargs...) =
    simulate_coherence_correlator_from_samples_qs(args...; side = :RHS, kwargs...)


# --------------------------------------------
# Commutator correlator from fixed stationary samples
# --------------------------------------------

function simulate_coherence_commutator_from_samples_qs(
    S_samples::AbstractVector{Int},
    M_samples::AbstractVector{Int},
    model::CoherenceModel;
    τ_max::Float64,
    t_mix::Float64,
    n_grid::Int = 1000,
    tail_mode::Symbol = :deterministic,
    rng::AbstractRNG = Random.default_rng(),
)
    τ_grid, C_LHS, stats_LHS = simulate_coherence_correlator_from_samples_qs(
        S_samples,
        M_samples,
        model;
        side = :LHS,
        τ_max = τ_max,
        t_mix = t_mix,
        n_grid = n_grid,
        tail_mode = tail_mode,
        rng = rng,
    )

    _, C_RHS, stats_RHS = simulate_coherence_correlator_from_samples_qs(
        S_samples,
        M_samples,
        model;
        side = :RHS,
        τ_max = τ_max,
        t_mix = t_mix,
        n_grid = n_grid,
        tail_mode = tail_mode,
        rng = rng,
    )

    C_comm = C_LHS .- C_RHS

    stats = (
        LHS = stats_LHS,
        RHS = stats_RHS,
    )

    return τ_grid, C_comm, C_LHS, C_RHS, stats
end


# --------------------------------------------
# Convenience wrappers that first sample the population steady state
# --------------------------------------------

function simulate_correlator_qs(
    pop_model::PopulationModel;
    side::Symbol,
    N_traj::Int,
    t_ss::Float64,
    τ_max::Float64,
    t_mix::Float64,
    n_grid::Int = 1000,
    S0::Int = pop_model.Jmax,
    M0::Int = pop_model.Jmax,
    tail_mode::Symbol = :deterministic,
    rng::AbstractRNG = Random.default_rng(),
)
    S_samples, M_samples = get_stationary_population_samples(
        pop_model,
        N_traj;
        t_ss = t_ss,
        S0 = S0,
        M0 = M0,
        rng = rng,
    )

    coh_model = CoherenceModel(pop_model)

    τ_grid, Cτ, stats = simulate_coherence_correlator_from_samples_qs(
        S_samples,
        M_samples,
        coh_model;
        side = side,
        τ_max = τ_max,
        t_mix = t_mix,
        n_grid = n_grid,
        tail_mode = tail_mode,
        rng = rng,
    )

    return τ_grid, Cτ, S_samples, M_samples, stats
end


function simulate_correlator_ergodic_qs(
    pop_model::PopulationModel;
    side::Symbol,
    N_traj::Int,
    t_ss::Float64,
    t_autocorrelation::Float64,
    τ_max::Float64,
    t_mix::Float64,
    n_grid::Int = 1000,
    S0::Int = pop_model.Jmax,
    M0::Int = pop_model.Jmax,
    tail_mode::Symbol = :deterministic,
    rng::AbstractRNG = Random.default_rng(),
)
    S_samples, M_samples = get_stationary_population_samples_ergodic(
        pop_model,
        N_traj;
        t_ss = t_ss,
        t_sample = t_autocorrelation,
        S0 = S0,
        M0 = M0,
        rng = rng,
    )

    coh_model = CoherenceModel(pop_model)

    τ_grid, Cτ, stats = simulate_coherence_correlator_from_samples_qs(
        S_samples,
        M_samples,
        coh_model;
        side = side,
        τ_max = τ_max,
        t_mix = t_mix,
        n_grid = n_grid,
        tail_mode = tail_mode,
        rng = rng,
    )

    return τ_grid, Cτ, S_samples, M_samples, stats
end


function simulate_commutator_qs(
    pop_model::PopulationModel;
    N_traj::Int,
    t_ss::Float64,
    τ_max::Float64,
    t_mix::Float64,
    n_grid::Int = 1000,
    S0::Int = pop_model.Jmax,
    M0::Int = pop_model.Jmax,
    tail_mode::Symbol = :deterministic,
    rng::AbstractRNG = Random.default_rng(),
)
    S_samples, M_samples = get_stationary_population_samples(
        pop_model,
        N_traj;
        t_ss = t_ss,
        S0 = S0,
        M0 = M0,
        rng = rng,
    )

    coh_model = CoherenceModel(pop_model)

    τ_grid, C_comm, C_LHS, C_RHS, stats = simulate_coherence_commutator_from_samples_qs(
        S_samples,
        M_samples,
        coh_model;
        τ_max = τ_max,
        t_mix = t_mix,
        n_grid = n_grid,
        tail_mode = tail_mode,
        rng = rng,
    )

    return τ_grid, C_comm, C_LHS, C_RHS, S_samples, M_samples, stats
end


function simulate_commutator_ergodic_qs(
    pop_model::PopulationModel;
    N_traj::Int,
    t_ss::Float64,
    t_autocorrelation::Float64,
    τ_max::Float64,
    t_mix::Float64,
    n_grid::Int = 1000,
    S0::Int = pop_model.Jmax,
    M0::Int = pop_model.Jmax,
    tail_mode::Symbol = :deterministic,
    rng::AbstractRNG = Random.default_rng(),
)
    S_samples, M_samples = get_stationary_population_samples_ergodic(
        pop_model,
        N_traj;
        t_ss = t_ss,
        t_sample = t_autocorrelation,
        S0 = S0,
        M0 = M0,
        rng = rng,
    )

    coh_model = CoherenceModel(pop_model)

    τ_grid, C_comm, C_LHS, C_RHS, stats = simulate_coherence_commutator_from_samples_qs(
        S_samples,
        M_samples,
        coh_model;
        τ_max = τ_max,
        t_mix = t_mix,
        n_grid = n_grid,
        tail_mode = tail_mode,
        rng = rng,
    )

    return τ_grid, C_comm, C_LHS, C_RHS, S_samples, M_samples, stats
end
