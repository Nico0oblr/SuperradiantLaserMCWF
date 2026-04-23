struct PopulationModel{T}
    N::Int
    Jmax::Int
    global_decay::T
    global_pump::T
    local_decay::T
    local_pump::T
    local_dephasing::T
end

mutable struct PopState
    S::Int
    M::Int
    t::Float64
end

mutable struct EventBuffer
    rates::Vector{Float64}
    dS::Vector{Int}
    dM::Vector{Int}
    n::Int
end

function EventBuffer(max_events::Int = 11)
    return EventBuffer(zeros(max_events), zeros(Int, max_events), zeros(Int, max_events), 0)
end

@inline function valid_state(S::Int, M::Int, Jmax::Int)
    return 0 <= S <= Jmax && abs(M) <= S
end

@inline function push_event!(buf::EventBuffer, rate::Float64, dS::Int, dM::Int)
    if rate <= 0.0
        return
    end
    n = buf.n + 1
    buf.rates[n] = rate
    buf.dS[n] = dS
    buf.dM[n] = dM
    buf.n = n
end


# =============================================
# Build active event list for current population
# =============================================

function build_events!(buf::EventBuffer, model::PopulationModel, st::PopState)
    buf.n = 0

    S = st.S
    M = st.M
    N = model.N
    Jmax = model.Jmax

    # global decay: (S,M) -> (S, M-1)
    if model.global_decay > 0
        if valid_state(S, M - 1, Jmax)
            push_event!(buf, model.global_decay * A_JM_minus2(S, M), 0, -1)
        end
    end

    # global pump: (S,M) -> (S, M+1)
    if model.global_pump > 0
        if valid_state(S, M + 1, Jmax)
            push_event!(buf, model.global_pump * A_JM_plus2(S, M), 0, +1)
        end
    end

    # local decay: (S,M) -> (S+1,M-1), (S,M-1), (S-1,M-1)
    if model.local_decay > 0
        if valid_state(S + 1, M - 1, Jmax)
            push_event!(buf, model.local_decay * P_JM_minus_plus2(S, M, N), +1, -1)
        end
        if valid_state(S, M - 1, Jmax)
            push_event!(buf, model.local_decay * P_JM_minus_02(S, M, N), 0, -1)
        end
        if valid_state(S - 1, M - 1, Jmax)
            push_event!(buf, model.local_decay * P_JM_minus_minus2(S, M, N), -1, -1)
        end
    end

    # local pump: (S,M) -> (S+1,M+1), (S,M+1), (S-1,M+1)
    if model.local_pump > 0
        if valid_state(S + 1, M + 1, Jmax)
            push_event!(buf, model.local_pump * P_JM_plus_plus2(S, M, N), +1, +1)
        end
        if valid_state(S, M + 1, Jmax)
            push_event!(buf, model.local_pump * P_JM_plus_02(S, M, N), 0, +1)
        end
        if valid_state(S - 1, M + 1, Jmax)
            push_event!(buf, model.local_pump * P_JM_plus_minus2(S, M, N), -1, +1)
        end
    end

    # local dephasing: (S,M) -> (S+1,M), (S,M), (S-1,M)
    if model.local_dephasing > 0
        if valid_state(S + 1, M, Jmax)
            push_event!(buf, model.local_dephasing * P_JM_z_plus2(S, M, N), +1, 0)
        end
        if valid_state(S - 1, M, Jmax)
            push_event!(buf, model.local_dephasing * P_JM_z_minus2(S, M, N), -1, 0)
        end
    end

    return nothing
end


# =====================
# One Gillespie step
# =====================

function gillespie_step!(rng::AbstractRNG, st::PopState, buf::EventBuffer, model::PopulationModel)
    build_events!(buf, model, st)

    n = buf.n
    if n == 0
        return false
    end

    rtot = 0.0
    @inbounds for k in 1:n
        rtot += buf.rates[k]
    end

    if rtot <= 0.0
        return false
    end

    st.t += -log(rand(rng)) / rtot

    x = rand(rng) * rtot
    acc = 0.0
    chosen = n
    @inbounds for k in 1:n
        acc += buf.rates[k]
        if x <= acc
            chosen = k
            break
        end
    end

    @inbounds begin
        st.S += buf.dS[chosen]
        st.M += buf.dM[chosen]
    end

    return true
end


# =====================
# Trajectory driver
# =====================

function run_population_trajectory(
    model::PopulationModel,
    observable;
    S0::Int = model.Jmax,
    M0::Int = model.Jmax,
    tmax::Float64 = Inf,
    max_steps::Int = typemax(Int),
    stop_condition = nothing,
    rng::AbstractRNG = Random.default_rng(),
    sizehint_steps::Int = 1024,
)
    st = PopState(S0, M0, 0.0)
    buf = EventBuffer()

    times = Float64[]
    vals = Vector{typeof(observable(st.S, st.M, model))}()

    sizehint!(times, sizehint_steps)
    sizehint!(vals, sizehint_steps)

    push!(times, st.t)
    push!(vals, observable(st.S, st.M, model))

    steps = 0
    while st.t < tmax && steps < max_steps
        if stop_condition !== nothing && stop_condition(st.S, st.M, st.t, model)
            break
        end

        ok = gillespie_step!(rng, st, buf, model)
        if !ok
            break
        end

        push!(times, st.t)
        push!(vals, observable(st.S, st.M, model))
        steps += 1
    end

    return times, vals
end

function simulate_population_ensemble(
    model::PopulationModel,
    observable,
    N_traj::Int;
    S0::Int = model.Jmax,
    M0::Int = model.Jmax,
    t_max::Float64,
    n_grid::Int = 1000,
    rng::AbstractRNG = Random.default_rng(),
)
    time_grid = collect(range(0.0, t_max, length = n_grid))

    # Infer observable type from initial state
    #Tobs = typeof(observable(S0, M0, model))
    avg_obs = zeros(Float64, n_grid)

    for n in 1:N_traj
        ts, vals = run_population_trajectory(
            model,
            observable;
            S0 = S0,
            M0 = M0,
            tmax = t_max,
            rng = rng,
        )

        idx = 1
        @inbounds for k in eachindex(time_grid)
            t_point = time_grid[k]
            while idx < length(ts) && ts[idx + 1] <= t_point
                idx += 1
            end
            avg_obs[k] += vals[idx]
        end
    end

    avg_obs ./= N_traj
    return time_grid, avg_obs
end

function get_stationary_population_samples(
    model::PopulationModel,
    N_traj::Int;
    t_ss::Float64,
    S0::Int = model.Jmax,
    M0::Int = model.Jmax,
    rng::AbstractRNG = Random.default_rng(),
)
    S_samples = Vector{Int}(undef, N_traj)
    M_samples = Vector{Int}(undef, N_traj)

    for n in 1:N_traj
        st = PopState(S0, M0, 0.0)
        buf = EventBuffer()

        while st.t < t_ss
            ok = gillespie_step!(rng, st, buf, model)
            if !ok
                break
            end
        end

        S_samples[n] = st.S
        M_samples[n] = st.M
    end

    return S_samples, M_samples
end

function get_stationary_population_samples_ergodic(
    model::PopulationModel,
    N_traj::Int;
    t_ss::Float64,
    t_sample::Float64,
    S0::Int = model.Jmax,
    M0::Int = model.Jmax,
    rng::AbstractRNG = Random.default_rng(),
)
    S_samples = Vector{Int}(undef, N_traj)
    M_samples = Vector{Int}(undef, N_traj)

    # Random sample times in the post-burn-in window, sorted ascending
    sample_times = t_ss .+ rand(rng, N_traj) .* t_sample
    p = sortperm(sample_times)
    sample_times = sample_times[p]

    st = PopState(S0, M0, 0.0)
    buf = EventBuffer()

    i = 1
    while i <= N_traj
        build_events!(buf, model, st)
        n = buf.n

        # If no more jumps are possible, the state stays frozen forever
        if n == 0
            while i <= N_traj
                S_samples[i] = st.S
                M_samples[i] = st.M
                i += 1
            end
            break
        end

        rtot = 0.0
        @inbounds for k in 1:n
            rtot += buf.rates[k]
        end

        if rtot <= 0.0
            while i <= N_traj
                S_samples[i] = st.S
                M_samples[i] = st.M
                i += 1
            end
            break
        end

        dt = -log(rand(rng)) / rtot
        t_next = st.t + dt

        # Current state is occupied on [st.t, t_next)
        while i <= N_traj && sample_times[i] < t_next
            if sample_times[i] >= t_ss
                S_samples[i] = st.S
                M_samples[i] = st.M
            end
            i += 1
        end

        # Perform the jump
        x = rand(rng) * rtot
        acc = 0.0
        chosen = n
        @inbounds for k in 1:n
            acc += buf.rates[k]
            if x <= acc
                chosen = k
                break
            end
        end

        @inbounds begin
            st.S += buf.dS[chosen]
            st.M += buf.dM[chosen]
            st.t = t_next
        end
    end

    # Undo sorting so output order matches original random draws
    invp = invperm(p)
    return S_samples[invp], M_samples[invp]
end

obs_exc(S, M, model) = M + model.N / 2
obs_intensity(S, M, model) = (S + M) * (S - M + 1)
obs_spin(S, M, model) = S