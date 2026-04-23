# ============================================
# Coherence-sector numerics (no public interface yet)
# Convention:
#   coherence label (S,M) represents |S,M-1><S,M|
# ============================================

# --------------------------------------------
# State / buffer
# --------------------------------------------

mutable struct CoherenceState
    S::Int
    M::Int
    t::Float64
    absorbed::Bool
end

mutable struct CoherenceEventBuffer
    rates::Vector{Float64}
    dS::Vector{Int}
    dM::Vector{Int}
    n::Int
    sink_rate::Float64
end

function CoherenceEventBuffer(max_events::Int = 11)
    CoherenceEventBuffer(zeros(max_events), zeros(Int, max_events), zeros(Int, max_events), 0, 0.0)
end

@inline function valid_coherence_state(S::Int, M::Int, Jmax::Int)
    # (S,M) labels |S,M-1><S,M|, so both M and M-1 must lie in [-S,S]
    return 0 <= S <= Jmax && abs(M) <= S && abs(M - 1) <= S
end

@inline function push_coherence_event!(buf::CoherenceEventBuffer, rate::Float64, dS::Int, dM::Int)
    if rate <= 0.0
        return
    end
    n = buf.n + 1
    buf.rates[n] = rate
    buf.dS[n] = dS
    buf.dM[n] = dM
    buf.n = n
end

# --------------------------------------------
# Model
# --------------------------------------------

struct CoherenceModel{T} <: AbstractSpinModel
    N::Int
    Jmax::Int
    global_decay::T
    global_pump::T
    local_decay::T
    local_pump::T
    local_dephasing::T
end

CoherenceModel(model::PopulationModel) = CoherenceModel(
    model.N,
    model.Jmax,
    model.global_decay,
    model.global_pump,
    model.local_decay,
    model.local_pump,
    model.local_dephasing,
)

# --------------------------------------------
# Bookkeeping factors for S rho_ss and Tr(S† rho')
# with convention (S,M) <-> |S,M-1><S,M|
# --------------------------------------------

@inline initial_coherence_weight(S::Int, M::Int) = sqrt(A_JM_minus2(S, M))
@inline coherence_readout(S::Int, M::Int) = sqrt(A_JM_minus2(S, M))

# --------------------------------------------
# Channel-by-channel helpers
#
# For a given population-sector branch rate r(S,M),
# coherence-sector:
#   internal rate  = sqrt(r(S,M) * r(S,M-1))
#   diagonal loss  = (r(S,M) + r(S,M-1))/2
#   sink contribution = diagonal - internal
# --------------------------------------------

@inline function coherence_internal_and_loss(rM::Float64, rMm1::Float64)
    internal = sqrt(rM * rMm1)
    loss = 0.5 * (rM + rMm1)
    return internal, loss
end

# --------------------------------------------
# Population-branch rates evaluated on a single side
# These are the r(S,M) objects entering the rule above
# --------------------------------------------



# --------------------------------------------
# Add one coherence-sector branch from one physical process
#
# If population branch maps
#   |S,M>   -> |S+dS, M+dM>
# then coherence label (S,M) = |S,M-1><S,M|
# maps internally to
#   (S+dS, M+dM)
# because:
#   |S,M-1><S,M| -> |S+dS,M-1+dM><S+dS,M+dM|
# which is exactly the same offset convention.
# --------------------------------------------

function add_coherence_branch!(
    buf::CoherenceEventBuffer,
    model::CoherenceModel,
    S::Int,
    M::Int,
    dS::Int,
    dM::Int,
    poprate_fun,
)
    rM = poprate_fun(S, M, model)
    rMm1 = poprate_fun(S, M - 1, model)

    internal, loss = coherence_internal_and_loss(rM, rMm1)

    # internal target in the same coherence block
    S2 = S + dS
    M2 = M + dM
    if valid_coherence_state(S2, M2, model.Jmax)
        push_coherence_event!(buf, internal, dS, dM)
        buf.sink_rate += loss - internal
    else
        # if the target coherence state is outside the tracked block,
        # the whole diagonal loss contributes to sink
        buf.sink_rate += loss
    end

    return nothing
end

# --------------------------------------------
# Build all coherence-sector outgoing rates
# --------------------------------------------

function build_coherence_events!(buf::CoherenceEventBuffer, model::CoherenceModel, st::CoherenceState)
    buf.n = 0
    buf.sink_rate = 0.0

    S = st.S
    M = st.M

    # global decay
    if model.global_decay > 0
        add_coherence_branch!(buf, model, S, M, 0, -1, poprate_gdec)
    end

    # global pump
    if model.global_pump > 0
        add_coherence_branch!(buf, model, S, M, 0, +1, poprate_gpump)
    end

    # local decay
    if model.local_decay > 0
        add_coherence_branch!(buf, model, S, M, +1, -1, poprate_ldec_plus)
        add_coherence_branch!(buf, model, S, M,  0, -1, poprate_ldec_0)
        add_coherence_branch!(buf, model, S, M, -1, -1, poprate_ldec_minus)
    end

    # local pump
    if model.local_pump > 0
        add_coherence_branch!(buf, model, S, M, +1, +1, poprate_lpump_plus)
        add_coherence_branch!(buf, model, S, M,  0, +1, poprate_lpump_0)
        add_coherence_branch!(buf, model, S, M, -1, +1, poprate_lpump_minus)
    end

    # local dephasing
    if model.local_dephasing > 0
        add_coherence_branch!(buf, model, S, M, +1, 0, poprate_deph_plus)
        add_coherence_branch!(buf, model, S, M, 0, 0, poprate_deph_0)        
        add_coherence_branch!(buf, model, S, M, -1, 0, poprate_deph_minus)
    end

    # numerical safety
    if buf.sink_rate < 0 && abs(buf.sink_rate) < 1e-12
        buf.sink_rate = 0.0
    end

    return nothing
end

# --------------------------------------------
# Diagnostics
# --------------------------------------------

function sum_coherence_internal_rates(buf::CoherenceEventBuffer)
    s = 0.0
    @inbounds for k in 1:buf.n
        s += buf.rates[k]
    end
    return s
end

function inspect_coherence_rates(model::CoherenceModel, S::Int, M::Int)
    st = CoherenceState(S, M, 0.0, false)
    buf = CoherenceEventBuffer()
    build_coherence_events!(buf, model, st)

    println("state = ", (S, M), "  meaning |$S,$(M-1)><$S,$M|")
    println("n internal = ", buf.n)
    for k in 1:buf.n
        println("  -> ", (S + buf.dS[k], M + buf.dM[k]), "   rate = ", buf.rates[k])
    end
    println("sink_rate = ", buf.sink_rate)
    println("sum_internal = ", sum_coherence_internal_rates(buf))
    println("total_out = ", sum_coherence_internal_rates(buf) + buf.sink_rate)
    return nothing
end

# ============================================
# One coherence-sector Gillespie step
# ============================================

function coherence_step!(
    rng::AbstractRNG,
    st::CoherenceState,
    buf::CoherenceEventBuffer,
    model::CoherenceModel,
)
    if st.absorbed
        return false
    end

    build_coherence_events!(buf, model, st)

    r_internal = sum_coherence_internal_rates(buf)
    rtot = r_internal + buf.sink_rate

    if rtot <= 0.0
        return false
    end

    st.t += -log(rand(rng)) / rtot

    x = rand(rng) * rtot
    acc = 0.0

    @inbounds for k in 1:buf.n
        acc += buf.rates[k]
        if x <= acc
            st.S += buf.dS[k]
            st.M += buf.dM[k]
            return true
        end
    end

    st.absorbed = true
    return true
end


# ============================================
# Simulate one coherence trajectory on a grid
# Returns the contribution to C(τ) from one sample
# ============================================

function simulate_coherence_sector(
    model::CoherenceModel,
    S0::Int,
    M0::Int,
    τ_grid::AbstractVector;
    rng::AbstractRNG = Random.default_rng(),
)
    st = CoherenceState(S0, M0, 0.0, false)
    buf = CoherenceEventBuffer()

    Tacc = Float64
    contrib = zeros(Tacc, length(τ_grid))

    w_in = initial_coherence_weight(S0, M0)

    if !valid_coherence_state(S0, M0, model.Jmax)
        return contrib
    end

    idx = 1
    while idx <= length(τ_grid)
        current_val = st.absorbed ? 0.0 : w_in * coherence_readout(st.S, st.M)

        ok = coherence_step!(rng, st, buf, model)

        if !ok
            @inbounds while idx <= length(τ_grid)
                contrib[idx] = current_val
                idx += 1
            end
            break
        end

        @inbounds while idx <= length(τ_grid) && τ_grid[idx] <= st.t
            contrib[idx] = current_val
            idx += 1
        end

        if st.absorbed
            @inbounds while idx <= length(τ_grid)
                contrib[idx] = 0.0
                idx += 1
            end
            break
        end
    end

    return contrib
end

function simulate_coherence_correlator_from_samples(
    S_samples::AbstractVector{Int},
    M_samples::AbstractVector{Int},
    model::CoherenceModel;
    τ_max::Float64,
    n_grid::Int = 1000,
    rng::AbstractRNG = Random.default_rng(),
)
    @assert length(S_samples) == length(M_samples)

    τ_grid = collect(range(0.0, τ_max, length = n_grid))
    Cτ = zeros(Float64, n_grid)

    for n in eachindex(S_samples)
        Cτ .+= simulate_coherence_sector(
            model,
            S_samples[n],
            M_samples[n],
            τ_grid;
            rng = rng,
        )
    end

    Cτ ./= length(S_samples)
    return τ_grid, Cτ
end

function simulate_correlator(
    pop_model::PopulationModel;
    N_traj::Int,
    t_ss::Float64,
    τ_max::Float64,
    n_grid::Int = 1000,
    S0::Int = pop_model.Jmax,
    M0::Int = pop_model.Jmax,
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

    τ_grid, Cτ = simulate_coherence_correlator_from_samples(
        S_samples,
        M_samples,
        coh_model;
        τ_max = τ_max,
        n_grid = n_grid,
        rng = rng,
    )

    return τ_grid, Cτ, S_samples, M_samples
end


function simulate_correlator_ergodic(
    pop_model::PopulationModel;
    N_traj::Int,
    t_ss::Float64,
    τ_max::Float64,
    n_grid::Int = 1000,
    S0::Int = pop_model.Jmax,
    M0::Int = pop_model.Jmax,
    rng::AbstractRNG = Random.default_rng(),
)
    S_samples, M_samples = get_stationary_population_samples_ergodic(
        pop_model,
        N_traj;
        t_ss = t_ss,
        t_sample = t_ss * 3,
        S0 = S0,
        M0 = M0,
        rng = rng,
    )

    coh_model = CoherenceModel(pop_model)

    τ_grid, Cτ = simulate_coherence_correlator_from_samples(
        S_samples,
        M_samples,
        coh_model;
        τ_max = τ_max,
        n_grid = n_grid,
        rng = rng,
    )

    return τ_grid, Cτ, S_samples, M_samples
end