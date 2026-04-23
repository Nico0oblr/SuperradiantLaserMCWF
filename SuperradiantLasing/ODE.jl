"""
state is a tuple (S,M)
params is a dictionary (global_decay, global_pump, local_pump, local_decay, local_dephasing)
gives the outgoing rates and the states they point to
"""
function outgoing_for_state(state, params, Jmax)
    J, M = state
    out = DefaultDict{Tuple{Rational, Rational}, Float64}(0.0)
    # global_decay M -> M-1 and J->J
    if sane(J,M-1,Jmax) out[(J, M-1)] += A_JM_minus(J,M) * sqrt(params["global_decay"]) end
    # global_decay M -> M+1 and J->J
    if sane(J, M+1,Jmax) out[(J, M+1)] += A_JM_plus(J,M) * sqrt(params["global_pump"]) end
    # local_decay M -> M-1 and J->(J+1,J,J-1)
    if sane(J+1, M-1,Jmax) out[(J+1, M-1)] += P_JM_minus_plus(J,M,2*Jmax) * sqrt(params["local_decay"]) end
    if sane(J-1,M-1,Jmax) out[(J-1,M-1)] += P_JM_minus_minus(J,M,2*Jmax) * sqrt(params["local_decay"]) end
    if sane(J,M-1,Jmax) out[(J,M-1)] += P_JM_minus_0(J,M,2*Jmax) * sqrt(params["local_decay"]) end
    # local_pump M -> M+1 and J->(J+1,J,J-1)
    if sane(J+1,M+1,Jmax) out[(J+1,M+1)] += P_JM_plus_plus(J,M,2*Jmax) * sqrt(params["local_pump"]) end
    if sane(J-1,M+1,Jmax) out[(J-1,M+1)] += P_JM_plus_minus(J,M,2*Jmax) * sqrt(params["local_pump"]) end
    if sane(J,M+1,Jmax) out[(J,M+1)] += P_JM_plus_0(J,M,2*Jmax) * sqrt(params["local_pump"]) end 
    # local_dephasing M -> M and J->(J+1,J,J-1)
    if sane(J+1, M,Jmax) out[(J+1, M)] += P_JM_z_plus(J,M,2*Jmax) * sqrt(params["local_dephasing"]) end
    if sane(J-1, M,Jmax) out[(J-1, M)] += P_JM_z_minus(J,M,2*Jmax) * sqrt(params["local_dephasing"]) end 
    if sane(J, M,Jmax) out[(J, M)] += P_JM_z_0(J,M,2*Jmax) * sqrt(params["local_dephasing"]) end
    return out
end

function incoming_for_state(state, params, Jmax)
    J, M = state
    inp = DefaultDict{Tuple{Rational, Rational}, Float64}(0.0)

    if sane(J,M+1, Jmax); inp[(J,M+1)] += A_JM_minus(J,M+1)* sqrt(params["global_decay"]); end
    if sane(J,M-1, Jmax); inp[(J,M-1)] += A_JM_plus(J, M-1) * sqrt(params["global_pump"]); end

    if sane(J-1, M+1, Jmax); inp[(J-1,M+1)] += P_JM_minus_plus(J-1,M+1, 2*Jmax)  * sqrt(params["local_decay"]); end
    if sane(J+1, M+1, Jmax); inp[(J+1,M+1)] += P_JM_minus_minus(J+1,M+1, 2*Jmax)  * sqrt(params["local_decay"]); end
    if sane(J, M+1, Jmax); inp[(J, M+1)] += P_JM_minus_0(J,M+1, 2*Jmax)   * sqrt(params["local_decay"]); end

    if sane(J-1, M-1, Jmax); inp[(J-1, M-1)] += P_JM_plus_plus(J-1, M-1, 2*Jmax) * sqrt(params["local_pump"]); end
    if sane(J+1, M-1, Jmax); inp[(J+1, M-1)] += P_JM_plus_minus(J+1,M-1, 2*Jmax) * sqrt(params["local_pump"]); end
    if sane(J, M-1, Jmax); inp[(J, M-1)] += P_JM_plus_0(J,    M-1, 2*Jmax)* sqrt(params["local_pump"]); end

    if sane(J-1, M, Jmax); inp[(J-1, M)] += P_JM_z_plus(J-1, M, 2*Jmax) * sqrt(params["local_dephasing"]); end
    if sane(J+1, M, Jmax); inp[(J+1, M)] += P_JM_z_minus(J+1, M, 2*Jmax) * sqrt(params["local_dephasing"]); end
    if sane(J, M, Jmax); inp[(J, M)] += P_JM_z_0(J, M, 2*Jmax) * sqrt(params["local_dephasing"]); end

    return inp
end

function sum_incoming_for_state(state, params, Jmax)
    tmp = incoming_for_state(state, params, Jmax)
    return sum(x[2] ^ 2 for x in tmp)
end

function sum_outgoing_for_state(state, params, Jmax)
    tmp = outgoing_for_state(state, params, Jmax)
    return sum(x[2] ^ 2 for x in tmp)
end

function check_inversion(target, params, Jmax; atol=1e-10)
    incoming = incoming_for_state(target, params, Jmax)
    ok = true
    for (src, w_in) in incoming
        out = outgoing_for_state(src, params, Jmax)
        w_out = get(out, target, 0.0)
        ok &= isapprox(w_in, w_out; atol=atol)
        if !ok
            @warn "Mismatch" src=src target=target w_in=w_in w_out=w_out
        end
    end
    return ok
end

"""
For diagonal evolution, define the rate matrix
"""
function set_up_rate_matrix(Jmax, params)
    basis = build_basis(Jmax)
    index_mapping = Dict(basis .=> 1:length(basis))
    sparse_constructor = DefaultDict{Tuple{Int64, Int64}, Float64}(0.0)

    for state in basis
        out_rate = sum_outgoing_for_state(state, params, Jmax)
        in_rates = incoming_for_state(state, params, Jmax)
        is = index_mapping[state]
        sparse_constructor[(is, is)] -= out_rate
        for (src, rate) in in_rates
            sparse_constructor[(is, index_mapping[src])] += rate ^ 2
        end
    end


    I = Int[]; J = Int[]; V = Float64[]
    sizehint!(I, length(sparse_constructor)); sizehint!(J, length(sparse_constructor)); sizehint!(V, length(sparse_constructor))
    for ((i, j), v) in sparse_constructor
        push!(I, i); push!(J, j); push!(V, v)
    end
        
    N = length(basis)
    out = sparse(I,J,V,N,N)
    dropzeros!(out)
    return basis, index_mapping, out, N
end

function initial_state(index_mapping, state, nstates)
    u = zeros(Float64, nstates)
    u[index_mapping[state]] = 1.0
    return u
end

# ============================================
# Coherence-sector deterministic helpers
# Convention:
#   state = (S,M)  <->  |S,M-1><S,M|
# ============================================

function build_coherence_basis(Jmax)
    return [(S, M) for S in 0:Jmax for M in -S:S if valid_coherence_state(S, M, Jmax)]
end

# returns:
#   internal::DefaultDict{Tuple{Int,Int},Float64}
#   sink_rate::Float64
function outgoing_for_coherence_state(state, params, Jmax)
    S, M = state
    N = 2 * Jmax

    model = CoherenceModel(
        N,
        Jmax,
        params["global_decay"],
        params["global_pump"],
        params["local_decay"],
        params["local_pump"],
        params["local_dephasing"],
    )

    internal = DefaultDict{Tuple{Int,Int}, Float64}(0.0)
    sink_rate = 0.0

    # helper: add one physical branch using arithmetic/geometric mean rule
    function add_branch!(dS, dM, poprate_fun)
        rM   = poprate_fun(S, M, model)
        rMm1 = poprate_fun(S, M - 1, model)

        internal_rate, loss_rate = coherence_internal_and_loss(rM, rMm1)

        S2, M2 = S + dS, M + dM
        if valid_coherence_state(S2, M2, Jmax)
            internal[(S2, M2)] += internal_rate
            sink_rate += loss_rate - internal_rate
        else
            sink_rate += loss_rate
        end
        return nothing
    end

    # global decay
    if params["global_decay"] > 0
        add_branch!(0, -1, poprate_gdec)
    end

    # global pump
    if params["global_pump"] > 0
        add_branch!(0, +1, poprate_gpump)
    end

    # local decay
    if params["local_decay"] > 0
        add_branch!(+1, -1, poprate_ldec_plus)
        add_branch!( 0, -1, poprate_ldec_0)
        add_branch!(-1, -1, poprate_ldec_minus)
    end

    # local pump
    if params["local_pump"] > 0
        add_branch!(+1, +1, poprate_lpump_plus)
        add_branch!( 0, +1, poprate_lpump_0)
        add_branch!(-1, +1, poprate_lpump_minus)
    end

    # local dephasing
    if params["local_dephasing"] > 0
        add_branch!(+1, 0, poprate_deph_plus)
        add_branch!( 0, 0, poprate_deph_0)
        add_branch!(-1, 0, poprate_deph_minus)
    end

    if sink_rate < 0 && abs(sink_rate) < 1e-12
        sink_rate = 0.0
    end

    return internal, sink_rate
end

function incoming_for_coherence_state(state, params, Jmax)
    S, M = state
    incoming = DefaultDict{Tuple{Int,Int}, Float64}(0.0)

    # each possible source is obtained by reversing one allowed branch
    candidate_sources = (
        ((S,     M + 1),  0, -1, :gdec),
        ((S,     M - 1),  0, +1, :gpump),

        ((S - 1, M + 1), +1, -1, :ldec_plus),
        ((S,     M + 1),  0, -1, :ldec_0),
        ((S + 1, M + 1), -1, -1, :ldec_minus),

        ((S - 1, M - 1), +1, +1, :lpump_plus),
        ((S,     M - 1),  0, +1, :lpump_0),
        ((S + 1, M - 1), -1, +1, :lpump_minus),

        ((S - 1, M    ), +1,  0, :deph_plus),
        ((S,     M    ),  0,  0, :deph_0),
        ((S + 1, M    ), -1,  0, :deph_minus),
    )

    N = 2 * Jmax
    model = CoherenceModel(
        N,
        Jmax,
        params["global_decay"],
        params["global_pump"],
        params["local_decay"],
        params["local_pump"],
        params["local_dephasing"],
    )

    for ((Ss, Ms), dS, dM, tag) in candidate_sources
        if !valid_coherence_state(Ss, Ms, Jmax)
            continue
        end

        rate = 0.0

        if tag === :gdec && params["global_decay"] > 0
            rM   = poprate_gdec(Ss, Ms, model)
            rMm1 = poprate_gdec(Ss, Ms - 1, model)
            rate, _ = coherence_internal_and_loss(rM, rMm1)

        elseif tag === :gpump && params["global_pump"] > 0
            rM   = poprate_gpump(Ss, Ms, model)
            rMm1 = poprate_gpump(Ss, Ms - 1, model)
            rate, _ = coherence_internal_and_loss(rM, rMm1)

        elseif tag === :ldec_plus && params["local_decay"] > 0
            rM   = poprate_ldec_plus(Ss, Ms, model)
            rMm1 = poprate_ldec_plus(Ss, Ms - 1, model)
            rate, _ = coherence_internal_and_loss(rM, rMm1)

        elseif tag === :ldec_0 && params["local_decay"] > 0
            rM   = poprate_ldec_0(Ss, Ms, model)
            rMm1 = poprate_ldec_0(Ss, Ms - 1, model)
            rate, _ = coherence_internal_and_loss(rM, rMm1)

        elseif tag === :ldec_minus && params["local_decay"] > 0
            rM   = poprate_ldec_minus(Ss, Ms, model)
            rMm1 = poprate_ldec_minus(Ss, Ms - 1, model)
            rate, _ = coherence_internal_and_loss(rM, rMm1)

        elseif tag === :lpump_plus && params["local_pump"] > 0
            rM   = poprate_lpump_plus(Ss, Ms, model)
            rMm1 = poprate_lpump_plus(Ss, Ms - 1, model)
            rate, _ = coherence_internal_and_loss(rM, rMm1)

        elseif tag === :lpump_0 && params["local_pump"] > 0
            rM   = poprate_lpump_0(Ss, Ms, model)
            rMm1 = poprate_lpump_0(Ss, Ms - 1, model)
            rate, _ = coherence_internal_and_loss(rM, rMm1)

        elseif tag === :lpump_minus && params["local_pump"] > 0
            rM   = poprate_lpump_minus(Ss, Ms, model)
            rMm1 = poprate_lpump_minus(Ss, Ms - 1, model)
            rate, _ = coherence_internal_and_loss(rM, rMm1)

        elseif tag === :deph_plus && params["local_dephasing"] > 0
            rM   = poprate_deph_plus(Ss, Ms, model)
            rMm1 = poprate_deph_plus(Ss, Ms - 1, model)
            rate, _ = coherence_internal_and_loss(rM, rMm1)

        elseif tag === :deph_0 && params["local_dephasing"] > 0
            rM   = poprate_deph_0(Ss, Ms, model)
            rMm1 = poprate_deph_0(Ss, Ms - 1, model)
            rate, _ = coherence_internal_and_loss(rM, rMm1)

        elseif tag === :deph_minus && params["local_dephasing"] > 0
            rM   = poprate_deph_minus(Ss, Ms, model)
            rMm1 = poprate_deph_minus(Ss, Ms - 1, model)
            rate, _ = coherence_internal_and_loss(rM, rMm1)
        end

        # only keep sources that actually land on the target state
        if rate > 0
            incoming[(Ss, Ms)] += rate
        end
    end

    return incoming
end

function sum_outgoing_for_coherence_state(state, params, Jmax)
    internal, sink_rate = outgoing_for_coherence_state(state, params, Jmax)
    return sum(x[2] for x in internal) + sink_rate
end

function set_up_coherence_rate_matrix(Jmax, params)
    basis = build_coherence_basis(Jmax)
    index_mapping = Dict(basis .=> 1:length(basis))
    sparse_constructor = DefaultDict{Tuple{Int64, Int64}, Float64}(0.0)

    for state in basis
        out_rate = sum_outgoing_for_coherence_state(state, params, Jmax)
        in_rates = incoming_for_coherence_state(state, params, Jmax)

        is = index_mapping[state]
        sparse_constructor[(is, is)] -= out_rate

        for (src, rate) in in_rates
            sparse_constructor[(is, index_mapping[src])] += rate
        end
    end

    I = Int[]; J = Int[]; V = Float64[]
    sizehint!(I, length(sparse_constructor))
    sizehint!(J, length(sparse_constructor))
    sizehint!(V, length(sparse_constructor))

    for ((i, j), v) in sparse_constructor
        push!(I, i); push!(J, j); push!(V, v)
    end

    Nbasis = length(basis)
    out = sparse(I, J, V, Nbasis, Nbasis)
    dropzeros!(out)
    return basis, index_mapping, out, Nbasis
end


function steady_state_from_rate_matrix(R)
    A = Matrix(R)
    n = size(A, 1)

    # replace last row by normalization condition
    A[end, :] .= 1.0
    b = zeros(Float64, n)
    b[end] = 1.0

    pss = A \ b
    return pss
end

function coherence_initial_from_population_ss(pop_basis, pop_ss, coh_index_mapping)
    c0 = zeros(Float64, length(coh_index_mapping))
    for (i, (S, M)) in enumerate(pop_basis)
        if haskey(coh_index_mapping, (S, M))
            c0[coh_index_mapping[(S, M)]] = sqrt(A_JM_minus2(S, M)) * pop_ss[i]
        end
    end
    return c0
end

function coherence_readout_vector(basis)
    return [sqrt(A_JM_minus2(S, M)) for (S, M) in basis]
end
