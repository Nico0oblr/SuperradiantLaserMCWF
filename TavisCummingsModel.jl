using SparseArrays
using QuantumOptics
using QuantumCumulants
using SpecialFunctions

"""
    block_ranges(space::SumBasis)

Given a direct sum basis `space`, returns a vector of UnitRange{Int},
where each element gives the global index range corresponding to each subspace.
"""
function block_ranges(space::SumBasis)
    ranges = Vector{UnitRange{Int}}(undef, length(space.shape))
    start = 1
    for (i, len) in enumerate(space.shape)
        ranges[i] = start:(start+len-1)
        start += len
    end
    return ranges
end


"""
    build_operator(N, matrix_element)

Constructs a sparse operator over the Dicke Hilbert space (direct sum of SpinBasis blocks)
for a system of N atoms. The input `matrix_element` is a function with signature

    matrix_element(N, J, M, Jp, Mp) -> ComplexF64

which returns the matrix element connecting state |J, M⟩ in one block to |Jp, Mp⟩ in another.
Only nonzero values are stored.

Returns a tuple (op, space) where `op` is the operator (sparse matrix) and `space` is
the constructed SumBasis.
"""
function build_operator(N::Int, matrix_element::Function)
    # Construct the Dicke Hilbert space as a direct sum of SpinBasis blocks.
    subs = dicke_subspaces(N)         # vector of SpinBasis for each allowed J
    space = directsum(subs...)          # the full SumBasis
    total_dim = length(space)           # total dimension of the Hilbert space
    op = spzeros(ComplexF64, total_dim, total_dim)
    
    # Obtain global index ranges for each block in the direct sum.
    br = block_ranges(space)
    
    # Loop over all pairs of blocks: source block (with spin J) and target block (with spin Jp)
    for (i, b_source) in enumerate(subs)
        J = b_source.spinnumber
        dim_source = length(b_source)   # dimension = 2J + 1
        for (j, b_target) in enumerate(subs)
            Jp = b_target.spinnumber
            dim_target = length(b_target)  # dimension = 2Jp + 1
            
            # Get the global index ranges for these blocks.
            range_source = br[i]
            range_target = br[j]
            
            # Loop over basis states in the source block.
            # Convention: basis index k corresponds to M = J - (k - 1)
            for k_source in 1:dim_source
                M = J - (k_source - 1)
                # And similarly for the target block:
                for k_target in 1:dim_target
                    Mp = Jp - (k_target - 1)
                    val = matrix_element(N, J, M, Jp, Mp)
                    if val != 0
                        global_source = range_source[k_source]
                        global_target = range_target[k_target]
                        op[global_target, global_source] = val
                    end
                end
            end
        end
    end
    return SparseOperator(space, op)
end

# Collective factors (used also in the individual channels):
"""
    A_JM_minus(J, M)

Collective lowering amplitude:
  A_{JM}⁻ = √[(J + M)(J – M + 1)]
"""
A_JM_minus(J, M) = sqrt((J + M) * (J - M + 1))

"""
    A_JM_plus(J, M)

Collective raising amplitude:
  A_{JM}⁺ = √[(J – M)(J + M + 1)]
"""
A_JM_plus(J, M) = sqrt((J - M) * (J + M + 1))

# ---------------------------
# Individual Decay (lowering M)
# ---------------------------

"""
    P_JM_minus_0(J, M, N)

Individual decay channel with no change in J:
  P_{JM}^{-,0} = √[(2+N)/(4J(J+1))] A_{JM}⁻.
"""
P_JM_minus_0(J, M, N) = sqrt((2 + N) / (4 * J * (J + 1))) * A_JM_minus(J, M)

"""
    P_JM_minus_minus(J, M, N)

Individual decay channel lowering J by 1 (“s = -1”):
  P_{JM}^{-,-} = -√[(N+2J+2)(J+M)(J+M-1)/(4J(2J+1))].
"""
P_JM_minus_minus(J, M, N) = -sqrt((N + 2 * J + 2) * (J + M) * (J + M - 1) / (4 * J * (2 * J + 1)))

"""
    P_JM_minus_plus(J, M, N)

Individual decay channel increasing J by 1 (“s = +1”):
  P_{JM}^{-,+} = √[(N-2J)(J-M+1)(J-M+2)/(4(J+1)(2J+1))].
"""
P_JM_minus_plus(J, M, N) = sqrt((N - 2 * J) * (J - M + 1) * (J - M + 2) / (4 * (J + 1) * (2 * J + 1)))

# ---------------------------
# Individual Pumping (raising M)
# ---------------------------

"""
    P_JM_plus_0(J, M, N)

Individual pumping channel with no change in J:
  P_{JM}^{+,0} = √[(2+N)/(4J(J+1))] A_{JM}⁺.
"""
P_JM_plus_0(J, M, N) = sqrt((2 + N) / (4 * J * (J + 1))) * A_JM_plus(J, M)

"""
    P_JM_plus_minus(J, M, N)

Individual pumping channel lowering J by 1 (“s = -1”):
  P_{JM}^{+,-} = √[(N+2J+2)(J-M)(J-M-1)/(4J(2J+1))].
"""
P_JM_plus_minus(J, M, N) = sqrt((N + 2 * J + 2) * (J - M) * (J - M - 1) / (4 * J * (2 * J + 1)))

"""
    P_JM_plus_plus(J, M, N)

Individual pumping channel increasing J by 1 (“s = +1”):
  P_{JM}^{+,+} = -√[(N-2J)(J+M+1)(J+M+2)/(4(J+1)(2J+1))].
"""
P_JM_plus_plus(J, M, N) = -sqrt((N - 2 * J) * (J + M + 1) * (J + M + 2) / (4 * (J + 1) * (2 * J + 1)))

# ---------------------------
# Individual Dephasing
# ---------------------------

"""
    P_JM_z_0(J, M, N)

Individual dephasing channel that leaves M unchanged:
  P_{JM}^{z,0} = √[(2+N)/(4J(J+1))] M.
"""
function P_JM_z_0(J, M, N)
    if J == 0
        return 0.0
    else
        return sqrt((2 + N) / (4 * J * (J + 1))) * M
    end
end
#P_JM_z_0(J, M, N) = sqrt((2 + N) / (4 * J * (J + 1))) * M

"""
    P_JM_z_minus(J, M, N)

Individual dephasing channel lowering J by 1 (“s = -”):
  P_{JM}^{z,-} = √[(N+2J+2)(J-M)(J+M)/(4J(2J+1))].
"""
P_JM_z_minus(J, M, N) = sqrt((N + 2 * J + 2) * (J - M) * (J + M) / (4 * J * (2 * J + 1)))

"""
    P_JM_z_plus(J, M, N)

Individual dephasing channel increasing J by 1 (“s = +”):
  P_{JM}^{z,+} = √[(N-2J)(J+1-M)(J+1+M)/(4(J+1)(2J+1))].
"""
P_JM_z_plus(J, M, N) = sqrt((N - 2 * J) * (J + 1 - M) * (J + 1 + M) / (4 * (J + 1) * (2 * J + 1)))


"""
d_N_J(N, J)

Returns the degeneracy factor for a Dicke state with N atoms and total spin J,
defined by

  d_N^J = (N! (2J+1)) / ((N/2 - J)! (N/2 + J + 1)!)
      = Γ(N+1) (2J+1) / (Γ(N/2 - J + 1) Γ(N/2 + J + 2))
"""
function d_N_J(N::Int, J::Real)
    return 1.0 #gamma(N + 1) * (2 * J + 1) / (gamma(N/2 - J + 1) * gamma(N/2 + J + 2))
end

allowed_J(N) = (N % 2 == 0) ? collect(0:1:N÷2) : collect(1//2:1:N÷2)


# Create a list of SpinBasis subspaces for the allowed total spins
function dicke_subspaces(N::Int)
    js = (N % 2 == 0) ? [j for j in 0:N÷2] : [j//2 for j in 1:2:N]
    return reverse([SpinBasis(j) for j in js])
end

# Full Dicke Hilbert space as a direct sum
function dicke_space(N::Int)
    return directsum(dicke_subspaces(N)...)
end

# -- Matrix Element Functions for Jump Operators --

# 1. Collective Decay (acts within the same J, lowering M by 1)
function matrix_element_collective_decay(N, J, M, Jp, Mp)
    if J == Jp && Mp == M - 1
        # Multiply by sqrt(degeneracy) for the source block.
        return sqrt(d_N_J(N, J)) * A_JM_minus(J, M)
    else
        return 0.0
    end
end

# 2. Individual Decay, s = 0 channel (no change in J, lowering M by 1)
function matrix_element_individual_decay_s0(N, J, M, Jp, Mp)
    if J == Jp && Mp == M - 1
        return sqrt(d_N_J(N, J)) * P_JM_minus_0(J, M, N)
    else
        return 0.0
    end
end

# 3. Individual Decay, s = -1 channel (J -> J-1, lowering M by 1)
function matrix_element_individual_decay_sminus(N, J, M, Jp, Mp)
    if Jp == J - 1 && Mp == M - 1
        return sqrt(d_N_J(N, J)) * P_JM_minus_minus(J, M, N)
    else
        return 0.0
    end
end

# 4. Individual Decay, s = +1 channel (J -> J+1, lowering M by 1)
function matrix_element_individual_decay_splus(N, J, M, Jp, Mp)
    if Jp == J + 1 && Mp == M - 1
        return sqrt(d_N_J(N, J)) * P_JM_minus_plus(J, M, N)
    else
        return 0.0
    end
end


# 1. Individual Pumping Operators:

# s = 0: No change in J, M → M+1.
function matrix_element_individual_pump_s0(N, J, M, Jp, Mp)
    if J == Jp && Mp == M + 1
        return sqrt(d_N_J(N, J)) * P_JM_plus_0(J, M, N)
    else
        return 0.0
    end
end

# s = -1: J → J - 1, M → M+1.
function matrix_element_individual_pump_sminus(N, J, M, Jp, Mp)
    if Jp == J - 1 && Mp == M + 1
        return sqrt(d_N_J(N, J)) * P_JM_plus_minus(J, M, N)
    else
        return 0.0
    end
end

# s = +1: J → J + 1, M → M+1.
function matrix_element_individual_pump_splus(N, J, M, Jp, Mp)
    if Jp == J + 1 && Mp == M + 1
        return sqrt(d_N_J(N, J)) * P_JM_plus_plus(J, M, N)
    else
        return 0.0
    end
end

# 2. Individual Dephasing Operators:

# s = 0: No change in J, M unchanged.
function matrix_element_individual_dephase_s0(N, J, M, Jp, Mp)
    if J == Jp && Mp == M
        return sqrt(d_N_J(N, J)) * P_JM_z_0(J, M, N)
    else
        return 0.0
    end
end

# s = -: J → J - 1, M unchanged.
function matrix_element_individual_dephase_sminus(N, J, M, Jp, Mp)
    if Jp == J - 1 && Mp == M
        return sqrt(d_N_J(N, J)) * P_JM_z_minus(J, M, N)
    else
        return 0.0
    end
end

# s = +: J → J + 1, M unchanged.
function matrix_element_individual_dephase_splus(N, J, M, Jp, Mp)
    if Jp == J + 1 && Mp == M
        return sqrt(d_N_J(N, J)) * P_JM_z_plus(J, M, N)
    else
        return 0.0
    end
end

QuantumOptics.sigmam(hilbert_space::SumBasis) = directsum(sigmam.(hilbert_space.bases)...)
QuantumOptics.sigmap(hilbert_space::SumBasis) = directsum(sigmap.(hilbert_space.bases)...)
QuantumOptics.sigmax(hilbert_space::SumBasis) = directsum(sigmax.(hilbert_space.bases)...)
QuantumOptics.sigmay(hilbert_space::SumBasis) = directsum(sigmay.(hilbert_space.bases)...)
QuantumOptics.sigmaz(hilbert_space::SumBasis) = directsum(sigmaz.(hilbert_space.bases)...)
QuantumOptics.identityoperator(hilbert_space::SumBasis) = directsum(identityoperator.(hilbert_space.bases)...)

function build_collapse_operators(N::Int, Γ_c, γ_l, κ_l, d_l)
    # Get the common Hilbert space for N atoms.
    common_space = dicke_space(N)
    
    # Initialize empty lists for operators and rates.
    op_list::Vector = []
    rate_list::Vector = []
    
    # 1. Collective decay operator (rate: Γ_c)
    push!(op_list, build_operator(N, matrix_element_collective_decay))
    push!(rate_list, Γ_c)
    
    # 2. Individual decay operators (rate: γ_l)
    push!(op_list, build_operator(N, matrix_element_individual_decay_s0))
    push!(rate_list, γ_l)
    
    push!(op_list, build_operator(N, matrix_element_individual_decay_sminus))
    push!(rate_list, γ_l)
    
    push!(op_list, build_operator(N, matrix_element_individual_decay_splus))
    push!(rate_list, γ_l)
    
    # 3. Individual pumping operators (rate: κ_l)
    push!(op_list, build_operator(N, matrix_element_individual_pump_s0))
    push!(rate_list, κ_l)
    
    push!(op_list, build_operator(N, matrix_element_individual_pump_sminus))
    push!(rate_list, κ_l)
    
    push!(op_list, build_operator(N, matrix_element_individual_pump_splus))
    push!(rate_list, κ_l)
    
    # 4. Individual dephasing operators (rate: d_l)
    push!(op_list, build_operator(N, matrix_element_individual_dephase_s0))
    push!(rate_list, d_l)
    
    push!(op_list, build_operator(N, matrix_element_individual_dephase_sminus))
    push!(rate_list, d_l)
    
    push!(op_list, build_operator(N, matrix_element_individual_dephase_splus))
    push!(rate_list, d_l)
    
    return common_space, op_list, rate_list
end

"""
    build_tensor_collapse_operators(N, Γ_c, γ_l, κ_l, d_l)

Constructs the full tensor‑product Hilbert space for N two‑level atoms and returns:
  - the tensor Hilbert space,
  - a list of collapse operators,
  - a list of their corresponding rates.
  
These operators can then be used to simulate the master equation in the full tensor space.
"""
function build_tensor_collapse_operators(N::Int, Γ_c, γ_l, κ_l, d_l)
    # Single-atom basis for a two-level system.
    b_atom = SpinBasis(1//2)
    # Full tensor-product Hilbert space:
    b_tensor = tensor([b_atom for i=1:N]...)
    
    # Define helper functions to embed single-atom operators into the full space:
    sm(i) = embed(b_tensor, i, sigmam(b_atom))
    sp(i) = embed(b_tensor, i, sigmap(b_atom))
    sz(i) = embed(b_tensor, i, sigmaz(b_atom))
    
    # 1. Collective decay operator: S_- = sum_i σ_-^(i)
    S_minus = sum(sm.(1:N))
    
    # 2. Individual decay operators: one per atom.
    individual_decay_ops = [ sm(i) for i in 1:N ]
    
    # 3. Individual pumping operators: one per atom.
    individual_pump_ops = [ sp(i) for i in 1:N ]
    
    # 4. Individual dephasing operators: one per atom.
    # (Often one uses sigma_z for dephasing; adjust if needed.)
    individual_dephase_ops = [ sz(i) for i in 1:N ]
    
    # Build lists for all collapse operators and their rates.
    op_list = []
    rate_list = []
    
    # Append collective decay operator:
    push!(op_list, S_minus)
    push!(rate_list, Γ_c)
    
    # Append individual decay operators:
    for op in individual_decay_ops
        push!(op_list, op)
        push!(rate_list, γ_l)
    end
    
    # Append individual pumping operators:
    for op in individual_pump_ops
        push!(op_list, op)
        push!(rate_list, κ_l)
    end
    
    # Append individual dephasing operators:
    for op in individual_dephase_ops
        push!(op_list, op)
        push!(rate_list, d_l)
    end
    
    return b_tensor, op_list, rate_list
end

function fully_excited_state_sum(N::Int)
    subs = dicke_subspaces(N)
    sum_space = directsum(subs...)
    psi_max = spinup(subs[1])
    psi_rest = [Ket(b, zeros(ComplexF64, length(b))) for b in subs[2:end]]
    ψ = directsum(psi_max, psi_rest...)
    return ψ
end

fully_excited_state_product(N::Int) = tensor([spinup(SpinBasis(1//2)) for i=1:N]...)

function fully_ground_state_sum(N::Int)
    subs = dicke_subspaces(N)
    sum_space = directsum(subs...)
    psi_max = spindown(subs[1])
    psi_rest = [Ket(b, zeros(ComplexF64, length(b))) for b in subs[2:end]]
    ψ = directsum(psi_max, psi_rest...)
    return ψ
end

fully_ground_state_product(N::Int) = tensor([spindown(SpinBasis(1//2)) for i=1:N]...)

function TavisCummingsModel(N, γ, r, ξ, κ, g, Ωcav, Ωspin; n_phot = 3)
    nb = n_phot
    bcav = FockBasis(nb)
    id_cav = identityoperator(bcav)
    sum_space, sum_op_list, sum_rate_list = build_collapse_operators(N, 0.0, γ, r, ξ)
    id_spin = identityoperator(sum_space)
    ad = id_spin ⊗ create(bcav)
    a = id_spin ⊗ destroy(bcav)

    sum_op_list = [op ⊗ id_cav for op in sum_op_list]
    push!(sum_op_list, a)
    push!(sum_rate_list, κ)

    Hint = g * (sigmap(sum_space) ⊗ destroy(bcav) + sigmam(sum_space) ⊗ create(bcav))
    Hdrive_cav = Ωcav * (identityoperator(sum_space) ⊗ destroy(bcav) + identityoperator(sum_space) ⊗ create(bcav))
    Hdrive_spin = Ωspin * (sigmap(sum_space) ⊗ identityoperator(bcav) + sigmam(sum_space) ⊗ identityoperator(bcav))
    Htot = Hint + Hdrive_spin + Hdrive_cav
    psi0_sum = fully_ground_state_sum(N)

    Sp_sum = sigmap(sum_space) ⊗ identityoperator(bcav)
    Sm_sum = sigmam(sum_space) ⊗ identityoperator(bcav)
    Sz_sum = sigmaz(sum_space) ⊗ identityoperator(bcav)
    return sum_space, bcav, psi0_sum ⊗ coherentstate(bcav, 0.0), Htot, sum_op_list, sum_rate_list, Sp_sum, Sm_sum, Sz_sum
end

using QuantumCumulants
using ModelingToolkit
using OrdinaryDiffEq

function TavisCummingsMeanfield(;
    order = 2,
    scaled = true,
    mix_choice = maximum
)
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha

    @cnumbers N γ r ξ κ g Ωcav Ωspin
    @qnumbers a::Destroy(h)

    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)

    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)
    l = Index(h, :l, N, ha)

    sm(k) = σ(1, 2, k)
    sp(k) = σ(2, 1, k)
    see(k) = σ(2, 2, k)

    H_int = Σ(g * (a' * sm(i) + a * sp(i)), i)
    H_drive_cav = Ωcav * (a' + a)
    H_drive_spin = Σ(Ωspin * (sp(i) + sm(i)), i)

    H = H_int + H_drive_cav + H_drive_spin

    J = [a, sm(i), sp(i), see(i) ]
    rates = [ κ, γ, r, ξ]
    ops = [a, sm(j),see(j)]

    eqs = meanfield(ops, H, J; rates = rates, order = order, mix_choice = mix_choice)
    eqs_c = complete(eqs)

    # Replaces sums by N * whatever
    eqs_sc = scaled ? scale(eqs_c) : eqs_c

    @named sys = System(eqs_sc)
    return (; sys, eqs = eqs_sc, eqs_c, h, a, σ, sm, sp, see, i, j, params = (; N, γ, r, ξ, κ, g, Ωcav, Ωspin))
end

function TavisCummingsMeanfieldModel(;
    order = 2,
    mix_choice = maximum,
    scaled = true,
    N = 1000,
    γ = 0.1,
    r = 4.0,
    ξ = 0.0,
    κ = 0.5,
    g = 0.5 / sqrt(2),
    Ωcav = 0.0,
    Ωspin = 1e-3,
)
    mf = TavisCummingsMeanfield(order = order, scaled = scaled, mix_choice = mix_choice)

    sys = mf.sys
    p = mf.params

    u0 = Dict(unknowns(sys) .=> zeros(ComplexF64, length(unknowns(sys))))

    p0 = Dict(
        p.N => N,
        p.γ => γ,
        p.r => r,
        p.ξ => ξ,
        p.κ => κ,
        p.g => g,
        p.Ωcav => Ωcav,
        p.Ωspin => Ωspin,
    )

    # Common observable handles.
    a = mf.a
    sm1 = mf.sm(1)
    sp1 = mf.sp(1)
    sm2 = mf.sm(2)
    sp2 = mf.sp(2)
    see1 = mf.see(1)

    # Collective observables for scaled cumulant equations.
    # These are not new symbolic variables, just convenient reconstruction rules.
    Sminus(sol) = N .* sol[sm1]
    Splus(sol) = N .* sol[sp1]
    Nee(sol) = N .* sol[see1]

    return mf, sys, u0, p0, (; a, sm1, sp1, sm2, sp2, see1, Sminus, Splus, Nee)
end

# ============================================================
# Parameter wrapper / convenience constructors
# ============================================================

Base.@kwdef struct TCParams
    N::Int = 20

    γ::Float64 = 0.1
    r::Float64 = 4.0
    ξ::Float64 = 0.0

    κ::Float64 = 0.5
    g::Float64 = 0.5 / sqrt(2)

    Ωcav::Float64 = 0.0
    Ωspin::Float64 = 1e-6

    n_phot::Int = 6
end

function TavisCummingsModel(p::TCParams)
    return TavisCummingsModel(
        p.N,
        p.γ,
        p.r,
        p.ξ,
        p.κ,
        p.g,
        p.Ωcav,
        p.Ωspin;
        n_phot = p.n_phot,
    )
end

function TavisCummingsMeanfieldModel(
    p::TCParams;
    order = 2,
    scaled = true,
    mix_choice = maximum,
)
    return TavisCummingsMeanfieldModel(;
        order = order,
        scaled = scaled,
        mix_choice = mix_choice,
        N = p.N,
        γ = p.γ,
        r = p.r,
        ξ = p.ξ,
        κ = p.κ,
        g = p.g,
        Ωcav = p.Ωcav,
        Ωspin = p.Ωspin,
    )
end

with_spin_drive(p::TCParams, Ωspin::Real) = TCParams(; 
    N = p.N, γ = p.γ, r = p.r, ξ = p.ξ,
    κ = p.κ, g = p.g, Ωcav = p.Ωcav,
    Ωspin = Float64(Ωspin), n_phot = p.n_phot,
)

with_cavity_drive(p::TCParams, Ωcav::Real) = TCParams(; 
    N = p.N, γ = p.γ, r = p.r, ξ = p.ξ,
    κ = p.κ, g = p.g, Ωcav = Float64(Ωcav),
    Ωspin = p.Ωspin, n_phot = p.n_phot,
)

without_drives(p::TCParams) = TCParams(; 
    N = p.N, γ = p.γ, r = p.r, ξ = p.ξ,
    κ = p.κ, g = p.g, Ωcav = 0.0,
    Ωspin = 0.0, n_phot = p.n_phot,
)

