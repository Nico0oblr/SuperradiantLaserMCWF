# Square root rates

# Collective factors (used also in the individual channels):

A_JM_minus(J, M) = sqrt((J + M) * (J - M + 1))
A_JM_plus(J, M) = sqrt((J - M) * (J + M + 1))

# ---------------------------
# Individual Decay (lowering M)
# ---------------------------

function P_JM_minus_0(J, M, N)
    if J == 0
        return 0.0
    end
    return sqrt((2 + N) / (4 * J * (J + 1))) * A_JM_minus(J, M)
end

function P_JM_minus_minus(J, M, N) 
    if J == 0
        return 0.0
    end
    return -sqrt((N + 2 * J + 2) * (J + M) * (J + M - 1) / (4 * J * (2 * J + 1)))
end

P_JM_minus_plus(J, M, N) = sqrt((N - 2 * J) * (J - M + 1) * (J - M + 2) / (4 * (J + 1) * (2 * J + 1)))

# ---------------------------
# Individual Pumping (raising M)
# ---------------------------

function P_JM_plus_0(J, M, N)
    if J == 0
        return 0.0
    end
    return sqrt((2 + N) / (4 * J * (J + 1))) * A_JM_plus(J, M)
end

function P_JM_plus_minus(J, M, N)
    if J == 0
        return 0.0
    end
    return sqrt((N + 2 * J + 2) * (J - M) * (J - M - 1) / (4 * J * (2 * J + 1)))
end

P_JM_plus_plus(J, M, N) = -sqrt((N - 2 * J) * (J + M + 1) * (J + M + 2) / (4 * (J + 1) * (2 * J + 1)))

# ---------------------------
# Individual Dephasing
# ---------------------------

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
function P_JM_z_minus(J, M, N)
    if J == 0
        return 0.0
    end
    sqrt((N + 2 * J + 2) * (J - M) * (J + M) / (4 * J * (2 * J + 1)))
end

"""
    P_JM_z_plus(J, M, N)

Individual dephasing channel increasing J by 1 (“s = +”):
  P_{JM}^{z,+} = √[(N-2J)(J+1-M)(J+1+M)/(4(J+1)(2J+1))].
"""
P_JM_z_plus(J, M, N) = sqrt((N - 2 * J) * (J + 1 - M) * (J + 1 + M) / (4 * (J + 1) * (2 * J + 1)))

# Real rates

function P_JM_z_minus2(J, M, N)
    if J == 0
        return 0.0
    end
    return (N + 2 * J + 2) * (J - M) * (J + M) / (4 * J * (2 * J + 1))
end

A_JM_minus2(J, M) = (J + M) * (J - M + 1)

P_JM_z_plus2(J, M, N) =
    (N - 2 * J) * (J + 1 - M) * (J + 1 + M) / (4 * (J + 1) * (2 * J + 1))

function P_JM_minus_02(J, M, N)
    if J == 0
        return 0.0
    end
    return ((2 + N) / (4 * J * (J + 1))) * A_JM_minus2(J, M)
end

function P_JM_minus_minus2(J, M, N)
    if J == 0
        return 0.0
    end
    return (N + 2 * J + 2) * (J + M) * (J + M - 1) / (4 * J * (2 * J + 1))
end

P_JM_minus_plus2(J, M, N) = (N - 2 * J) * (J - M + 1) * (J - M + 2) / (4 * (J + 1) * (2 * J + 1))

#####################

A_JM_minus2(J, M) = (J + M) * (J - M + 1)
A_JM_plus2(J, M)  = (J - M) * (J + M + 1)

function P_JM_z_minus2(J, M, N)
    if J == 0
        return 0.0
    end
    return (N + 2 * J + 2) * (J - M) * (J + M) / (4 * J * (2 * J + 1))
end

P_JM_z_plus2(J, M, N) =
    (N - 2 * J) * (J + 1 - M) * (J + 1 + M) / (4 * (J + 1) * (2 * J + 1))

function P_JM_minus_02(J, M, N)
    if J == 0
        return 0.0
    end
    return ((2 + N) / (4 * J * (J + 1))) * A_JM_minus2(J, M)
end

function P_JM_minus_minus2(J, M, N)
    if J == 0
        return 0.0
    end
    return (N + 2 * J + 2) * (J + M) * (J + M - 1) / (4 * J * (2 * J + 1))
end

P_JM_minus_plus2(J, M, N) =
    (N - 2 * J) * (J - M + 1) * (J - M + 2) / (4 * (J + 1) * (2 * J + 1))

function P_JM_plus_02(J, M, N)
    if J == 0
        return 0.0
    end
    return ((2 + N) / (4 * J * (J + 1))) * A_JM_plus2(J, M)
end

function P_JM_plus_minus2(J, M, N)
    if J == 0
        return 0.0
    end
    return (N + 2 * J + 2) * (J - M) * (J - M - 1) / (4 * J * (2 * J + 1))
end

P_JM_plus_plus2(J, M, N) = (N - 2 * J) * (J + M + 1) * (J + M + 2) / (4 * (J + 1) * (2 * J + 1))

function P_JM_z_02(J, M, N)
    if J == 0
        return 0.0
    else
        return (2 + N) / (4 * J * (J + 1)) * M ^ 2
    end
end

######### Model helpers

abstract type AbstractSpinModel end

# global decay: (S,M) -> (S,M-1)
@inline poprate_gdec(S::Int, M::Int, model::AbstractSpinModel) = model.global_decay * A_JM_minus2(S, M)
# global pump: (S,M) -> (S,M+1)
@inline poprate_gpump(S::Int, M::Int, model::AbstractSpinModel) = model.global_pump * A_JM_plus2(S, M)
# local decay branches
@inline poprate_ldec_plus(S::Int, M::Int, model::AbstractSpinModel) = model.local_decay * P_JM_minus_plus2(S, M, model.N)
@inline poprate_ldec_0(S::Int, M::Int, model::AbstractSpinModel) = model.local_decay * P_JM_minus_02(S, M, model.N)
@inline poprate_ldec_minus(S::Int, M::Int, model::AbstractSpinModel) = model.local_decay * P_JM_minus_minus2(S, M, model.N)
# local pump branches
@inline poprate_lpump_plus(S::Int, M::Int, model::AbstractSpinModel) = model.local_pump * P_JM_plus_plus2(S, M, model.N)
@inline poprate_lpump_0(S::Int, M::Int, model::AbstractSpinModel) = model.local_pump * P_JM_plus_02(S, M, model.N)
@inline poprate_lpump_minus(S::Int, M::Int, model::AbstractSpinModel) = model.local_pump * P_JM_plus_minus2(S, M, model.N)
# local dephasing branches
@inline poprate_deph_plus(S::Int, M::Int, model::AbstractSpinModel) = model.local_dephasing * P_JM_z_plus2(S, M, model.N)
@inline poprate_deph_0(S::Int, M::Int, model::AbstractSpinModel) = model.local_dephasing * P_JM_z_02(S, M, model.N)
@inline poprate_deph_minus(S::Int, M::Int, model::AbstractSpinModel) = model.local_dephasing * P_JM_z_minus2(S, M, model.N)