using OrdinaryDiffEq

# --------------------------
# Mean-field helper quantities
# --------------------------

cavity_cooperativity(g, κ, γ, γp) = g^2 / (κ * (γ + γp))

function pump_thresholds(N, g, κ, γ)
    x = g^2 * N / κ
    disc = x * (x / 4 - 2γ)
    if disc < 0
        return (NaN, NaN)
    end
    root = sqrt(disc)
    center = x / 2 - γ
    return (center - root, center + root)
end

steady_state_sz(g, κ, γ, γp) = 1 / (2 * cavity_cooperativity(g, κ, γ, γp))

function steady_state_alpha2(N, g, κ, γ, γp)
    Cp = cavity_cooperativity(g, κ, γ, γp)
    return (1 / κ) * (N * (γp - γ) / 2 - (γ + γp) / (2 * Cp))
end

function steady_state_s2(N, g, κ, γ, γp)
    α2 = steady_state_alpha2(N, g, κ, γ, γp)
    return (κ^2 / g^2) * α2
end

function has_lasing_steady_state(N, g, κ, γ, γp)
    α2 = steady_state_alpha2(N, g, κ, γ, γp)
    return isfinite(α2) && α2 > 0
end

# --------------------------
# Mean-field ODEs
# State vector u = [Re(α), Im(α), Re(s), Im(s), sz]
# --------------------------

function mf_superradiant_laser!(du, u, p, t)
    N, g, κ, γ, γp = p

    α = u[1] + 1im * u[2]
    s = u[3] + 1im * u[4]
    sz = u[5]

    dα = -κ * α - 1im * g * s
    ds = -(γ + γp) * s + 2im * g * α * sz
    dsz = -2 * (γ + γp) * sz + 1im * g * (conj(α) * s - α * conj(s)) + N * (γp - γ)

    du[1] = real(dα)
    du[2] = imag(dα)
    du[3] = real(ds)
    du[4] = imag(ds)
    du[5] = real(dsz)
    return nothing
end

# --------------------------
# Convenience wrappers
# --------------------------

function solve_mf_superradiant_laser(N, g, κ, γ, γp;
    α0 = 0.0 + 0.0im,
    s0 = 0.0 + 0.0im,
    sz0 = -N / 2,
    tspan = (0.0, 100.0),
    saveat = nothing,
    solver = Tsit5()
)
    u0 = [real(α0), imag(α0), real(s0), imag(s0), sz0]
    p = (N, g, κ, γ, γp)
    prob = ODEProblem(mf_superradiant_laser!, u0, tspan, p)
    return solve(prob, solver; saveat = saveat, reltol = 1e-9, abstol = 1e-9)
end

# Extract useful observables from a solution
alpha_of_u(u) = u[1] + 1im * u[2]
s_of_u(u) = u[3] + 1im * u[4]
sz_of_u(u) = u[5]

photon_number(u) = abs2(alpha_of_u(u))
dipole_coherence(u) = abs2(s_of_u(u))
inversion(u) = sz_of_u(u)

# helper: actual lasing-window check
function in_lasing_regime(N, g, κ, γ, γp)
    γp_minus, γp_plus = pump_thresholds(N, g, κ, γ)
    return isfinite(γp_minus) && isfinite(γp_plus) && (γp_minus < γp < γp_plus)
end