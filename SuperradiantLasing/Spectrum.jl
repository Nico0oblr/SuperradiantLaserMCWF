using FFTW
using LsqFit

function spectrum_from_correlator(τ_grid, Cτ, ω_grid; η=0.0)
    Δτ = τ_grid[2] - τ_grid[1]
    Sω = zeros(Float64, length(ω_grid))

    @inbounds for i in eachindex(ω_grid)
        ω = ω_grid[i]
        s = 0.0
        for j in eachindex(τ_grid)
            s += Cτ[j] * exp(-η * τ_grid[j]) * cos(ω * τ_grid[j])
        end
        Sω[i] = 2 * Δτ * s
    end
    return Sω
end

using FFTW


function spectrum_even_fft(τ_grid, Cτ; η=0.0)
    N = length(τ_grid)
    Δτ = τ_grid[2] - τ_grid[1]

    @assert all(isapprox(τ_grid[j+1] - τ_grid[j], Δτ; rtol=1e-10, atol=1e-12) for j in 1:N-1)
    @assert isapprox(τ_grid[1], 0.0; atol=1e-12)

    Cp = Cτ .* exp.(-η .* τ_grid)

    # even extension: 0,1,2,...,N-1,N-2,...,1
    Ceven = vcat(Cp, Cp[end-1:-1:2])

    M = length(Ceven)
    F = fft(Ceven) * Δτ

    # angular-frequency grid
    ω = 2π .* vcat(0:M÷2, -((M-1)÷2):-1) ./ (M * Δτ)

    p = sortperm(ω)
    ω = ω[p]
    Sω = real.(F[p])

    return ω, Sω
end

# --------------------------------
# FWHM of a single-peaked spectrum
# Assumes ω_grid is sorted and Sω >= 0
# --------------------------------
function spectrum_fwhm(ω_grid, Sω)
    y = Sω ./ maximum(Sω)
    i0 = argmax(y)
    half = 0.5

    # left crossing
    iL = nothing
    for i in i0:-1:2
        if y[i-1] <= half <= y[i]
            iL = i
            break
        end
    end

    # right crossing
    iR = nothing
    for i in i0:length(y)-1
        if y[i] >= half >= y[i+1]
            iR = i
            break
        end
    end

    if isnothing(iL) || isnothing(iR)
        return NaN
    end

    # linear interpolation
    ωL = ω_grid[iL-1] + (half - y[iL-1]) * (ω_grid[iL] - ω_grid[iL-1]) / (y[iL] - y[iL-1])
    ωR = ω_grid[iR] + (half - y[iR]) * (ω_grid[iR+1] - ω_grid[iR]) / (y[iR+1] - y[iR])

    return ωR - ωL
end

# --------------------------------
# Simple linear fit in gamma
# linewidth ≈ a + b*gamma
# --------------------------------
lin_model(x, p) = p[1] .+ p[2] .* x

function fit_linewidth_vs_gamma(gammas, widths)
    p0 = [minimum(widths), 1.0]
    fit = curve_fit(lin_model, gammas, widths, p0)
    return fit
end