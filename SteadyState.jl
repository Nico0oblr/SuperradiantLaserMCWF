module SteadyState

using SparseArrays
using LinearAlgebra
using DataStructures
using OrdinaryDiffEq
using PyPlot

include("SuperradiantLasing/SuperradiantLasing.jl")
using .SuperradiantLasing


export active_basis_from_samples,
       pad_basis_by_outgoing_shells,
       set_up_rate_matrix_on_basis,
       population_from_samples,
       steady_by_implicit_euler,
       evolve_rate_ode,
       plot_population_heatmap,
       plot_population_scatter,
       stationary_diagnostics,
       dicke_counts_sparse,
       plot_dicke_counts_sparse

# ------------------------------------------------------------
# Active basis construction
# ------------------------------------------------------------

function active_basis_from_samples(S_samples, M_samples)
    return unique(Tuple.(zip(round.(Int, S_samples), round.(Int, M_samples))))
end

function pad_basis_by_outgoing_shells(active_basis, params, Jmax; nshells=1)
    Ω = Set(active_basis)

    for _ in 1:nshells
        new_states = Set{eltype(Ω)}()

        for state in Ω
            out = outgoing_for_state(state, params, Jmax)
            for dst in keys(out)
                push!(new_states, dst)
            end
        end

        union!(Ω, new_states)
    end

    return collect(Ω)
end

# ------------------------------------------------------------
# Restricted rate matrix
# Assumes outgoing_for_state(state, params, Jmax) is defined externally.
# Convention: dp/dt = R p
# ------------------------------------------------------------

function set_up_rate_matrix_on_basis(active_basis, params, Jmax; close_boundary=true)
    basis = collect(active_basis)
    index_mapping = Dict(basis .=> eachindex(basis))

    sparse_constructor = DefaultDict{Tuple{Int, Int}, Float64}(0.0)
    leakage = zeros(Float64, length(basis))

    for src in basis
        js = index_mapping[src]
        out = outgoing_for_state(src, params, Jmax)

        internal_out_rate = 0.0
        total_out_rate = 0.0

        for (dst, amp) in out
            rate = amp^2
            total_out_rate += rate

            if haskey(index_mapping, dst)
                is = index_mapping[dst]
                sparse_constructor[(is, js)] += rate
                internal_out_rate += rate
            else
                leakage[js] += rate
            end
        end

        diag_rate = close_boundary ? internal_out_rate : total_out_rate
        sparse_constructor[(js, js)] -= diag_rate
    end

    I = Int[]
    J = Int[]
    V = Float64[]

    sizehint!(I, length(sparse_constructor))
    sizehint!(J, length(sparse_constructor))
    sizehint!(V, length(sparse_constructor))

    for ((i, j), v) in sparse_constructor
        push!(I, i)
        push!(J, j)
        push!(V, v)
    end

    n = length(basis)
    R = sparse(I, J, V, n, n)
    dropzeros!(R)

    return basis, index_mapping, R, leakage
end

# ------------------------------------------------------------
# Initial distribution from samples
# ------------------------------------------------------------

function population_from_samples(basis, idx, S_samples, M_samples)
    p0 = zeros(Float64, length(basis))

    for state in Tuple.(zip(round.(Int, S_samples), round.(Int, M_samples)))
        if haskey(idx, state)
            p0[idx[state]] += 1.0
        end
    end

    if sum(p0) == 0
        error("No samples lie inside the provided basis.")
    end

    p0 ./= sum(p0)
    return p0
end

# ------------------------------------------------------------
# Pseudo-time steady state by implicit Euler
# ------------------------------------------------------------

function steady_by_implicit_euler(R, p0;
    dt = 1e4,
    nsteps = 5,
    renormalize = true,
    clip_negative = true,
)
    n = length(p0)
    A = sparse(I, n, n) - dt * R
    F = factorize(A)

    p = copy(p0)

    for _ in 1:nsteps
        p = F \ p

        if clip_negative
            p[p .< 0] .= 0.0
        end

        if renormalize
            p ./= sum(p)
        end
    end

    return p
end

# ------------------------------------------------------------
# ODE fallback
# ------------------------------------------------------------

function evolve_rate_ode(R, p0, tspan;
    alg = TRBDF2(),
    reltol = 1e-8,
    abstol = 1e-12,
    saveat = nothing,
)
    function rhs!(du, u, p, t)
        mul!(du, R, u)
        return nothing
    end

    prob = ODEProblem(rhs!, p0, tspan)

    sol = solve(
        prob,
        alg;
        reltol = reltol,
        abstol = abstol,
        saveat = saveat,
        save_everystep = saveat === nothing,
    )

    return sol
end

# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------

function stationary_diagnostics(R, p; leakage=nothing)
    out = Dict{String, Float64}()

    out["sum"] = sum(p)
    out["min"] = minimum(p)
    out["max"] = maximum(p)
    out["residual_l1"] = norm(R * p, 1)
    out["residual_linf"] = norm(R * p, Inf)
    out["columnsum_linf"] = norm(vec(sum(R, dims=1)), Inf)

    if leakage !== nothing
        out["mean_leakage"] = dot(leakage, p)
        out["max_leakage"] = maximum(leakage)
    end

    return out
end

# ------------------------------------------------------------
# Plotting: population heatmap on basis
# ------------------------------------------------------------

function plot_population_heatmap(
    basis,
    p;
    logscale = true,
    eps = 1e-30,
    vmin = nothing,
    vmax = nothing,
    cmap = "magma",
    outfile = nothing,
)
    Svals = first.(basis)
    Mvals = last.(basis)

    Smin, Smax = minimum(Svals), maximum(Svals)
    Mmin, Mmax = minimum(Mvals), maximum(Mvals)

    grid = fill(NaN, Smax - Smin + 1, Mmax - Mmin + 1)

    for (i, (S, M)) in enumerate(basis)
        val = logscale ? log10(p[i] + eps) : p[i]
        grid[S - Smin + 1, M - Mmin + 1] = val
    end

    fig, ax = subplots(figsize=(4.5, 4.0))

    kwargs = Dict{Symbol, Any}(
        :origin => "lower",
        :aspect => "auto",
        :interpolation => "nearest",
        :cmap => cmap,
    )

    if vmin !== nothing
        kwargs[:vmin] = vmin
    end
    if vmax !== nothing
        kwargs[:vmax] = vmax
    end

    im = ax.imshow(grid'; kwargs...)

    ax.set_xlabel("S")
    ax.set_ylabel("M")

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(logscale ? L"\log_{10} p" : L"p")

    fig.tight_layout()

    if outfile !== nothing
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
    end

    return fig, ax, grid
end

function plot_population_scatter(
    basis,
    p;
    N = nothing,
    logscale = true,
    eps = 1e-30,
    markersize = 1.0,
    cmap = "magma",
    outfile = nothing,
)
    Svals = first.(basis)
    Mvals = last.(basis)

    if N === nothing
        x = Svals
        y = Mvals
        xlabel = "S"
        ylabel = "M"
    else
        x = Svals ./ N
        y = Mvals ./ N
        xlabel = L"S/N"
        ylabel = L"M/N"
    end

    cvals = logscale ? log10.(p .+ eps) : p

    fig, ax = subplots(figsize=(4.5, 4.0))

    sc = ax.scatter(
        x,
        y;
        c = cvals,
        s = markersize,
        cmap = cmap,
        marker = "s",
        linewidths = 0,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(logscale ? L"\log_{10} p" : L"p")

    fig.tight_layout()

    if outfile !== nothing
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
    end

    return fig, ax
end

# ------------------------------------------------------------
# Sparse sample-count plotting
# ------------------------------------------------------------

function dicke_counts_sparse(S_samples, M_samples; N::Int)
    Jmax = N ÷ 2

    rows = Int[]
    cols = Int[]

    sizehint!(rows, length(S_samples))
    sizehint!(cols, length(S_samples))

    for k in eachindex(S_samples, M_samples)
        S = round(Int, S_samples[k])
        M = round(Int, M_samples[k])

        if 0 <= S <= Jmax && abs(M) <= S
            push!(rows, M + Jmax + 1)
            push!(cols, S + 1)
        end
    end

    vals = ones(Int, length(rows))
    return sparse(rows, cols, vals, 2Jmax + 1, Jmax + 1)
end

function plot_dicke_counts_sparse(
    S_samples,
    M_samples;
    N::Int,
    outfile = nothing,
    logscale::Bool = true,
    markersize = 1.0,
    cmap = "magma",
)
    Jmax = N ÷ 2
    C = dicke_counts_sparse(S_samples, M_samples; N=N)

    rows, cols, vals = findnz(C)

    M = rows .- Jmax .- 1
    S = cols .- 1

    s = S ./ N
    m = M ./ N

    colorvals = logscale ? log10.(vals .+ 1.0) : Float64.(vals)
    cbar_label = logscale ? L"\log_{10}(1+\mathrm{counts})" : "counts"

    fig, ax = subplots(figsize=(4.5, 4.0))

    ss = range(0.0, 0.5; length=500)
    ax.fill_between(ss, ss, 0.5; color="0.85", linewidth=0, zorder=0)
    ax.fill_between(ss, -0.5, -ss; color="0.85", linewidth=0, zorder=0)

    sc = ax.scatter(
        s,
        m;
        c = colorvals,
        s = markersize,
        cmap = cmap,
        marker = "s",
        linewidths = 0,
        zorder = 2,
    )

    ax.plot(ss, ss; color="black", linestyle="--", linewidth=0.8)
    ax.plot(ss, -ss; color="black", linestyle="--", linewidth=0.8)
    ax.plot([0.5, 0.5], [-0.5, 0.5]; color="black", linestyle="--", linewidth=0.8)

    ax.set_xlim(0.0, 0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_aspect(0.5)
    ax.set_xlabel(L"S/N")
    ax.set_ylabel(L"M/N")

    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(cbar_label)

    fig.tight_layout()

    if outfile !== nothing
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
    end

    return fig, ax, C
end

end # module