function integer_hist(x::Vector{Int})
    xmin = minimum(x)
    xmax = maximum(x)
    counts = zeros(Int, xmax - xmin + 1)
    @inbounds for xi in x
        counts[xi - xmin + 1] += 1
    end
    xs = collect(xmin:xmax)
    return xs, counts
end

function plot_integer_step(xs, p; label = "")
    plot(xs, p; drawstyle = "steps-mid", linewidth = 1.5, label = label)
end