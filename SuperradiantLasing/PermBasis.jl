using DataStructures
using SparseArrays
using LinearAlgebra


function build_basis(J)
    S0 = iseven(2*J) ? 0 : 1/2
    return [(S,M) for S=S0:J for M=-S:S]
end

function sane(J,M,Jmax)
    if J > Jmax return false end
    if abs(M) > J return false end
    if J < 0 return false end
    return true
end 

valid(s,m) = abs(m) <= s || isapprox(s, m, atol=1e-5)
mvals(N) = range(-N/2, N/2)
svals(N) = collect(0:N÷2)