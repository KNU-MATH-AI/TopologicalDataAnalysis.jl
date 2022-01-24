function cswap!(M::Matrix{T}, j1::T,j2::T) where T <: Integer
    j1 == j2 && return M
    M[:,[j1,j2]] = M[:,[j2,j1]]
    return M
end
function rswap!(M::Matrix{T}, i1::T,i2::T) where T <: Integer
    i1 == i2 && return M
    M[[i1,i2],:] = M[[i2,i1],:]
    return M
end

function gauss!(M::Matrix{T}) where T <: Integer
    Q1R = M[1,2:end] .// first(M)
    M[:,2:end] = (denominator.(Q1R)' .* M[:,2:end]) - (numerator.(Q1R)' .* M[:,1])
    Q1C = M[2:end,1] .// first(M)
    M[2:end,:] = (M[2:end,:] .* denominator.(Q1C)) - (M[1,:]' .* numerator.(Q1C))
    return M
end

# function pullfirst!(M::Matrix{T}, d::T) where T <: Integer
#     idx = findfirst(m -> abs(m)==d, M)
#     rswap!(M, 1, idx[1])
#     cswap!(M, 1, idx[2])
# end

function nonzeroabs(x::Real)
    if x > 0
        return x
    elseif x < 0
        return -x
    else
        return Inf
    end
end

function pullfirst!(M::Matrix{T}) where T <: Integer
    nonzeroM = nonzeroabs.(M)
    d = minimum(nonzeroM)
    idx = findfirst(m -> (m==d), nonzeroM)
    rswap!(M, 1, idx[1])
    cswap!(M, 1, idx[2])
end

"""
    archimedes(a::Real,b::Real)

Return the minimul/maximul integer number such that ak > b or ak < b.
"""
function archimedes(a::T,b::T) where T <: Real
    if abs(a) > abs(b)
        a,b = b,a
    end
    c = b ÷ a
    return c + sign(c)
end

function issmithable(D, M)
    M11 = first(M)
    return last(D) == abs(M11) || iszero(mod.(M, M11))
end

# ---

M = [2 3 1 5; 3 -1 -5 2; 3 0 6 9; -2 -2 4 0]; X = copy(M)
D = Int64[]

m,n = size(M)
pullfirst!(M)
push!(D, abs(first(M)))
gauss!(M)
M = M[2:end, 2:end]

#---

m,n = size(M)
pullfirst!(M)
if mod(M[1,2],M[1,1]) == 0
    M = Matrix(M')
end
k = archimedes(-M[1,2],M[1,1])
M[:,2] += k*M[:,1]
pullfirst!(M)

if last(D) == abs(first(M))
    push!(D, abs(first(M)))
    gauss!(M)
    M = M[2:end, 2:end]
end

# ---

m,n = size(M)
pullfirst!(M)
if mod(M[1,2],M[1,1]) == 0
    M = Matrix(M')
end
k = archimedes(-M[1,2],M[1,1])
M[:,2] += k*M[:,1]
pullfirst!(M)
if issmithable(D, M)
    push!(D, abs(first(M)))
    gauss!(M)
    M = M[2:end, 2:end]
end

# ---

m,n = size(M)
if m == 1
    pullfirst!(M)
    push!(D, abs(first(M)))
end


# ---
using SmithNormalForm, LinearAlgebra
@time answer = smith(X)
answer |> diagm |> diag