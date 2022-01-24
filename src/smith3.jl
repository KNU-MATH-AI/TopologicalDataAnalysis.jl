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

function nonzero_abs(x::Real)
    if x > 0
        return x
    elseif x < 0
        return -x
    else
        return Inf
    end
end
function nonzero_argmin(M::Matrix{T}) where T <: Integer
    nonzero_M = nonzero_abs.(M)
    d = minimum(nonzero_M)
    idx = findfirst(m -> (m==d), nonzero_M)
    return idx
end
function pullfirst!(M::Matrix{T}) where T <: Integer
    idx = nonzero_argmin(M)
    rswap!(M, 1, idx[1])
    cswap!(M, 1, idx[2])
end

"""
    archimedes(a::Real,b::Real)

Return the minimul/maximul integer number such that ak ≥ b or ak ≤ b.
"""
function archimedes(a::T,b::T) where T <: Real
    if abs(a) > abs(b)
        a,b = b,a
    end
    c = b ÷ a
    if (b % a) == 0
        return c
    else
        return c + sign(c)
    end
end

"""
    issmithable(D::Array{T}, M::Matrix{T}) where T <: Real

Check the first component of given matrix M is minimul or devide all components of M.
"""
function issmithable(D::Array{T}, M::Matrix{T}) where T <: Real
    M11 = first(M)
    return (last(D) == abs(M11)) || iszero(M11) || iszero(mod.(M, M11))
end

function smith2(M::Matrix{T}) where T <: Integer
    D = [1]
    m,n = size(M)
    if (m > n)
        M = Matrix(M')
        m,n = n,m
    else
        M = copy(M)
    end
    return smith!(D, M)[2:end]
end

function smith!(D::Array{T}, M::Matrix{T}) where T <: Integer
    m,n = size(M)
    @assert m ≤ n

    if iszero(M)
        return append!(D, zeros(m))
    end
    
    pullfirst!(M)
    # if m == 1
    #     return push!(D, abs(first(M)))
    # end
    if issmithable(D, M)
        push!(D, abs(first(M)))
        gauss!(M)
        M = M[2:end, 2:end]
        return smith!(D, M)
    end

    while !issmithable(D, M)
        M11 = first(M)
        if (M[1,2] % M11) != 0
            k = archimedes(M11, M[1,2])
            M[:,2] -= k*M[:,1]
        else 
            row1 = M[1,:] .% M11
            if !iszero(row1)
                j = findlast(x -> x != 0, row1)
                cswap!(M, 2, j)
                k = archimedes(M11, M[1,2])
                M[:,2] -= k*M[:,1]
            else
                column1 = M[:,1] .% M11
                if !iszero(column1)
                    i = findlast(x -> x != 0, column1)
                    rswap!(M, 2, i)
                    k = archimedes(M11, M[2,1])
                    M[2,:] -= k*M[1,:]
                else
                    element = M .% M11
                    i,j = nonzero_argmin(element).I
                    rswap!(M,2,i)
                    cswap!(M,2,j)
                    M[1,:] += M[2,:]
                    k = archimedes(M11, M[1,2])
                    M[:,2] -= k*M[:,1]
                end
            end
        end
        pullfirst!(M)
    end

    push!(D, abs(first(M)))
    gauss!(M)
    M = M[2:end, 2:end]
    return smith!(D, M)
end

using SmithNormalForm, LinearAlgebra

M = [2 3 3 5; 3 -1 -5 2; 3 0 6 9; -2 -2 4 0]; X = copy(M)
@time smith2(X)
@time answer = X |> smith |> diagm |> diag

X = rand(0:10, 100,100)
@time smith2(X)
@time answer = X |> smith |> diagm |> diag

prod(smith2(X) .== X |> smith |> diagm |> diag)

Y = rand(0:10, 10,1000)
@time smith2(Y)
@time answer = Y |> smith |> diagm |> diag