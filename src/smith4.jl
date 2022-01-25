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
    if iszero(first(M))
        @error "The First Element M₁₁ is zero, $(M[1,1]) $(M[1,2]) \n $(M[2,1]) $(M[2,2])"
    end
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
    return (last(D) == abs(M11)) || iszero(mod.(M, M11))
end

function firstrow(M::Matrix{T}) where T <: Integer
    return first(eachrow(M))
end
function firstcol(M::Matrix{T}) where T <: Integer
    return first(eachcol(M))
end

function ccycle!(M::Matrix{T}) where T <: Integer
    @views M = hcat(M[:,2:end], M[:,[1]])
    return M
end
function rcycle!(M::Matrix{T}) where T <: Integer
    @views M = vcat(M[2:end,:], M[[1],:])
    return M
end

function pull!(M::Matrix{T}) where T <: Integer
    return nothing
end


function smith4(M::Matrix{T}) where T <: Integer
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
    # Now M has nonzero element at least one.

    if m == 1
        return push!(D, abs(M[nonzero_argmin(M)]))
    end

    while !issmithable(D, M)

        while iszero(first(M))
            # println("1")
            j = findfirst(!iszero, firstrow(M))
            # j = nonzero_argmin(firstrow(M))
            if isnothing(j) # It means the firstrow of M is a zero vector.
                if iszero(firstcol(M))
                    ccycle!(M)
                end
                rcycle!(M)
            else
                cswap(M, 1,j)
            end
        end
        M11 = first(M)
        @assert !iszero(first(M)) # Now the M₁₁ is nonzero.

        while mod(M[1,2], first(M)) == 0
            # println(M)
            # println("?")
            j = findfirst(M1j -> (mod(M1j, first(M)) != 0), firstrow(M))
            if isnothing(j)
                if m > 2
                    gauss!(M)
                    M[2,1] = first(M)
                    rcycle!(M)
                else
                    M[[1,2],:] = M[[2,1],:]
                end
            else
                cswap!(M, 2,j)
            end
        end
        # @assert mod(M[1,2], M11) != 0

        if abs(first(M)) > abs(M[1,2])
            cswap!(M, 1,2)
            continue
        end

        k = (M[1,2] ÷ first(M))
        M[:,2] -= k*M[:,1]
        
        if abs(first(M)) > abs(M[1,2])
            cswap!(M, 1,2)
        end
        # println(M)
    end

    push!(D, abs(first(M)))
    gauss!(M)
    M = M[2:end, 2:end]
    return smith!(D, M)
end

using SmithNormalForm, LinearAlgebra

M = [2 3 3 5; 3 -1 -5 2; 3 0 6 9; -2 -2 4 0]; X = copy(M)
@time smith4(X)
@time answer = X |> smith |> diagm |> diag

X = rand(0:10, 100,100)
@time smith4(X)
@time answer = X |> smith |> diagm |> diag

X = rand(0:10, 100,100)
smith4(X) == (answer = X |> smith |> diagm |> diag)