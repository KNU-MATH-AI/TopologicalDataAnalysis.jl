function det22(M::Matrix{T}) where T <: Integer
    return (M[1,1]*M[2,2] == M[1,2]*M[2,1])
end
function gauss!(M::Matrix{T}) where T <: Integer
    Q = M[1,2:end] .// first(M)
    M[:,2:end] = (denominator.(Q)' .* M[:,2:end]) - (numerator.(Q)' .* M[:,1])
    return M
end

function cswap!(M::Matrix{T}, j1::T,j2::T) where T <: Integer
    M[:,[j1,j2]] = M[:,[j2,j1]]
    return M
end
function rswap!(M::Matrix{T}, i1::T,i2::T) where T <: Integer
    M[[i1,i2],:] = M[[i2,i1],:]
    return M
end

function ccycle!(M::Matrix{T}) where T <: Integer
    M = hcat(M[:,2:end], M[:,1])
    return M
end
function rcycle!(M::Matrix{T}) where T <: Integer
    M = vcat(M[2:end,:], M[1,:])
    return M
end


# TODO: If m << n or n << m, there is no reason let M given size
function smith(M::Matrix{T}) where T <: Integer
    D_ = Int64[]
    m,n = size(M)
    if (m > n)
        M = M'
        m,n = n,m
    else
        M = copy(M)
    end
    return smith!(D_, M)
end

# TODO: We don't know what happen in case of M has nullity
function smith!(D_::Array{T}, M::Matrix{T}) where T <: Integer
    m,n = size(M)
    @assert m ≤ n

    if iszero(M)
        return append!(D_, zeros(m))
    end
    if m == 1
        push!(D_, first(M))
        return D_
    end
    
    if first(M) == 0
        if iszero(M[:,1])
            ccycle!(M)
            return smith!(D_, M)
        end
        if iszero(M[1,:])
            rcycle!(M)
            return smith!(D_, M)
        end
        # ↑ Efficient Implementation

        # ↓ Division Error Handling
        #   I geuss for-loof is more easy to read.
        for j in 2:n
            if !iszero(M[1,j])
                cswap!(M, 1,j)
                break
            end
        end
    end

    for j in 2:n
        if det22(M)
            cswap!(M, 2,j)
            continue
        else
            break
        end
    end

    push!(D_, first(M))
    gauss!(M)
    M = M[2:end, 2:end]
    return smith!(D_, M)
end