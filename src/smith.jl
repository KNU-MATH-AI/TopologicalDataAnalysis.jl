function gauss!(M::Matrix{T}) where T <: Integer
    Q = M[1,2:end] .// M[1,1]
    M[:,2:end] = (denominator.(Q)' .* M[:,2:end]) - (numerator.(Q)' .* M[:,1])
    return M
end

function rswap!(M::Matrix{T}, i1::T,i2::T) where T <: Integer
    M[[i1,i2],:] = M[[i2,i1],:]
    return M
end

# TODO: If m << n or n << m, there is no reason let M given size
function smith(M::Matrix{T}) where T <: Integer
    D_ = Int64[]
    m,n = size(M)
    m > n && M = M'
    return smith!(D_, M)
end

# TODO: We don't know what happen in case of M has nullity
function smith!(D_::Array{T}, M::Matrix{T}) where T <: Integer
    print('.')

    m = size(M)[1] # row size
    if m == 1 push!(D_, M[1,1]); return D_ end

    # TODO: We only have row swap operation
    for i in 2:m
        if (M[1,1] == 0)
            rswap!(M, 1,i)
            continue
        elseif (M[1,1]*M[2,2] == M[1,2]*M[2,1])
            rswap!(M, 2,i)
            continue
        end
        break
    end
    push!(D_, M[1,1])
    gauss!(M)
    return smith!(D_, M[2:end, 2:end])
    # @error "There is a singularity!"
end