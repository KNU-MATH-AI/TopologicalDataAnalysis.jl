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
function pullfirst!(M::Matrix{T}, idx) where T <: Integer
    rswap!(M, 1, idx[1])
    cswap!(M, 1, idx[2])
end

function eratosthenes(n::T) where T <: Integer
    if n == 1 return [[1]], [] end
    if n == 2 return [[1], [2]], [2] end
    primes = [2]
    factorized = [[1], [2]]
    for k ∈ 3:n
        m = k
        for p ∈ primes
            if m % p == 0
                m ÷= p
                temp = [p; factorized[m]]
                push!(factorized, temp)
                break
            end
        end
        if length(factorized) != k
            push!(primes, k)
            push!(factorized, [k])
        end
    end
    return factorized, primes
end

function eratosthenes!(factorized::Vector{Vector{T}}, primes::Vector{T}, n::T) where T <: Integer
    old_n = length(factorized)
    if n ≤ old_n return factorized, primes end
    for k ∈ (old_n + 1):n
        m = k
        for p ∈ primes
            if m % p == 0
                m ÷= p
                temp = [p; factorized[m]]
                push!(factorized, temp)
                break
            end
        end
        if length(factorized) != k
            push!(primes, k)
            push!(factorized, [k])
        end
    end
    return factorized, primes
end

F, P = eratosthenes(10)
eratosthenes!(F, P, 100)
static_nfactor = length.(F)

function nfactor(k::T) where T <: Integer
    if k < 0 k = abs(k) end
    if k == 0 return typemax(T) end
    if k == 1 return 0 end
    return static_nfactor[k]
end

M = [
   -6 3 2 2 5
    2 2 0 2 7
   -2 5 2 4 3
   -7 3 7 8 9
]

idx = argmin(nfactor.(M))
pullfirst!(M, idx)

nfactor.(M[begin,:])
nfactor.(M[:,begin])