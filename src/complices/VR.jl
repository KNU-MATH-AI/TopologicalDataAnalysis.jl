using DataFrames
using SparseArrays
using NearestNeighbors

function diam(M_dist, vertices)
    submatrix = M_dist[vertices, vertices]
    if (submatrix .|> iszero |> sum) > length(vertices)
        return Inf
    else
        return maximum(submatrix)
    end
end

function vietoris_rips(point_colud::Matrix; max_dims = 3, max_epsilon = 2.0)
    # Making Distance Matrix
    d, n = size(point_colud)
    kdtree = KDTree(point_colud)
    vrtx, dist = knn(kdtree, point_colud, n, true)
    
    # This spart matrix is much more efficient when max_epsilon is taken as small value.
    M_dist = spzeros(Float16, n,n)
    for i in 1:n
        distance = dist[i]
        vertices = vrtx[i][distance .≤ max_epsilon]
        for (idxj, j) in enumerate(vertices)
            M_dist[i,j] = distance[idxj]
        end
    end

    # Initialization
    filtered_complex = DataFrame(appearance = Float16[], simplex = Array{Int8, 1}[])
    Δ = Dict()
    hash_cache = []

    for k ∈ 0:max_dims
        push!(Δ, k::Int64 => Array{Array{Int8, 1}, 1}())
        if k == 0
            for j ∈ 1:n
                splx = [j]
                push!(Δ[k], splx)
                push!(filtered_complex, [0, splx])
            end
            continue
        end
        for face in Δ[k-1]
            for v ∈ setdiff(1:maximum(face), face)
                candy = sort([v; face])
                hash_candy = hash(candy)
                if !(hash_candy ∈ hash_cache)
                    push!(hash_cache, hash_candy)
                    δ = diam(M_dist, candy)
                    if 0 < δ ≤ max_epsilon
                        push!(Δ[k], candy)
                        push!(filtered_complex, [δ, candy])
                    end
                end
            end
        end
    end
    sort!(filtered_complex, :appearance)

    appearance = filtered_complex.appearance
    δ1 = unique(appearance)

    degree = zeros(Int64, length(appearance))
    for (deg, δ) in enumerate(δ1)
        degree[appearance .== δ] .= deg-1
    end
    filtered_complex[!,:degree] = degree

    # Δ = sort(Δ) # actually we don't need this
    return filtered_complex[:,[:appearance, :degree, :simplex]]
end
VR = vietoris_rips