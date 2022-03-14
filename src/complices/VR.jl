include("abstraction.jl")
include("../utilities/tDATAsets.jl")
using DataFrames
using SparseArrays

function diam(vertices)
    submatrix = M_dist[vertices, vertices]
    if (submatrix .|> iszero |> sum) > length(vertices)
        return Inf
    else
        return maximum(M_dist[vertices, vertices])
    end
end

Random.seed!(0)
M2n = rand(Sphere(2), 10)
d, n = size(M2n)
ε = 2.0
# dims = n-1
dims = 8

# Initialization
kdtree = KDTree(M2n)
vrtx, dist = knn(kdtree, M2n, dims+1, true)
M_dist = spzeros(Float64, n,n)
for i in 1:n
  for (idxj, j) in enumerate(vrtx[i])
    M_dist[i,j] = dist[i][idxj]
  end
end

df = DataFrame(ID = UInt64[], Dimension = Int64[],
               Appearance = Float64[], Simplex = Set{Int64}[])
Δ = Dict()
for k ∈ 0:dims
    push!(Δ, k::Int64 => Set{Set{Int64}}())
    if k == 0
        for j ∈ 1:n
            splx = Set([j])
            push!(Δ[k], splx)
            push!(df, [hash(splx), 0, 0, splx])
        end
        continue
    end
    for face in Δ[k-1]
        for v ∈ 1:n
            if v ∉ face
                candy = [collect(face) ; v]
            else
                continue
            end
            δ = diam(candy)
            if 0 < δ ≤ ε
                splx = Set(candy)
                push!(Δ[k], splx)
                push!(df, [hash(splx), k, δ, splx])
            end
        end
    end
    # println(k)
    # # k = 6
    # δ = last.(dist)
    # id_VR = findall(δ .≤ ε)
    # if !(id_VR |> isempty)
    #     for id_splx ∈ id_VR
    #         splx = Set(vrtx[id_splx])
    #         push!(Δ[k], Set(splx))
    #     end
    # end
    # map(pop!, vrtx)
    # map(pop!, dist)
end
Δ = sort(Δ)
sort!(df, :Appearance, rev = true)
unique!(df, :ID)
show(stdout, "text/plain", df)
sort!(df, :Dimension, rev = true)

print(1)
# for k in dims:-1:0
# # k = 0
#     vert = view.(prevert, Ref((0:k+1).+1))
#     dist = getindex.(predist, k+1)
#     bit_VR = dist.≤2ε
#     preΔ[k+1] = vert[bit_VR] .|> sort |> unique
#     for β in preΔ[k]
#         # flag_Δ = false
#         # for α in preΔ[k+1]
#         #     # Skipping if β already in Δ
#         #     # If we have {[1], [2], [3], [1,2]},
#         #     # then we want {[3], [1,2]}
#         #     if β.V ⊆ α.V
#         #         flag_Δ = true
#         #         break
#         #     end
#         # end
#         # if flag_Δ continue end
#         push!(Δ[k], β)
#     end

#     # Setting next dimension or escape
#     if preΔ[k+1] |> isempty
#         break
#     else
#         Δ[k+1] = Int64[]
#     end
# end
#    Δ = sort(Δ)
# preΔ = sort(preΔ)


x, y = eachrow(M2n);
p1 = scatter(x, y, text = 1:n, size = (600, 600), color = :white)
for (idx, δ) ∈ enumerate(reduce(∪, dist) |> sort)
    p2 = plot(p1)
    for t = 1:10
        p2 = circle!(p2, (x[t], y[t]), δ/2, label = :none)
    end
    png(p2, "D:/TDAexmaple/$idx-δ=$δ.png")
end
