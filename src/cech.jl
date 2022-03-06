using BoundingSphere
using NearestNeighbors
using Plots

struct AbstractSimplicialComplex
    Simplex::Set{Int64}
end
const ASC = AbstractSimplicialComplex

include("tDATAsets.jl")
M2n = rand(Sphere(2), 10)
d, n = size(M2n)
x,y = eachrow(M2n); scatter(x,y, text = 1:n)
ε = 0.3



kdtree = KDTree(M2n)
knn0 = knn(kdtree, M2n, 2)
bit_Δ0 = first.(knn0[2]) .≥ ε
Δ0 = [[k] for k in 1:n][bit_Δ0]
preΔ1 = knn0[1][.!bit_Δ0] .|> sort |> unique
v2 = reduce(union, preΔ1)

Δ1 = []
preΔ2 = []
for (idx_taboo, v0v1) in enumerate(preΔ1)
    flag_Δ1 = false
    for Δ2 in preΔ2
        if v0v1 ⊆ Δ2
            flag_Δ1 = true
            continue
        end
    end
    if flag_Δ1 continue end
    for v in v2[idx_taboo:end]
        if v in v0v1
            continue
        end
        v0v1v2 = [v0v1; v]
        println("candidate: ", v0v1v2)
        c, r = boundingsphere([M2n[:,idx] for idx in v0v1v2])
        if r ≤ ε && v0v1v2 == inrange(kdtree, c, r + eps(), true)
            flag_Δ1 = true
            println(inrange(kdtree, c, r + eps(), true))
            push!(preΔ2, v0v1v2)
            break
        end
    end
    push!(Δ1, v0v1)
end

Δk = copy(Δ1)
preΔkpp = []

cech(x) = x .^ 2