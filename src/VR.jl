using NearestNeighbors
using Plots
using DynamicsPlots
using Random

struct AbstractSimplicialComplex
    Simplex::Set{Int64}
end
const ASC = AbstractSimplicialComplex

include("tDATAsets.jl")
Random.seed!(0)
M2n = rand(Sphere(2), 10)
d, n = size(M2n)
ε = 0.2
maxdim = 6

# Initialization
kdtree = KDTree(M2n)

Δ = []
preΔ = []
push!(preΔ, [[j] for j in 1:n])

k = 4
push!(Δ, [])
knnk = knn(kdtree, M2n, k+1, true)
push!(preΔ, knnk[1][last.(knnk[2]) .≤ 2ε] .|> sort |> unique)
for β in preΔ[k]
    flag_Δ = false
    for α in preΔ[k+1]
        if β ⊆ α
            flag_Δ = true
            break
        end
    end
    if flag_Δ continue end
    push!(Δ[k], β)
    println(Δ[k])
end
Δ
preΔ




x,y = eachrow(M2n); p1 = scatter(x,y, text = 1:n, size = (600,600), color = :white)