using NearestNeighbors
using Plots
using DynamicsPlots
using Random

struct AbstractSimplicialComplex
    Simplex::Set{Int64}
end
const ASC = AbstractSimplicialComplex

include("../utilities/tDATAsets.jl")
Random.seed!(0)
M2n = rand(Sphere(2), 10)
d, n = size(M2n)
ε = 0.2
maxdim = 4

# Initialization
kdtree = KDTree(M2n)
knnk = knn(kdtree, M2n, maxdim+1, true)
prevert = knnk[1] #pre-vertex
prediam = knnk[2] #pre-diameter
map(popfirst!,prediam)

Δ = []
preΔ = []
push!(preΔ, [[j] for j in 1:n])

for k in 1:maxdim
    vert = view.(prevert, Ref(1:(k+1)))
    diam = getindex.(prediam, k)
    push!(Δ, [])
    push!(preΔ, vert[diam .≤ 2ε] .|> sort |> unique)
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
end
Δ
preΔ


knnk = knn(kdtree, M2n, k+1, true)

knnk[1]

x,y = eachrow(M2n); p1 = scatter(x,y, text = 1:n, size = (600,600), color = :white)