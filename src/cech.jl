using BoundingSphere
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
M2n = rand(Sphere(2), 20)
d, n = size(M2n)
ε = 0.3
maxdim = 6

# Initialization
kdtree = KDTree(M2n)
knn0 = knn(kdtree, M2n, 2)
# bit_Δ0 = first.(knn0[2]) .≥ ε
# Δ0 = [[j] for j in 1:n][bit_Δ0]
begin
cache_r = Float64[]    
Δ = []
preΔ = []
# push!(preΔ, knn0[1][.!bit_Δ0] |> sort .|> sort |> unique)
push!(preΔ, [[j] for j in 1:n])
for k in 1:maxdim
    if preΔ[k] |> isempty break end
    push!(Δ, [])
    push!(preΔ, [])
    # flatΔ = reduce(∪, preΔ[k]) |> sort
    flatΔ = 1:n
    for (idx_taboo, β) in enumerate(preΔ[k])
        flag_Δ = false
        for α in preΔ[k+1]
            if β ⊆ α
                flag_Δ = true
                break
            end
        end
        if flag_Δ continue end
        for v in flatΔ[idx_taboo:end]
            if v in β
                continue
            end
            # println("candy - β: ", β, " ∪  {", v, "}")
            candy = [β; v] |> sort # `candy` is a alias of 'candidate'.
            c, r = boundingsphere([M2n[:,idx] for idx in candy])
            checked = inrange(kdtree, c, r + eps(), true)
            if r ≤ ε && (length(candy) < length(checked) || candy == checked)
                flag_Δ = true
                push!(preΔ[k+1], candy)
                push!(cache_r, r)
                continue
            end
        end
        if flag_Δ continue end
        push!(Δ[k], β)
        println(Δ[k])
    end
end
end
# Δ0
Δ
preΔ[end]
cache_r |> unique |> sort
union(reduce(∪, Δ0), reduce(∪, reduce(∪, Δ))) |> sort

candy = [5,8,11,12]
c, r = boundingsphere([M2n[:,idx] for idx in candy])
inrange(kdtree, c, r + eps(), true)

color = [:yellow, :gray, :orange, :red, :black, :purple]
x,y = eachrow(M2n); p1 = scatter(x,y, text = 1:n, size = (600,600), color = :white)
for i in 1:length(Δ)
    for j in Δ[i]
        c, r = boundingsphere([M2n[:,k] for k in j])
        circle!(c,r, color = color[i], label = :none)
    end
end
circle!((0,0), ε, color = :blue, legend= :none)

c, r = boundingsphere([M2n[:,k] for k in [4,15,20]])
circle!(c, r, color = :green, legend= :none)
circle!(c, ε, color = :blue, legend= :none)
