@time using NearestNeighbors
@time using BoundingSphere
@time using Plots
@time using Random
@time using DynamicsPlots

struct Ball
    c::Vector{Float64}
    r::Float64
end

function Point2Matrix(P)
    M = zeros(length(P[1]), length(P))
    for (k, p) in enumerate(P)
        M[:,k] = p
    end
    return M
end

L2(M) = sqrt.(sum(abs2, M, dims = 1))
function ∂(c,r,P)
    boundary = vec(((P .- Ref(c)) |> Point2Matrix |> L2) .≈ r)
    closure = (boundary) .|| vec(((P .- Ref(c)) |> Point2Matrix |> L2) .≤ r)
    return boundary, closure
end
    
function simplex!(pc, Q)
    m = length(Q)
    for i in 1:m
        for j in 1:i
            plot!(pc,
            [Q[i][1],Q[j][1]],
            [Q[i][2],Q[j][2]],
            color = :black,alpha = 0.5)
        end
    end
    return pc
end

function MinSphere(P::AbstractMatrix)
    d,n = size(P)
    return boundingsphere([P[:,i] for i in 1:n]) # Welzl algorithm,
    # https://github.com/JuliaFEM/BoundingSphere.jl/blob/master/src/geometry.jl
    # Note that the calculation is unefficient in BoundingSphere.jl/src/geometry.jl
    # We may replace this kdtree.
end
function MinSphere(P::AbstractVector)
    return boundingsphere(P)
end

function isMaxCell(Q, P, ε)
    if MinSphere(Q)[2] < ε
        return true
    end
    print("!")
    Qᶜ = P[P .∉ Ref(Q)]
    for pᵢ in Qᶜ
        if MinSphere(vcat(Q,[pᵢ]))[2] ≥ ε
            # pop!(Q)
            return false
        end
    end
    return true
end

function AllComplices(P, verbose = 0)
    c, ε = MinSphere(P)
    C = [(c,ε)]
    Complices = []
    i = 1

    while ε > 0
        push!(Complices, ε => copy(C))
        bndr, clsr = ∂(c,ε,P)
        Q = P[bndr]
        R = P[clsr]
        popat!(C,i)

        if verbose == 1 println(ε) end
        for q in Q
            check = isMaxCell(setdiff(R, [q]), R, ε)
            if check
                push!(C, MinSphere(setdiff(R, [q])))
            end
        end
        unique!(C)

        i = argmax([c[2] for c in C])
        c = C[i][1]
        ε = C[i][2]
    end
    return Complices
end
