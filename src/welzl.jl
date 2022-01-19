Point = Vector{Float64}

struct Ball
    c::Point
    r::Float64
end

x = Ball([1,3,3,2,1],2)
x.c
x.r


norm(Point) = sum(abs2, Point) |> sqrt

x = [1,2,3]
t = [4,5,6]

x*t

using BoundingSphere

P = [randn(2) for _ in 1:10000]
@time x = boundingsphere(P)


using DataFrames

