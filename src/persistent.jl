include("cech.jl")

default(size = (400,400), legend = :none)
n = 10
Random.seed!(0)
P = [randn(2) for _ in 1:n]
pc = pointcloud(P)

Complices = AllComplices(P)

kdtree = KDTree(Point2Matrix(P))
inrange(kdtree, Complices[1][2][1][1], Complices[1][2][1][2])

function vertices(kdtree, )

epsilon = Complices[2][1]
Simplilces = Complices[2][2]

for Simplex ∈ Simplilces
    c = Simplex[1]
    δ = Simplex[2]
end

(x = 1 => y=2)

Complices[1][2][1][1]