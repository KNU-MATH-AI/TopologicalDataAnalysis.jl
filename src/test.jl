include("cech.jl")

default(size = (400,400), legend = :none)
n = 100
Random.seed!(0)
P = [randn(2) for _ in 1:n]
pc = pointcloud(P)

@time Complices = AllComplices(P)

Complices[1][2]
Complices[2][2]
Complices[3][2]
Complices[4][2]
Complices[5][2]
Complices[end-20][2]
Complices[end-2][2]
Complices[end-1][2]

for C in Complices
    println(length(C[2]))
end

asdf = @animate for (k, complex) in enumerate(reverse(Complices))
    pc = pointcloud(P, title = "$k, ε = $(trunc(complex[1],digits = 2))", text = 1:10)
    for (c,r) = complex[2]
        circle!(pc,c,r)
    end
end; gif(asdf, fps = 2)

kdtree = KDTree(Point2Matrix(P))

Complices[end][2]

inrange(kdtree, Complices[end][2][1][1], Complices[end][2][1][2]+eps())

Complices[end][2][1] in Complices[end-1][2]