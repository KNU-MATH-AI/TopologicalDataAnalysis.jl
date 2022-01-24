include("../src/smith.jl")

X = [1 2 2 7;
     1 2 6 1;
     4 7 1 0]; Y = copy(X)
gauss!(X)
rswap!(X,1,2)


@time smith(Y)
@time smith(Y)