# include("../src/smith.jl")

X = [1 2 2 7;
     1 2 6 1;
     4 7 1 10]
Y = copy(X)

# gauss!(X)
# rswap!(X,1,2)

include("../src/smith3.jl")

# Y = rand(0:3, 100, 100)
@time smith2(Y)
# @time S_N_F(Y)

using SmithNormalForm, LinearAlgebra
include("../src/boundary_map_Alex.jl")
Y = bound_matrix_5(8,3,5)
@time result = snf(Y, inverse = false);
# diag(result[3])