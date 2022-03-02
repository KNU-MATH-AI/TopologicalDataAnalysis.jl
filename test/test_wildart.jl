using Dates
using LinearAlgebra
using SparseArrays
include("../src/wildart.jl")

# using Logging
# debug_logger = ConsoleLogger(stderr, Logging.Debug)
# global_logger(debug_logger)


# X = [1 2 2 7;
#      1 2 6 1;
#      4 7 1 10]
# Y = copy(X)
# snf(Y)[3]
# using SmithNormalForm
# smith(Y)

include("../src/boundary_map_Alex.jl")
bigY = bound_matrix_4(9,4,4) 
result = snf(bigY)

# 170141142824549539861647102733501356222 < typemax(Int128)
# @time for i in 1:10000
#     rand(BigInt) ÷ rand(BigInt)
# end

# M = rand(0:9,3,3)

# convert.(Real, M)

# x = -9223372036854775808; typeof(x)
# y = Int128(x); typeof(y)
# div(x,-1)
# div(y,-1)

# div(92233720368547758080,-1)

# typemax(BigInt)
# -9223372036854775808
# typemax(-9223372036854775808)

# convert.(Int128,rand(0:9, 4,4))