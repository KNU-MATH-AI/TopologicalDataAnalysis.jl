include("../src/homology/zomorodian.jl")

filtered_complex = DataFrame([
0 [[1]]
0 [[2]]
1 [[3]]
1 [[4]]
1 [[1, 2]]
1 [[2, 3]]
2 [[3, 4]]
2 [[1, 4]]
3 [[1, 3]]
4 [[1, 2, 3]]
5 [[1, 3, 4]]
], ["degree", "simplex"])

p_intervals = zomorodian(filtered_complex)