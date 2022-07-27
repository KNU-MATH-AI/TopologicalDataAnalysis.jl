include("../src/homology/zomorodian.jl")

filtered_complex = DataFrame([
0 Set([1])
0 Set([2])
1 Set([3])
1 Set([4])
1 Set([1, 2])
1 Set([2, 3])
2 Set([3, 4])
2 Set([1, 4])
3 Set([1, 3])
4 Set([1, 2, 3])
5 Set([1, 3, 4])
], ["degree", "simplex"])

p_intervals = zomorodian(filtered_complex)