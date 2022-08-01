include("../src/homology/zomorodian.jl")

filtered_complex = DataFrame([
0 [Int8[1]]
0 [Int8[2]]
1 [Int8[3]]
1 [Int8[4]]
1 [Int8[1, 2]]
1 [Int8[2, 3]]
2 [Int8[3, 4]]
2 [Int8[1, 4]]
3 [Int8[1, 3]]
4 [Int8[1, 2, 3]]
5 [Int8[1, 3, 4]]
], ["degree", "simplex"])

p_intervals = zomorodian(filtered_complex)