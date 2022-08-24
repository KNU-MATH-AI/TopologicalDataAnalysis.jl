include("../src/PersistenceDiagram.jl")
include("../src/utilities/tDATAsets.jl")
include("../src/complices/VR.jl")
include("../src/homology/zomorodian.jl")

using Random
using Plots

Random.seed!(0)
point_colud = rand(Sphere(2), 20)

@time fc = VR(point_colud)
p_intervals = zomorodian(fc)

δ_table = unique(fc[:,[1,2]], :degree)
δ_dict = Dict(δ_table.degree .=> δ_table.appearance)

Infval = 2

p = persistencediagram(p_intervals, fc, Infval)