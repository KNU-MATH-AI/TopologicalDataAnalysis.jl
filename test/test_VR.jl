include("../src/utilities/tDATAsets.jl")
include("../src/complices/VR.jl")
include("../src/homology/zomorodian.jl")

using Random
using Plots
using Ripserer

Random.seed!(0)
point_colud = rand(Sphere(2), 15)

@time fc = VR(point_colud, max_epsilon = 2.0)
p_intervals = zomorodian(fc)

δ_table = unique(fc[:,[1,2]], :degree)
δ_dict = Dict(δ_table.degree .=> δ_table.appearance)

result_rips = ripserer([(x,y) for (x,y) in eachcol(point_colud)]; dim_max=3)
p1 = plot(result_rips)

for k in 0:2
    Infval = 2
    # plot([0,Infval], [0,Infval], color = :black, size = (400,400), legend = :bottomright, label = :none)
    scatter!(p1,
        get.(Ref(δ_dict), first.(p_intervals[k]), Infval),
        get.(Ref(δ_dict), last.(p_intervals[k]), Infval),
        markeralpha = 0.5,
        markershape = :x,
        label = "H$k"
    )
end
p1