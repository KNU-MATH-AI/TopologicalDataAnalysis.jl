using Distributions
using Plots
using NearestNeighbors
using Ripserer

include("../src/utilities/utils.jl")
include("../src/complices/VR.jl")
include("../src/homology/zomorodian.jl")
include("../src/utilities/tDATAsets.jl")

N = 20
dim_ = 3
# point_cloud = noisy_circle_data(N, 0)
point_cloud = rand(Sphere(dim_), N)
# scatter(point_cloud[1,:], point_cloud[2,:])

VR_KNU = vietoris_rips(point_cloud)
p_intervals = zomorodian(VR_KNU)

degree=3
δ_table = unique(VR_KNU[:,[1,2]], :degree)
δ_dict = Dict(δ_table.degree .=> δ_table.appearance)

point_cloud_ripserer = [(point_cloud[1,i], point_cloud[2,i], point_cloud[3,i]) for i ∈ 1:N]
VR_ripserer = ripserer(point_cloud_ripserer, dim_max=3)

p1 = plot(VR_ripserer)
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

# plot(VR_ripserer)
# barcode(VR_ripserer)

p_intervals
show(stdout, "text/plain", VR_ripserer[1])
show(stdout, "text/plain", VR_ripserer[1])