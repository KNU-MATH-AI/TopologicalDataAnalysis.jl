using Distributions
using Plots
using NearestNeighbors
using Ripserer

include("../src/utilities/utils.jl")
include("../src/complices/VR.jl")
include("../src/homology/zomorodian.jl")

N = 20
point_cloud = noisy_circle_data(N, 0)
scatter(point_cloud[1,:], point_cloud[2,:])

VR_KNU = vietoris_rips(point_cloud)
p_intervals = zomorodian(VR_KNU)

point_cloud_ripserer = [(point_cloud[1,i], point_cloud[2,i]) for i ∈ 1:N]
VR_ripserer = ripserer(point_cloud_ripserer, dim_max=3)

plot(VR_ripserer)
barcode(VR_ripserer)