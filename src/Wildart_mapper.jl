# Example of Wildart
using Plots
using Pkg

Pkg.add("TDA")

using TDA

# generate and plot some dataset
X = hcat(TDA.circlepoints(500, 0.5, noise=0.05)...)'
plot(X[1,:], X[2,:], seriestype=:scatter)

# define Mapper filter function for dataset: f(x) = ||x.x - p.x||
fltfn = (data)->vec(mapslices(p->p[1]-minimum(data[1,:]), data, dims=1))
# plot data colored by filter function values
plot(X[1,:], X[2,:], label="", zcolor=fltfn(X), seriestype=:scatter, ms=2)

# call Mapper algorithm with the particular filter function.
mpr = TDA.mapper(X, filter=fltfn, seed=0, intervals=5, overlap=0.2)

# plot topological layout - mapper graph (by default circular layout is used)
plot(mpr, c=:viridis)

# use `constant_layout` for positioning Mapper graph vertices
# at centers of cover patches
plot(mpr, c=:viridis, complex_layout=TDA.constant_layout)