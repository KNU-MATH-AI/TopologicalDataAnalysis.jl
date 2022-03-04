using StatsPlots 
using Clustering 
using Distances 
using Distributions

a = rand(Uniform(0,.8), 2, 20)
b = rand(Uniform(1.2,2), 2, 20)
c = a + cat(1.2ones(1,20),zeros(1,20), dims=1)
X = cat(a,b,c, dims=2)

scatt = scatter(X[1,:], X[2,:], label=false)

D_a = pairwise(Euclidean(), X, X)
SL = hclust(D_a, linkage=:single)
plot(SL)

propertynames(Hclust)

propertynames(SL)
SL.merges
SL.heights
SL.order
SL.linkage
SL.height
SL.labels
SL.merge
SL.method

plot(SL.heights)
SL.heights[56]
SL.heights[57]
SL.heights[58]
SL.heights[59]

kmeans_ = kmeans(reshape(SL.heights, (1,59)), 2)
propertynames(kmeans_)
@show kmeans_

cutree(SL, k=2)