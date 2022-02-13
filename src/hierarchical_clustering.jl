using StatsPlots 
using Clustering 
using Distances 
using Distributions
using Plots


a = rand(Uniform(0,1), 2, 20)
b = rand(Uniform(1,2), 2, 20)
X = cat(a,b, dims=2)
scatt = scatter(X[1,:], X[2,:], label=false)

D_a = pairwise(Euclidean(), X, X)

SL = hclust(D_a, linkage=:single)
plot_SL = plot(SL)

p = plot(scatt, plot_SL, size=(800,400))

#=============================================================================#
a = [ 1 5 0.5 2.5 3;
      5 2  4  1   6]

SL = hclust(pairwise(Euclidean(), a, a), linkage=:single)
CL = hclust(pairwise(Euclidean(), a, a), linkage=:complete)
GA = hclust(pairwise(Euclidean(), a, a), linkage=:average)


p_sl = plot(SL, title="SL")
p_cl = plot(CL, title="CL")
p_ga = plot(GA, title="GA")

l = @layout()
p = plot(p_sl, p_cl, p_ga, layout=(1,3), size=(728,300))
savefig(p, "hierarchical_clustering.png")


a = rand(Uniform(-1,1), 2, 25)
scatt = scatter(a[1,:], a[2,:], label=false)
savefig(scatt, "julia_hclust_scatter.png")

D_a = pairwise(Euclidean(), a, a)
SL = hclust(D_a, linkage=:single)
plot_SL = plot(SL)

p = plot(scatt, plot_SL, size=(800,400))
savefig(p, "julia_hclust.png")


a = rand(2, 10)
D_a = pairwise(Euclidean(), a, a)
SL = hclust(D_a, linkage=:single)
dendrogram = plot(SL)
savefig(dendrogram, "julia_dendrogram.png")

p = plot(scatt, plot_SL, size=(800,400))
savefig(p, "julia_hclust.png")


a = rand(Uniform(0,1), 2, 20)
b = rand(Uniform(1,2), 2, 20)
X = cat(a,b, dims=2)
scatt = scatter(X[1,:], X[2,:], label=false)

D_a = pairwise(Euclidean(), X, X)

SL = hclust(D_a, linkage=:single)
plot_SL = plot(SL)

p = plot(scatt, plot_SL, size=(800,400))
savefig(p, "julia_hclust.png")

a = rand(Uniform(0,.8), 2, 20)
b = rand(Uniform(1.2,2), 2, 20)
c = a + cat(1.2ones(1,20),zeros(1,20), dims=1)
X = cat(a,b, dims=2)

N = 200
r = rand(Uniform(0.8,1.2), N)
theta = rand(Uniform(0,2π), N)
X = zeros(2,N)
for i ∈ 1:N
      X[:,i] = [r[i]*cos(theta[i]), r[i]*sin(theta[i])]
end
scatt = scatter(X[1,:], X[2,:], label=false)

D_a = pairwise(Euclidean(), X, X)
SL = hclust(D_a, linkage=:single)
plot_SL = plot(SL)
p = plot(scatt, plot_SL, size=(800,400))
savefig(p, "circle_clust.png")