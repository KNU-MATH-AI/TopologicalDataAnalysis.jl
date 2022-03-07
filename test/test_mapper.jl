using StatsPlots 
using Clustering 
using Distances 
using Distributions
using

a = rand(Uniform(0,.8), 2, 20)
b = rand(Uniform(1.2,2), 2, 20)
c = a + cat(1.2ones(1,20),zeros(1,20), dims=1)
X = cat(a,b,c, dims=2)

scatt = scatter(X[1,:], X[2,:], label=false)

D_a = pairwise(Euclidean(), X, X)
SL = hclust(D_a, linkage=:single)
CL = hclust(D_a, linkage=:complete)
AV = hclust(D_a, linkage=:average)
plot(SL)
plot(CL)
plot(AV)

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
kmeans_ = kmeans(reshape(CL.heights, (1,59)), 2)
propertynames(kmeans_)
@show kmeans_

cutree(SL, k=2)

D_SL = pairwise(Euclidean(), SL.heights)
dbscan(reshape(SL.heights, (1,59)), 0.8maximum(SL.heights))
dbscan(reshape(SL.heights, (1,59)), 0.05maximum(SL.heights))
dbscan(reshape(SL.heights, (1,59)), 0.2maximum(SL.heights))

dbscan_ = dbscan(reshape(CL.heights, (1,59)), 0.1maximum(CL.heights))
dbscan_[1]
dbscan_[2]
dbscan(reshape(CL.heights, (1,59)), 0.08maximum(CL.heights))
dbscan(reshape(CL.heights, (1,59)), 0.2maximum(CL.heights))

cutree(CL, h=CL.heights[58])


Data = noisy_circle_data(200, 0.05)

function proj_1d(X)
    min = minimum(X[1, :])
    output = X[1, :] .- min
    
    return output
end

proj_1d(Data)

struct Interval
    a::Real
    b::Real
end

function ∈(x::Real, interval::Interval)
    if interval.a ≤ x ≤ interval.b
        return true
    else
        return false
    end
end

I = Interval(0.2 ,4)
J = Interval(-1,3.6)

function ∩(I::Interval, J::Interval)
    min = maximum(I.a, J.a)
    max = minimum(I.b, J.b)

    return Interval(min, max)
end

function intersection(I::Interval, J::Interval)
    min = maximum([I.a, J.a])
    max = minimum([I.b, J.b])

    return Interval(min, max)
end

∩(I,J)
intersection(I,J)


function mapper(X, filter, interval, overlab)
    filtered_X = filter(X)

    min_X = minimum(filtered_X)
    max_X = maximum(filtered_X)
    filter_range = Interval(min_X, max_X)

    length = (max_X - min_X)/interval
    ϵ = overlab*length

    Intervals = [intersection(Interval(min_X +(i-1)*length - ϵ, min_X +i*length + ϵ), filter_range) for i ∈ 1:interval]

    # Intervals = [(min_X , min_X + length + ϵ)]
    # for i ∈ 2:intervals-1
    #     push!(Intervals, (min_X + (i-1)*length - ϵ , min_X + i*length + ϵ))
    # end
    # push!(Intervals, (min_X + (intervals-1)*length - ϵ , min_X + intervals*length))

    return Intervals

end

mapper(Data, proj_1d, 5, 0.2)