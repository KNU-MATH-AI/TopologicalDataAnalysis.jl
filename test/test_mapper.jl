using StatsPlots 
using Clustering 
using Distances 
using Distributions
using DataFrames

import Base: ∈
import Base: ∩

include("../src/utilities/utils.jl")

# a = rand(Uniform(0,.8), 2, 20)
# b = rand(Uniform(1.2,2), 2, 20)
# c = a + cat(1.2ones(1,20),zeros(1,20), dims=1)
# X = cat(a,b,c, dims=2)

# scatt = scatter(X[1,:], X[2,:], label=false)

# D_a = pairwise(Euclidean(), X, X)
# SL = hclust(D_a, linkage=:single)
# CL = hclust(D_a, linkage=:complete)
# AV = hclust(D_a, linkage=:average)
# plot(SL)
# plot(CL)
# plot(AV)

# propertynames(Hclust)

# propertynames(SL)
# SL.merges
# SL.heights
# SL.order
# SL.linkage
# SL.height
# SL.labels
# SL.merge
# SL.method

# plot(SL.heights)
# SL.heights[56]
# SL.heights[57]
# SL.heights[58]
# SL.heights[59]

# kmeans_ = kmeans(reshape(SL.heights, (1,59)), 2)
# kmeans_ = kmeans(reshape(CL.heights, (1,59)), 2)
# propertynames(kmeans_)
# @show kmeans_

# cutree(SL, k=2)

# D_SL = pairwise(Euclidean(), SL.heights)
# dbscan(reshape(SL.heights, (1,59)), 0.8maximum(SL.heights))
# dbscan(reshape(SL.heights, (1,59)), 0.05maximum(SL.heights))
# dbscan(reshape(SL.heights, (1,59)), 0.2maximum(SL.heights))

# dbscan_ = dbscan(reshape(CL.heights, (1,59)), 0.1maximum(CL.heights))
# dbscan_[1]
# dbscan_[2]
# dbscan(reshape(CL.heights, (1,59)), 0.08maximum(CL.heights))
# dbscan(reshape(CL.heights, (1,59)), 0.2maximum(CL.heights))

# cutree(CL, h=CL.heights[58])
# x = Interval
# x = Interval(1,2)

# struct OrderedPair
#     x::Real
#     y::Real
#     OrderedPair(x,y) = x > y ? error("out of order") : new(x,y)
# end
# OrderedPair(3,1)

# I = Interval(0.2, 4.0)
# -1 ∈ I
# 1 ∈ I
# J = Interval(-1., 3.6)

# # function intersection(I::Interval, J::Interval)
# #     min = maximum([I.a, J.a])
# #     max = minimum([I.b, J.b])

# #     return Interval(min, max)
# # end

# ∩(I,J)

begin
using StatsPlots 
using Clustering 
using Distances 
using Distributions
using DataFrames

import Base: ∈
import Base: ∩

include("../src/utilities/utils.jl")
Data = noisy_circle_data(500, 0.05)

function proj_1d(X)
    min = minimum(X[1, :])
    output = X[1, :] .- min
    
    return output
end

struct Interval{T<:Real}
    a::T
    b::T
    # Interval(a,b) = a > b ? error("out of order") : new(a,b)
end

mutable struct MapperInformations
    data_number
    original_coordinate
    original_subset
    clustered_subset
    filtered_coordinate
    interval_subset
end

function ∈(x::T, interval::Interval) where T<:Real
    if interval.a ≤ x ≤ interval.b
        return true
    else
        return false
    end
end

function ∩(I::Interval, J::Interval)
    min = maximum([I.a, J.a])
    max = minimum([I.b, J.b])

    return Interval(min, max)
end

function filtered_coordinate(x, p)
    x.filtered_coordinate = p
end

#Y_ = MapperInformations(:None, :None, :None, :None, :None)
end

function mapper(X, filter, interval, overlab)
    
    original_X = []
    for i in 1:size(X,2)
        push!(original_X, X[:,i])
    end
    # findfirst(==(Info[500].original_coordinate), original_X)

    Info = [MapperInformations(i, X[:, i], Int64[], Tuple{Int64, Int64}[], :None, :None) for i ∈ 1:length(X[1,:])]

    filtered_X = filter(X)
    filtered_coordinate.(Info, filtered_X)
    Info
    
    min_X = minimum(filtered_X)
    max_X = maximum(filtered_X)
    filter_range = Interval(min_X, max_X)

    length_ = (max_X - min_X)/interval
    ϵ = overlab*length_

    Intervals = [∩(Interval(min_X +(i-1)*length_ - ϵ, min_X +i*length_ + ϵ), filter_range) for i ∈ 1:interval]

    function interval_test(Info, Intervals)
        M = length(Info)
        N = length(Intervals)
        
        X_subset = [[] for i ∈ 1:N]
        for i ∈ 1:M
            for j ∈ 1:N
                if Info[i].filtered_coordinate ∈ Intervals[j]
                    Info[i].interval_subset = Intervals[j]
                    Info[i].original_subset = push!(Info[i].original_subset, j)
                    push!(X_subset[j], Info[i].original_coordinate)
                end
            end
        end

        return X_subset
    end

    X_subset = interval_test(Info, Intervals)

    Info
    
    function domain_clustering(Info, X_subset)
        
        X_subset_clustering = []

        for i ∈ 1:length(X_subset)
            # i=3
            temp = zeros(length(X_subset[i][1]), length(X_subset[i]))
            for j ∈ 1:length(X_subset[i])
                temp[:, j] = X_subset[i][j]
            end
            scatter(temp[1,:], temp[2,:])               

            D = pairwise(Euclidean(), temp, temp)
            
            SL = hclust(D, linkage=:single)
            plot(SL)
            
            dbsacn_temp = dbscan(reshape(SL.heights, (1,length(SL.heights))), 0.2maximum(SL.heights))
            cutree_temp = cutree(SL, h=SL.heights[dbsacn_temp[1].core_indices[end]])
            
            X_subset_clustering_temp = [[] for j ∈ 1:maximum(cutree_temp)]
            
            for n ∈ 1:length(X_subset[i])
                clustering_index = cutree_temp[n]
                data_index = findfirst(==(X_subset[i][n]), original_X)
                push!(X_subset_clustering_temp[clustering_index], data_index)                
                Info[data_index].clustered_subset = push!(Info[data_index].clustered_subset, (i, clustering_index))
            end

            push!(X_subset_clustering, X_subset_clustering_temp)
            
        end
        return X_subset_clustering
    end

    X_subset_clustering = domain_clustering(Info, X_subset)

    # function myshowall(io, x, limit = false) 
    #     println(io, summary(x), ":")
    #     Base.print_matrix(IOContext(io, :limit => limit), x)
    # end
    # myshowall(Base.stdout, Info, false)

    # X_subset_clustering[1]
    # X_subset_clustering[2]
    # X_subset_clustering[3]
    # X_subset_clustering[4]
    # X_subset_clustering[5]
    
    link = []
    
    for i ∈ 1:length(X_subset)-1
        for j ∈ 1:length(X_subset_clustering[i])
            for k ∈ 1:length(X_subset_clustering[i+1])
                if length(∩(X_subset_clustering[i][j], X_subset_clustering[i+1][k])) > 0
                    push!(link, [[i,j], [i+1,k]])
                end
            end
        end
    end
    
    return link, X_subset_clustering, Info
end

function coor_mean(X_subset_clustering, Info)
    center_ = [ [] for i ∈ 1:length(X_subset_clustering)]
    for i ∈ 1:length(X_subset_clustering)
        for j ∈ 1:length(X_subset_clustering[i])
            temp = [Info[v].original_coordinate for v ∈ X_subset_clustering[i][j]]
            center_vector = sum(temp)./length(temp)
            push!(center_[i], center_vector)
        end
    end

    return center_
end

function mapper_plot(link, X_subset_clustering, Info, xlim, ylim)
    center_ = coor_mean(X_subset_clustering, Info)        

    linked_center = zeros(length(link), 4)
    for i ∈ 1:length(link)
        linked_center[i, 1:2] = center_[link[i][1][1]][link[i][1][2]]
        linked_center[i, 3:4] = center_[link[i][2][1]][link[i][2][2]]
    end
    linked_center

    p = plot([linked_center[1,1], linked_center[1,3]], [linked_center[1,2], linked_center[1,4]], color=:black, xlim=xlim, ylim=ylim, label="")
    for i ∈ 2:length(link)
        plot!([linked_center[i,1], linked_center[i,3]], [linked_center[i,2], linked_center[i,4]], color=:black, label="")
    end

    for i ∈ 1:length(X_subset_clustering)
        for j ∈ 1:length(X_subset_clustering[i])
            scatter!([center_[i][j][1]], [center_[i][j][2]], markersize=25*(452.106/181)*sqrt(length(X_subset_clustering[i][j])/500), label="", )
        end
    end

    return p
end


d = scatter(Data[1,:], Data[2,:])
link, X_subset_clustering, Info = mapper(Data, proj_1d, 5, 0.2)
p = mapper_plot(link, X_subset_clustering, Info)

plot(d, p, size=(1000,400))

include("../src/utilities/tDATAsets.jl")

Data = rand(Sphere(2), 500)
Data = cat(rand(Infty(), 500), rand(Sphere(2), 500).+1, dims=2)

d = scatter(Data[1,:], Data[2,:], zcolor=proj_1d(Data), color=:jet)

plot(Data[1,:], Data[2,:], zcolor=proj_1d(X), seriestype=:scatter)
link, X_subset_clustering, Info = mapper(Data, proj_1d, 5, 0.3)
p = mapper_plot(link, X_subset_clustering, Info, (-1.5, 1.5), (-1.5,1.5))

plot(d, p, size=(1000,400))


d = scatter(Data[1,:], Data[2,:])
link, X_subset_clustering, Info = mapper(Data, proj_1d, 10, 0.2)
p = mapper_plot(link, X_subset_clustering, Info)

plot(d, p, size=(1000,400))