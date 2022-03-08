struct Interval{T<:Real}
    a::T
    b::T
    # Interval(a,b) = a > b ? error("out of order") : new(a,b)
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

mutable struct MapperInformations
    data_number
    original_coordinate
    original_subset
    clustered_subset
    filtered_coordinate
    interval_subset
end

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

    X_subset = interval_test(Info, Intervals)

    Info

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

    center_ = coor_mean(X_subset_clustering, Info)

    center_number = 0
    for i ∈ 1:length(X_subset)
        for j ∈ 1:length(X_subset_clustering[i])
            center_number +=1
        end
    end

    center_matrix = zeros(2, center_number)
    k_temp = 1
    for i ∈ 1:length(X_subset)
        for j ∈ 1:length(X_subset_clustering[i])
            center_matrix[:,k_temp] = center_[i][j]
            k_temp += 1
        end
    end

    zcolor_temp = filter(cat(Data, center_matrix, dims=2))

    zcolor = [[] for i ∈ 1:length(X_subset)]
    k_temp = 1
    for i ∈ 1:length(X_subset)
        for j ∈ 1:length(X_subset_clustering[i])
            push!(zcolor[i], zcolor_temp[size(Data)[2]+k_temp])
            k_temp +=1
        end
    end
    
    return link, X_subset_clustering, Info, center_, zcolor
end

function mapper_plot(link, X_subset_clustering, center_, zcolor, xlim, ylim)

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
            scatter!([center_[i][j][1]], [center_[i][j][2]], markersize=25*(452.106/181)*sqrt(length(X_subset_clustering[i][j])/500), zcolor=zcolor[i][j], label="", color=:jet, clim=(0,Inf))
            length_temp = length(X_subset_clustering[i][j])
            annotate!([center_[i][j][1]], [center_[i][j][2]], "$length_temp")
        end
    end

    return p
end