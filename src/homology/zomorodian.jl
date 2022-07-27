using NearestNeighbors
using ProgressBars

"""
dim(σ::Array{Int8,1})

Return the length of given set minus 1. We may assume that σ is a simplex with integer vertices.
"""
dim(σ::Array{Int8,1}) = length(σ) - 1


"""
deg(T::DataFrame, σ::Array{Int8,1})

Return the degree of simplex σ in the given filtered_complex T. T must have the hased values of simplices.
"""
function deg(T::DataFrame, σ::Array{Int8,1})
    return T.degree[findfirst(T.ID .== hash(σ))]
end

"""
∂(σ::Array{Int8,1})

The boundary operator of a simplex. It returns the array of sets which is {σ ∖ s : s ∈ σ}.
"""
function ∂(σ::Array{Int8,1})
    return [Array{Int8,1}(σ[Not(t)]) for t ∈ 1:length(σ)]
end

"""
maxindex(T::DataFrame, chain)

Find the maximual index of chain in filtered complex T. 'chain' is an array of simplices.
"""
function maxindex(T::DataFrame, chain)
    return (T.simplex .∈ Ref(chain)) |> findall |> maximum
end

"""
REMOVEPIVOTROWS!(T::DataFrame, σ::Array{Int8,1})

Return a differentiated chain without unmarked simplex.
"""
function REMOVEPIVOTROWS!(T::DataFrame, σ::Array{Int8,1})
    k = dim(σ); d = ∂(σ)
    d = d[d .∈ Ref(T[T.marked,:simplex])] # Remove unmarked terms in $d$
    while !(d |> isempty)
        i = maxindex(T, d)
        if T[i,:slot] |> isempty break end
        d = symdiff(d, T[i,:slot])
        # print("d in empty")
    end
    return d
end

"""
zomorodian(filtered_complex::DataFrame)

Return a list of P-intervals and print how it coputed.

     - filtered_complex: A dataframe with two columns, :degree and :simplex.

"""
function zomorodian(filtered_complex::DataFrame)
    T = copy(filtered_complex)
    m = nrow(T)
    T[!, :"marked"] .= false
    T[!, :"slot"] .= [[]]
    T[!, :"J"] .= 0
    T[!,:"ID"] = hash.(filtered_complex.simplex)

    L_ = Dict([k => [] for k = 0:maximum(dim.(T.simplex))])
    for j ∈ ProgressBar(1:m)
        σʲ = T[j,:simplex]
        d = REMOVEPIVOTROWS!(T, σʲ)
        if d |> isempty
            T[j,:marked] = true
        else
            i = maxindex(T, d)
            σⁱ = T[i,:simplex]
            k = dim(σⁱ)
            T[i,[:J,:slot]] = (j-1), d
            push!(L_[k], (deg(T, σⁱ), deg(T, σʲ)))
        end
    end
    for j ∈ 1:m
        σʲ = T[j,:simplex]
        if (T[j,:marked]) && (T[j,:slot] |> isempty) && (T[j,:J] |> iszero)
            k = dim(σʲ)
            push!(L_[k], (deg(T, σʲ), Inf))
        end
    end

    # println(T)
    return L_
end
