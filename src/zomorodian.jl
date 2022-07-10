using DataFrames

data = DataFrame([
0 "a"
0 "b"
1 "c"
1 "d"
1 "ab"
1 "bc"
2 "cd"
2 "ad"
3 "ac"
4 "abc"
5 "acd"
], ["epsilon", "simplex"])
T = copy(data)
T[!, :"marked"] .= false
T[!, :"slot"] .= [[]]
T[!, :"J"] .= 0

dimK = 2
m = nrow(T)

test_σ = "abc"
dim(σ) = length(σ)

function deg(σ)
    return T.epsilon[findfirst(T.simplex .== σ)]
end
deg(test_σ)

function ∂(σ)
    k = dim(σ)
    return [σ[(1:k)[Not(t)]] for t = 1:k]
end
∂(test_σ)

function maxindex(chain)
    return (T.simplex .∈ Ref(chain)) |> findall |> maximum
end
maxindex(∂(test_σ))


function REMOVEPIVOTROWS(σ)
    k = dim(σ); d = ∂(σ)
    d = d[d .∈ Ref(T[T.marked,:simplex])] # Remove unmarked terms in $d$
    while !(d |> isempty)
        i = maxindex(d)
        if T[i,:slot] |> isempty break end
        d = symdiff(d, T[i,:slot])
        # print("d in empty")
    end
    return d
end
REMOVEPIVOTROWS(test_σ)

# coumpute_intervals
# L_ = Dict(); for k = 0:dimK L_[k] = [] end

T = copy(data)
T[!, :"marked"] .= false
T[!, :"slot"] .= [[]]
T[!, :"J"] .= 0

L_ = [[] for k = 0:dimK]
for j0 = 0:(m-1)
    j = j0+1
    σʲ = T[j,:simplex]
    d = REMOVEPIVOTROWS(σʲ)
    if d |> isempty
        T[j,:marked] = true
    else
        i = maxindex(d); k = dim(σʲ)
        σⁱ = T[i,:simplex]
        T[i,[:J,:slot]] = j0,d
        L_[k] = L_[k] ∪ [(deg(σⁱ), deg(σʲ))]
    end
end
for j0 = 0:(m-1)
    j = j0+1
    σʲ = T[j,:simplex]
    if (T[j,:marked]) && (T[j,:slot] |> isempty) && (T[j,:J] |> iszero)
        k = dim(σʲ); L_[k] = L_[k] ∪ [(deg(σʲ), Inf)]
        print("j: $j")
    end
end
L_