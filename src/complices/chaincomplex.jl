# function ⋃(Sets...)
#     return reduce(∪, Sets)
# end
const ⋃ = union

function ∂(Simplex)
    # if !(Simplex |> issorted) @warn "Simplex $Simplex is not sorted" end
    faces = [setdiff(Simplex, v) for v in reverse(Simplex)]
    return faces
end

function faces(Simplex::Set)
    return Set([setdiff(Simplex, v) for v in Simplex])
end
faces(Set([3,5,3,2]))
∂([1,2,3])
∂([1,3,2])