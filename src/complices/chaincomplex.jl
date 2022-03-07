function ∂(Simplex)
    # if !(Simplex |> issorted) @warn "Simplex $Simplex is not sorted" end
    faces = [setdiff(Simplex, v) for v in reverse(Simplex)]
    return faces
end

∂([1,2,3])
∂([1,3,2])