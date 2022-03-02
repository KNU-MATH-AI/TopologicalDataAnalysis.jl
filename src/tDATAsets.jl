"""
This functions are julia implementations of the functions in scikit-tda/tadasets.

Reference: https://github.com/scikit-tda/tadasets
"""

abstract type TopologicalDataType end
struct Infty <: TopologicalDataType end
struct SwissRoll <: TopologicalDataType end
struct Sphere <: TopologicalDataType
    d::Int64
end
struct Torus <: TopologicalDataType
    R::Float64
    r::Float64
end

function Sphere()
    return Sphere(2)
end
function Torus()
    return Torus(4.0, 1.0)
end

import Base.rand
"""
Infty

Sample random points on lemniscate.

Reference: https://math.stackexchange.com/a/7989/459895
"""
function rand(S::Infty, n::Int64; noise = 0.01)
    t = 2π * rand(Float64, n)
    x = cos.(t)
    y = sin.(t) .* cos.(t)
    return vcat(x',y') + noise * randn(2, n)
end

"""
Swiss Roll
Reference: http://people.cs.uchicago.edu/~dinoj/manifold/swissroll.html
"""
function rand(S::SwissRoll, n::Int64; noise = 0.01)
    t = 2π * rand(Float64, n)
    y = 2π * rand(Float64, n)
    x = t .* cos.(t)
    z = t .* sin.(t)
    return vcat(x',y',z') + noise * randn(3, n)
end


"""
d-dimensional Sphere

Reference: https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
"""
function rand(S::Sphere, n::Int64; noise = 0.01)
    direction = randn(S.d, n)
    boundary = direction ./ sqrt.(sum(abs2, direction, dims = 1))
    return boundary + noise * randn(S.d, n)
end

"""
Torus
Reference: https://en.wikipedia.org/wiki/Torus#Geometry
"""
function rand(S::Torus, n::Int64; noise = 0.01)
    u = 2π * rand(Float64, n)
    v = 2π * rand(Float64, n)
    x = (S.R .+ S.r * cos.(v)) .* cos.(u)
    y = (S.R .+ S.r * cos.(v)) .* sin.(u)
    z = S.r * sin.(v)
    return vcat(x',y',z') + noise * randn(3, n)
end