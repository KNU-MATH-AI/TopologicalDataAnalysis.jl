include("../src/tDATAsets.jl")

using Plots
default(size=(400, 400))

InftySymbol = rand(Infty(), 100)
x,y = eachrow(InftySymbol)
scatter(x,y, size = (400, 200))

S1 = rand(Sphere(), 100)
x,y = eachrow(S1)
scatter(x,y)

Plots.plotly()
default(msw = 0, ma = 0.5, ms = 1, mc = :black)

S2 = rand(Sphere(3), 1000)
x,y,z = eachrow(S2)
scatter(x,y,z)

SR = rand(SwissRoll(), 1000)
x,y,z = eachrow(SR)
scatter(x,y,z)

T = rand(Torus(), 1000, noise = 0.0)
x,y,z = eachrow(T)
scatter(x,y,z, lims = (-6,6), legend = :none)