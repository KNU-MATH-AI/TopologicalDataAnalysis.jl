include("../src/PersistenceDiagram.jl")

using DataFrames

data = DataFrame(x = 0:0.1:1, y=0:0.1:1)

PD(data = data, complex = "Vietoris-Rips")

PD(data = data, complex = "Cech", ɛ = 0.4)