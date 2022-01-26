using Ripserer, Plots, DataFrames

function PD(data, Backend)
    if Backend=="Vietoris-Rips"
        data = Tuple.(eachrow(data));
        plot(ripserer(data))
    end
end