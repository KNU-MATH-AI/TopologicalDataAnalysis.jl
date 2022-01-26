using Ripserer, ComputationalHomology, Plots

function PD(data, Backend)

    if Backend=="Vietoris-Rips"
        data = Tuple.(eachrow(data));
        plot(ripserer(data))
    end

    if Backend=="Cech"
        data = Matrix(data);
        ɛ = 0.01;
        cplx, w = čech(data, ɛ);
        pd = cplx |> filtration |> diagram;
        plot(pd, ms=4)
    end
    
end