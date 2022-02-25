using Ripserer, ComputationalHomology, Plots

function PD(; 
    data, 
    complex, 
    dim_max=1,
    cutoff=0,
    modulus=2,
    field=Mod{modulus},
    verbose=false,
    alg=:cohomology,
    reps=alg == :cohomology ? false : 1:dim_max,
    implicit=alg != :homology)

    if complex=="Vietoris-Rips"
        result_rips = ripserer(Tuple.(eachrow(data)), 
        dim_max=dim_max, 
        cutoff=cutoff, 
        modulus=modulus, 
        field=field, 
        verbose=verbose,
        alg=alg, 
        reps=reps, 
        implicit=implicit)
        return plot(result_rips)
    end   

    if complex=="Cech"
        data = Matrix(data);
        ɛ = 0.01;
        cplx, w = čech(data, ɛ);
        pd = cplx |> filtration |> diagram;
        return plot(pd, ms=4)
    end
    
end