using Ripserer, PersistenceDiagrams, Plots
include("ComputationalHomology.jl")

function PD(; 
    data, 
    complex, 
    dim_max=1,
    cutoff=0,
    verbose=false,
    alg=:cohomology,
    reps=alg == :cohomology ? false : 1:dim_max,
    implicit=alg != :homology,
    weights = true,
    ɛ = 1,
    filtration_weights = true
    )

    if complex=="Vietoris-Rips"
        pd = ripserer(Tuple.(eachrow(data)), 
        dim_max=dim_max, 
        cutoff=cutoff, 
        verbose=verbose,
        alg=alg, 
        reps=reps, 
        implicit=implicit)
        return plot(pd)
    end       

    if complex=="Cech"
        cplx, w = ComputationalHomology.čech(Matrix(data), ɛ, filtration_weights);
        pd = ComputationalHomology.diagram(ComputationalHomology.filtration(cplx));

        temp = NTuple{2}[];
        pd_converted = Array{PersistenceDiagram, 1}();

        for dim = 0:length(pd)-1
            for int_num = 1:length(pd[dim])
                push!(temp, (pd[dim][int_num].b, pd[dim][int_num].d))
            end
            push!(pd_converted, PersistenceDiagram(PersistenceInterval.(temp)));
        end

        return plot(pd_converted)
    end
end


