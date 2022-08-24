using Plots

function persistencediagram(p_intervals, filtered_complex, Infval)

    δ_table = unique(filtered_complex[:,[1,2]], :degree)
    δ_dict = Dict(δ_table.degree .=> δ_table.appearance)
    
    lim = maximum(δ_table.appearance)

    fig = plot([0,lim], [0, lim],
        xlabel = "birth",
        ylabel = "death",
        label =:none,
        color = "black"
        )
    
    for k in 0:Infval
        if !isempty(get.(Ref(δ_dict), first.(p_intervals[k]), Infval))
            fig = scatter!(fig,
            get.(Ref(δ_dict), first.(p_intervals[k]), Infval),
            get.(Ref(δ_dict), last.(p_intervals[k]), Infval),
            markeralpha = 1,
            markershape = :x,
            label = "H$k",
            legend=:outertopright
            )
        end
    end
    return fig
end


