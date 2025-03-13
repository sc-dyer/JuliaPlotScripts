using CairoMakie, CSV, DataFrames
colourchoice = 2
include("PlotDefaults.jl")

function filterSampleName(df::DataFrame,sampleName::String)
    return filter([:Sample]=>sam -> occursin(sampleName,sam),df)
end


function plot_classification(filein,savedir,group1,group2)
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 2.0)
    limits!(ax,0.0,2.0,0.0,1.0)

    lines!(ax,[0.5,0.5],[0.0,1.0],color = :black)
    lines!(ax,[1.5,1.5],[0.0,1.0],color = :black)
    lines!(ax,[0.0,2.0],[0.5,0.5],color = :black)

    ampdata = DataFrame(CSV.File(filein))
    for i in 1:lastindex(group1)
        sample = group1[i]
        ampfiltered = filterSampleName(ampdata,sample)
        if i <=7
            scatter!(ax,ampfiltered[!,:C_classify],ampfiltered[!,:A_classify],color = Cycled(i),label = sample, marker =:circle,markersize = 10,strokewidth=0.5)
    
        else
            scatter!(ax,ampfiltered[!,:C_classify],ampfiltered[!,:A_classify],color = Cycled(i),label = sample, marker =:utriangle,markersize = 10,strokewidth=0.5)

        end
    end
    for i in 1:lastindex(group2)
        sample = group2[i]
        ampfiltered = filterSampleName(ampdata,sample)
        scatter!(ax,ampfiltered[!,:C_classify],ampfiltered[!,:A_classify],color = Cycled(i),label = sample, marker =:rect,markersize = 10,strokewidth=0.5)
    end

    ax.xlabel = "C"
    ax.ylabel = "A"
    fig[1,2] = Legend(fig,ax)
    save(joinpath(savedir,"AmpProbe.svg"),fig)
end

mystuff = ["23SD03C-1","23SD03C-3","23SD03E","22SD55B3","22SD55D","22SD55E","23SD20B","23SD20C","23SD20F","23SD20G"]
trondstuff = ["2M02061A","2M02061B","2M0506-16_2","2M0506-16_3B"]

plot_classification("../../Probe/AmphiboleCationSums.csv","ProbePlots/",mystuff,trondstuff)