using CairoMakie
using DataFrames
using CSV
const MY_CATS = ["Host","Selvage","Pegmatite"]
const SLAGSTAD_CATS = ["Granodiorite mesosome","Stromatic leucosome","Slagstad pegmatite","Diorite mesosome", "Patch leucosome"]
const PLOT_VARS = ["Al2O3" "CaO"; "Na2O" "K2O"; "Fe2O3(T)" "MgO"]
colourchoice=2
include("PlotDefaults.jl")
set_theme!(myTheme)
# @recipe(TernaryText, x, y, z) do scene
#     Attributes(
#         text = "", 
#         marker = :circle,
#         markersize = 8,
#     )
# end

function variationdiagrams(filename,savedir)
    fig = Figure(;size = (800,1000));
    gdf = DataFrame(CSV.File(filename))
    for i in 1:3
        for j in 1:2
            ax = Axis(fig[i, j], aspect = 1.0)
            ax.ylabel = PLOT_VARS[i,j]
            for k in 1:lastindex(SLAGSTAD_CATS)
                cat = SLAGSTAD_CATS[k]
              
                gdfiltered = filter(:Category => c -> c == cat, gdf)
                scatter!(ax,gdfiltered[!,:SiO2],gdfiltered[!,PLOT_VARS[i,j]],color = Cycled(k), label = cat,markersize = 10,marker=:rect,strokewidth = 0.5)
            end
        
            for k in 1:lastindex(MY_CATS)
                cat = MY_CATS[k]
                
                gdfiltered = filter(:Category => c -> c == cat, gdf)
                scatter!(ax,gdfiltered[!,:SiO2],gdfiltered[!,PLOT_VARS[i,j]], color = Cycled(k),label = cat,markersize = 10,marker=:circle,strokewidth = 0.5)
                text!(ax,gdfiltered[!,:SiO2],gdfiltered[!,PLOT_VARS[i,j]],text = gdfiltered[!,:Sample],fontsize=6)
            end 

            if i != 3
                hidexdecorations!(ax,ticks = false,minorticks=false)
            else
                ax.xlabel = "SiO2"
            end
        end
    end
    for i in 1:2
        colsize!(fig.layout,i,Aspect(1,1.0))
    end
    rowgap!(fig.layout,0)
    save(joinpath(savedir,"variationdiagrams_labelled.svg"),fig)

end


variationdiagrams("../../Geochem/FluidFlux_BRC.csv","TernPlots/")


