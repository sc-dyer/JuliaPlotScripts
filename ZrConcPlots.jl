using CairoMakie
using DataFrames
using CSV
include("PlotDefaults.jl")

function plotZrConc(fileName::String, fig::Figure,label)


    ax = Axis(fig,yscale = log10,
    xminorticksize = 0,
    xminortickwidth = 0,
    yminorticksize = 0,
    yminortickwidth =0)
    ylims!(ax,10E-4,10E3)#Not sure why this makes limits 10E-3 and 10E3
    zrData = DataFrame(CSV.File(fileName))

    ax.xticks = (1:lastindex(names(zrData)),names(zrData))
    
    x = Array{Int64}([])
    y = Array{Float64}([])
    for i=1:lastindex(names(zrData))
        yAdd = collect(skipmissing(zrData[!,i]))
        xAdd = [i for dum = 1:lastindex(yAdd)]
        nvals = length(yAdd)
        
        x = vcat(x,xAdd)
        y = vcat(y,yAdd)
        if nvals > 0
            maxY = maximum(yAdd)
            text!(ax,i,maximum(maxY),text="n = "*string(nvals),color=:black,fontsize=14,offset = (0,3), align = (:center,:bottom))
        end
        # #An attempt to add empty columns for minerals with values BDL
        # if length(yAdd) == 0
        #     x = vcat(x,i)
        #     y = vcat(y,0.0)
        # end
    end

    boxplot!(ax,x,y,color = myColours[1])#color = x, colormap = myColours)#, color = myColours[1])
    # ax.ylabel = "Average Zr concentration (μg/g)"

    text!(ax, 0.05,0.87, text = label,fontsize=36,space = :relative,offset = (4, 0))
    return ax

end


set_theme!(myTheme)
update_theme!(;aspect = 1, size = (1000,1000))
# update_theme!(;Axis = (;ylims = (10E-2,10E3)))
#                     xminorticksize = 0,
#                     xminortickwidth = 0)
fig = Figure()
sam20SD06 = plotZrConc("../../LAICPMS/2023_09_20/20SD06_ZrTable.csv",fig,"20SD06")

sam20SD17A = plotZrConc("../../LAICPMS/2023_09_20/20SD17A_ZrTable.csv",fig,"20SD17A")
sam20SD06.ylabel = "Average Zr concentration (μg/g)"
hideydecorations!(sam20SD17A, ticks = false, minorticks = false)
sam21SD68 = plotZrConc("../../LAICPMS/2023_09_20/21SD68_ZrTable.csv",fig,"21SD68")
sam22SD55E = plotZrConc("../../LAICPMS/2023_09_20/22SD55E_ZrTable.csv",fig,"22SD55E")
sam21SD68.ylabel = "Average Zr concentration (μg/g)"
hideydecorations!(sam22SD55E, ticks = false, minorticks = false)
fig[1,1] = sam20SD06
fig[1,2] = sam20SD17A
fig[2,1] = sam21SD68
fig[2,2] = sam22SD55E
labels = ["20SD06" "20SD17A";"21SD68" "22SD55E"]

CairoMakie.save("ZrPlots/ZrConc.svg",fig)

# CairoMakie.save("ZrPlots//20SD17A.svg",sam20SD17A)

# CairoMakie.save("ZrPlots//21SD68.svg",sam21SD68)

# CairoMakie.save("ZrPlots//22SD55E.svg",sam22SD55E)