using CairoMakie
#using AlgebraOfGraphics
# using MakieCore
using DataFrames
using CSV
using Measurements
using Statistics
myColours = ["#FF8CA6","#6D8FFF","#F9C52A","#A9EFA5", "#E6A7FB","#000000","#FF0000","#008000","#0000FF","#4B0082","#FF8C00"]
REES = ["La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"]

# MakieCore.convert_arguments(x::Vector{String},y::Vector{Float32}) = 
function plotElem(lineData::DataFrame, scanSpeed::Int, elem::String,ax::Axis,fig::Figure,fileName::String)
    empty!(ax)
    lineData[!,:Distance] = lineData[!,Symbol("Elapsed Time")]*scanSpeed
    x = Vector(lineData[!,:Distance])
    y1 = DataFrame(lineData[!,Regex(elem)])
    y2 = Matrix(y1)
    y = vec(y2)
    ax.xlabel = "Distance (μm)"
    ax.ylabel = elem*" (μg/g)"
  
    xMax = round(maximum(lineData[!,:Distance]),sigdigits=2)
    xInterval = 20
    if xMax <120
        xInterval = 10
        if xMax <= 60
            xInterval = 5
            if xMax <= 25
                xInterval = 2
            end
        end
    end    

    ax.xticks = 0:xInterval:xMax
    lines!(x,y,color = myColours[1])
    save("ZrnLines//"*fileName*".svg",fig)
end

function plotRatio(lineData::DataFrame, scanSpeed::Int, elemN::String,elemD::String,ax::Axis,fig::Figure,fileName::String)
    empty!(ax)
    lineData[!,:Distance] = lineData[!,Symbol("Elapsed Time")]*scanSpeed
    x = Vector(lineData[!,:Distance])
    y1a = DataFrame(lineData[!,Regex(elemN)])
    y1b = Matrix(y1a)
    y1c = vec(y1b)
    y2a = DataFrame(lineData[!,Regex(elemD)])
    y2b = Matrix(y2a)
    y2c = vec(y2b)

    y = y1c./y2c
    ax.xlabel = "Distance (μm)"
    ax.ylabel = elemN*"/"*elemD
  
    xMax = round(maximum(lineData[!,:Distance]),sigdigits=2)
    xInterval = 20
    if xMax <120
        xInterval = 10
        if xMax <= 60
            xInterval = 5
            if xMax <= 25
                xInterval = 2
            end
        end
    end    

    ax.xticks = 0:xInterval:xMax
    lines!(x,y,color = myColours[1],linewidth=2)
    save("ZrnLines//"*fileName*".svg",fig)

end

function NaNIfNegative(x::Number)
    if x < 0
        return NaN
    end

    return x
end

#SampleData should be averaged with format: Element, Concentration, Error
#Comparison will be e.g. Chondrite concentrations
#Will try to match the corresponding elements
function plotSpider(sampleData::DataFrame,comparisonData::DataFrame,normName::String,ax::Axis,fig::Figure,fileName::String)

    x = Array{Int64}([])
    elemList = Array{String}([])
    elemRatio = Array{Float64}([])
    elemError = Array{Float64}([])
    count = 1
    for el in eachrow(comparisonData)
        for em in eachrow(sampleData)
            if el[:Element] == em[:Element]
                push!(x,count)
                push!(elemList,el[:Element])
                push!(elemRatio,em[:Concentration]/el[:Concentration])
                push!(elemError,em[:Error]/el[:Concentration])
                count += 1
            end
        end
    end

    
    empty!(ax)
    ax.ylabel ="sample/"*normName
    ax.yscale = log10
    lines!(x,elemRatio,linewidth = 2)
    # errorbars!(x,elemRatio,elemError)
    lines!(x,elemRatio + elemError,linewidth = 1,linestyle = :dash)

    lines!(x,NaNIfNegative.(elemRatio - elemError),linewidth = 1,linestyle = :dash)
    ax.xticks = (1:lastindex(elemList),elemList)
    
    # ax.ymin = -25
    # ax.ymax = 225
    save("SpiderPlots//"*fileName*".svg",fig)
    
end

#Takes the average across the given interval (or across the whole interval)
function averageLineData(lineData::DataFrame;start::Float64=0.0, stop::Float64=-1.0)

    if stop < 0
        stop = lineData(lastindex(lineData),Symbol("Elapsed Time"))
    end

    averageElems = DataFrame(Element=[],Concentration=[],Error=[])
    selection = filter([Symbol("Elapsed Time")]=> time -> start<=time<=stop,lineData)
    for i=3:lastindex(names(selection))
        elem = match(r"[A-z]+",names(selection)[i]).match
        avgConc = mean(selection[!,i])
        stdDev = std(selection[!,i])
        push!(averageElems,(elem,avgConc,stdDev))
    end

    return averageElems
end
fig = Figure(font = "B612",fontsize = 18)
ax = Axis(fig[1,1],
    xlabelsize = 28,
    ylabelsize = 28,
    xminorticksvisible = true,
    yminorticksvisible = true,
    xticksmirrored = true,
    yticksmirrored = true,
    xtickalign = 1,
    ytickalign = 1,
    xminortickalign = 1,
    yminortickalign = 1,
    xticksize = 10,
    xtickwidth = 2,
    yticksize = 10,
    ytickwidth = 2,
    xminorticksize = 5,
    xminortickwidth = 2,
    yminorticksize = 5,
    yminortickwidth = 2,
    xgridvisible = false,
    ygridvisible = false
)
# ax.xlabel = "Distance (μm)"

zrn37 = DataFrame(CSV.File(raw"C:\Users\Sabas\OneDrive - University of Waterloo\Documents\Waterloo\LAICPMS\2023_06_22\20SD06-2_Zrn37.csv"))
zrn41 = DataFrame(CSV.File(raw"C:\Users\Sabas\OneDrive - University of Waterloo\Documents\Waterloo\LAICPMS\2023_06_22\20SD06-2_Zrn41.csv"))
zrn10 = DataFrame(CSV.File(raw"C:\Users\Sabas\OneDrive - University of Waterloo\Documents\Waterloo\LAICPMS\2023_06_22\20SD06-2_Zrn10.csv"))
# ax.yticks = 0:25:150
# ax.limits = (nothing, nothing, 0,150)

# plotElem(zrn37,2,"Yb",ax,fig,"Zrn37_Yb")
# plotElem(zrn41,2,"Yb",ax,fig,"Zrn41_Yb")
# plotElem(zrn10,2,"Yb",ax,fig,"Zrn10_Yb")

# ax.yticks = 0:10:100
# ax.limits = (nothing, nothing, 0,100)
# plotElem(zrn37,2,"U",ax,fig,"Zrn37_U")
# plotElem(zrn41,2,"U",ax,fig,"Zrn41_U")
# plotElem(zrn10,2,"U",ax,fig,"Zrn10_U")

# plotElem(zrn37,2,"Hf",ax,fig,"Zrn37_Hf")
# ax.yticks = 0:10:100
# ax.limits = (nothing, nothing, 0,100)
# plotRatio(zrn37,2,"Lu","Hf",ax,fig,"Zrn37_LuHf")
# plotRatio(zrn41,2,"Lu","Hf",ax,fig,"Zrn41_LuHf")
# plotRatio(zrn10,2,"Lu","Hf",ax,fig,"Zrn10_LuHf")
# plotRatio(zrn37,2,"Th","U",ax,fig,"Zrn37_ThU")
# plotRatio(zrn41,2,"Th","U",ax,fig,"Zrn41_ThU")
# fig
ax.limits = (nothing, nothing, 0.001,1000)
zrn37rim = averageLineData(zrn37,start=85.0,stop=120.0)
zrn37core1 = averageLineData(zrn37,start=40.0,stop=80.0)
chondrite = DataFrame(CSV.File("Chondrite.csv"))
plotSpider(zrn37rim,chondrite,"chondrite",ax,fig,"Zrn37rim")
plotSpider(zrn37core1,chondrite,"chondrite",ax,fig,"Zrn37core1")