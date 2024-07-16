using CairoMakie
using Isoplot
using DataFrames
using CSV
using Statistics
using FileIO
using NativeFileDialog
using JLD2
using SkipNan
include("PlotDefaults.jl")

const FIG_SIZE = (1200,900)

mutable struct TransectData
    teData::DataFrame
    UPbData::DataFrame
    UPbAvgs::DataFrame
    teAvgs::DataFrame
    teNorm::DataFrame
    elemList::Array{String}
end

TransectData() = TransectData(DataFrame(),DataFrame(),DataFrame(),DataFrame(),DataFrame(),Array{String}([]))

function replot(teDF::DataFrame,upbDF::DataFrame,selectionDF::DataFrame,normDF::DataFrame,normname,scanspeed,t1,t2,savedir;xtranslate = 0.0, xstretch = 0.0)
    instanceData = TransectData()
    instanceData.teData = teDF
    elemList = Array{String}([])
    for i=3:lastindex(names(teDF))
        elemName = match(r"[A-z]+",names(instanceData.teData)[i]).match
        push!(elemList,elemName)
    end
    instanceData.elemList = elemList

    instanceData.teData[!,:Distance] = instanceData.teData[!,"Elapsed Time"]*scanspeed

    instanceData.UPbData = DataFrame(upbDF)
    instanceData.UPbData[!,:Distance] = instanceData.UPbData[!,"Elapsed Time"]*scanspeed
    stretch = 1 + xstretch/100
    instanceData.UPbData[!,:DistanceMod] = instanceData.UPbData[!,:Distance].*stretch .+ xtranslate

    sFig, sAx = initSpiderFig(normname)
    isoFig, isoAx = initIsoplotFig(t1 = t1, t2 = t2)


    instanceData.UPbAvgs = DataFrame(Selection = Array{String}([]), x1 = Array{Float64}([]),x2=Array{Float64}([]),Analyses = Array{UPbAnalysis{Float64}}([]),
                                Pb207U235=Array{Float64}([]), Pb207U235_se1=Array{Float64}([]), Pb206U238=Array{Float64}([]), Pb206U238_se1=Array{Float64}([]),
                                Pb207Pb206=Array{Float64}([]),Pb207Pb206_se1=Array{Float64}([]),cor_0735_0638=Array{Float64}([]),cor_0706_0638=Array{Float64}([]))

    instanceData.teAvgs = DataFrame(Selection = Array{String}([]), x1 = Array{Float64}([]),x2=Array{Float64}([]))
    for elem in instanceData.elemList
        instanceData.teAvgs[!,Symbol(elem*"_mean")] = Array{Float64}([])
        instanceData.teAvgs[!,Symbol(elem*"_1std")] = Array{Float64}([])
    end
   
    instanceData.teNorm = DataFrame(Selection = Array{String}([]), x1 = Array{Float64}([]),x2=Array{Float64}([]))
    for elem in normDF[!,:Element]
        instanceData.teNorm[!,elem] = Array{Float64}([])
    end


    
    
    for row in eachrow(selectionDF)
        sIndex = parse(Int,match(r"\d+",row[:Selection]).match)
        x1 = row[:x1]
        x2 = row[:x2]
        elemAvg, teRow = averageLineTE(instanceData,start = x1,stop=x2)
        
        push!(instanceData.teAvgs,["S"*string(sIndex);teRow])
        elemRatio = plotSpider!(elemAvg,normDF,sAx,sFig,sIndex)
        push!(instanceData.teNorm,["S"*string(sIndex);x1;x2;elemRatio])
        

        sel = "S"*string(sIndex)
        selAnal, selRaw = averageUPb(instanceData,start=x1,stop=x2)

        plotUPb!(isoAx,selAnal,sIndex)
        push!(instanceData.UPbAvgs,[sel;x1;x2;selAnal;selRaw])
        xmin, xmax, ymin, ymax = datalimits(instanceData.UPbAvgs[!,:Analyses])
        limits!(isoAx,xmin,xmax,ymin,ymax)
    end

    if isdir(savedir)
        CSV.write(savedir*"/TraceElementData.csv",instanceData.teData)
        CSV.write(savedir*"/UPbData.csv",instanceData.UPbData)
        CSV.write(savedir*"/TraceElementsAvgs.csv",instanceData.teAvgs)
        CSV.write(savedir*"/TraceElementsNormalized.csv",instanceData.teNorm)
        CSV.write(savedir*"/UPbAvgs.csv",instanceData.UPbAvgs)
        CairoMakie.save(savedir*"/Spider.svg",sFig)
        CairoMakie.save(savedir*"/Isoplot.svg",isoFig)
    end
end

function NaNIfNotPositive(x::Number)
    if x <= 0
        return NaN
    end

    return x
end

function averageLineTE(instanceData::TransectData;start::AbstractFloat=0.0, stop::AbstractFloat=-1.0)

    if stop < 0
        stop = instanceData.teData[end,:Distance]
    end

    averageElems = DataFrame(Element=[],Concentration=[],Error=[])
    teRow = [start,stop]
    selection = filter([:Distance]=> x -> start<=x<=stop,instanceData.teData)
    for i=3:lastindex(names(selection))
        elem = match(r"[A-z]+",names(selection)[i]).match
        avgConc = mean(selection[!,i])
        stdDev = std(selection[!,i])
        if elem != "Distance"
            push!(averageElems,(elem,avgConc,stdDev))
            teRow = [teRow;avgConc;stdDev]
        end
    end
    
  
    return averageElems, teRow
end

function averageUPb(instanceData::TransectData;start::AbstractFloat=0.0,stop::AbstractFloat = -1.0)

    if stop < 0
        stop = instanceData.UPbData[end,:Distance]
    end

    
    selection = filter([:DistanceMod]=> x -> start<=x<=stop,instanceData.UPbData)
    Pb206U238 = mean(selection[!,Symbol("Final Pb206/U238")])
    Pb206U238_sd1 = std(selection[!,Symbol("Final Pb206/U238")])
    Pb206U238_se1 = Pb206U238_sd1/sqrt(nrow(selection))
    
    Pb207U235 = mean(selection[!,Symbol("Final Pb207/U235")])
    Pb207U235_sd1 = std(selection[!,Symbol("Final Pb207/U235")])
    Pb207U235_se1 = Pb207U235_sd1/sqrt(nrow(selection))

    Pb207Pb206 = mean(selection[!,Symbol("Final Pb207/Pb206")])
    Pb207Pb206_sd1 = std(selection[!,Symbol("Final Pb207/Pb206")])
    Pb207Pb206_se1 = Pb207Pb206_sd1/sqrt(nrow(selection))

    cor_0706_0638 = cov(selection[!,Symbol("Final Pb207/Pb206")],selection[!,Symbol("Final Pb206/U238")])/(Pb207Pb206_sd1*Pb206U238_sd1)
    cor_0735_0638 = cov(selection[!,Symbol("Final Pb207/U235")],selection[!,Symbol("Final Pb206/U238")])/(Pb207U235_sd1*Pb206U238_sd1)
    
    return UPbAnalysis(Pb207U235, Pb207U235_se1, Pb206U238, Pb206U238_se1, cor_0735_0638), 
            [Pb207U235, Pb207U235_se1, Pb206U238, Pb206U238_se1,Pb207Pb206,Pb207Pb206_se1,cor_0735_0638,cor_0706_0638] #Add Pb-Pb data with update of Isoplot.jl
end


function initSpiderFig(normName::String)
    fig = Figure(size = FIG_SIZE)
    ax = Axis(fig[1,1],
            yticksmirrored=true,
            yscale = log10,
            xminorticksize = 0,
            xminortickwidth = 0)
    ax.ylabel = "sample/"*normName
    
    return fig,ax
end

function plotSpider!(sampleData::DataFrame,comparisonData::DataFrame,ax::Axis,fig::Figure,selectionIndex::Integer)

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
    # push!(instanceData.teNorm,elemRatio)

    lines!(ax,x,NaNIfNotPositive.(elemRatio),linewidth = 4, color = Cycled(selectionIndex), label = "S"*string(selectionIndex))
    
    # lines!(ax,x,NaNIfNotPositive.(elemRatio + elemError),linewidth = 2,linestyle = :dash, color = Cycled(selectionIndex))

    # lines!(ax,x,NaNIfNotPositive.(elemRatio - elemError),linewidth = 2,linestyle = :dash, color = Cycled(selectionIndex))
    ax.xticks = (1:lastindex(elemList),elemList)
    [delete!(leg) for leg in fig.content if leg isa Legend]
    fig[1,2] = Legend(fig,ax,framevisible = false)
    return elemRatio
end

function initIsoplotFig(;showConcordia::Bool = true, t1::Number = 0, t2::Number = Isoplot.t🜨)
    fig = Figure(size = FIG_SIZE)
    ax = Axis(fig[1,1])

    if showConcordia
        concordiacurve!(ax,t1,t2,fontsize = 15)
    end
    ax.xlabel =rich(superscript("207"),"Pb/",superscript("235"),"U")
    ax.ylabel =rich(superscript("206"),"Pb/",superscript("238"),"U")
    return fig,ax
end

function plotUPb!(ax::Axis,analysis::UPbAnalysis,selectionIndex::Integer)

    
    plot!(ax,analysis, color = Cycled(selectionIndex),alpha = 0.6)
    
    text!(ax,analysis.μ[1],analysis.μ[2],text="S"*string(selectionIndex),fontsize = 15)
end

println("Pick TE dataset")
newTE = pick_file(;filterlist="csv")
println("Pick UPb dataset")
newUPb = pick_file(;filterlist="csv")
println("Pick previously exported datasets")
oldExports = pick_file(;filterlist="csv")
println("Pick folder to save")
savedir =  pick_folder()
normFile = "Chondrite.csv"

set_theme!(myTheme)
teDF = DataFrame(CSV.File(newTE))
upbDF = DataFrame(CSV.File(newUPb))
selDF = DataFrame(CSV.File(oldExports))
upbDF = DataFrame(CSV.File(newUPb))
normDF = DataFrame(CSV.File(normFile))

println("What is the x-offset (μm)?")
offset = parse(Float64,readline())
println("What is the stretch factor (%)?")
stretch = parse(Float64,readline())

replot(teDF,upbDF,selDF,normDF,"chondrite",2,800, 1600,savedir;xtranslate = offset, xstretch = stretch)