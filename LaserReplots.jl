# using CairoMakie
using GLMakie
using Isoplot
using DataFrames
using CSV
using Statistics
using LinearRegression
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
    # @show savedir
    if isdir(savedir)
        CSV.write(savedir*"/TraceElementData.csv",instanceData.teData)
        CSV.write(savedir*"/UPbData.csv",instanceData.UPbData)
        CSV.write(savedir*"/TraceElementsAvgs.csv",instanceData.teAvgs)
        CSV.write(savedir*"/TraceElementsNormalized.csv",instanceData.teNorm)
        CSV.write(savedir*"/UPbAvgs.csv",instanceData.UPbAvgs)
        CairoMakie.save(savedir*"/Spider.png",sFig)
        CairoMakie.save(savedir*"/Isoplot.png",isoFig)
    else
        throw(ArgumentError("$savedir not found"))
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

function r_squared(y,y_model)
    return 1 - sum((y .- y_model).^2)/sum((y.-mean(y)).^2)
end

function rescale_UTh!(teData, upbData)
    #First adjust x slightly so that they match between UPb and TEs
    #Then perform linear regression using U_approx as independent and TEs as dependent
    #Scaling factor will be the equation of the line
    xTE = teData[!,:Distance]
    xUPb = upbData[!,:DistanceMod]
    U_TE_match = Array{Float64}([])
    Th_TE_match = Array{Float64}([])
    weights = Array{Float64}([])
    for i in 1:lastindex(xUPb)
        
        xdiff = abs(xUPb[i]-xTE[1])
        j = 2
        xdiffnext = abs(xUPb[i]-xTE[j])
        while xdiffnext < xdiff && j < lastindex(xTE)
            xdiff = xdiffnext
            j += 1
            xdiffnext = abs(xUPb[i]-xTE[j])
        end

        U_TE_match = push!(U_TE_match,teData[j,Regex("U\\d+")][1])
        Th_TE_match = push!(Th_TE_match,teData[j,Regex("Th\\d+")][1])

        if teData[j,Regex("U\\d+")][1] == 0
            push!(weights, 50.0)
        else
            push!(weights, 1.0)
        end
    end

    @show length(U_TE_match)
    regression = linregress(upbData[!,:Approx_U_PPM],U_TE_match,weights)
    slope = LinearRegression.slope(regression)
    y_intercept = LinearRegression.bias(regression)
    
    upbData[!,:U_rescaled] = upbData[!,:Approx_U_PPM].*slope .+ y_intercept
    println("R^2 U = " *string(r_squared(U_TE_match,upbData[!,:U_rescaled])))
    upbData[!,:Th_rescaled] = upbData[!,:U_rescaled]./upbData[!,"Final U/Th"]
    println("R^2 Th = " *string(r_squared(Th_TE_match,upbData[!,:Th_rescaled])))

end
# println("Pick TE dataset")
# newTE = pick_file(;filterlist="csv")
# println("Pick UPb dataset")
# newUPb = pick_file(;filterlist="csv")
# println("Pick previously exported datasets")
# oldExports = pick_file(;filterlist="csv")
# println("Pick folder to save")
# savedir =  pick_folder()
# normFile = "Chondrite.csv"

set_theme!(myTheme)
# teDF = DataFrame(CSV.File(newTE))
# upbDF = DataFrame(CSV.File(newUPb))
# selDF = DataFrame(CSV.File(oldExports))
# normDF = DataFrame(CSV.File(normFile))

# println("What is the x-offset (μm)?")
# offset = parse(Float64,readline())
# println("What is the stretch factor (%)?")
# stretch = parse(Float64,readline())

# replot(teDF,upbDF,selDF,normDF,"chondrite",2,800, 1600,savedir;xtranslate = offset, xstretch = stretch)


newTEs = []
newUPbs = []
oldExports = []
savedirs = []
normFile = "chondrite.csv"
offsets = []
stretches = []
normDF = DataFrame(CSV.File(normFile))

# push!(newTEs, "ZrnLaserPlots/20SD06/Zrn10/20SD06-2_Zrn10_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD06/Zrn10/20SD06_Zrn10_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD06/Zrn10/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD06/Zrn10/Exports_91500-com")
# push!(offsets,0.0)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/20SD06/Zrn26/20SD06-2_Zrn26_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD06/Zrn26/20SD06_Zrn26_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD06/Zrn26/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD06/Zrn26/Exports_91500-com")
# push!(offsets,-2.3)
# push!(stretches,2.1)

# push!(newTEs, "ZrnLaserPlots/20SD06/Zrn27/20SD06-2_Zrn27_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD06/Zrn27/20SD06_Zrn27_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD06/Zrn27/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD06/Zrn27/Exports_91500-com")
# push!(offsets,2.8)
# push!(stretches,10.4)

# push!(newTEs, "ZrnLaserPlots/20SD06/Zrn33/20SD06-2_Zrn33_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD06/Zrn33/20SD06_Zrn33_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD06/Zrn33/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD06/Zrn33/Exports_91500-com")
# push!(offsets,-3.7)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/20SD06/Zrn37/20SD06-2_Zrn37_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD06/Zrn37/20SD06_Zrn37_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD06/Zrn37/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD06/Zrn37/Exports_91500-com")
# push!(offsets,4.4)
# push!(stretches,2.7)


# push!(newTEs, "ZrnLaserPlots/20SD06/Zrn43/20SD06-2_Zrn43_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD06/Zrn38/20SD06_Zrn38_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD06/Zrn38/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD06/Zrn38/Exports_91500-com")
# push!(offsets,-2.6)#Offset to match old profile
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/20SD06/Zrn40/20SD06-2_Zrn40_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD06/Zrn40/20SD06_Zrn40_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD06/Zrn40/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD06/Zrn40/Exports_91500-com")
# push!(offsets,-2.7)
# push!(stretches,0.6)

# push!(newTEs, "ZrnLaserPlots/20SD06/Zrn41/20SD06-2_Zrn41_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD06/Zrn41/20SD06_Zrn41_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD06/Zrn41/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD06/Zrn41/Exports_91500-com")
# push!(offsets,5.3)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/20SD06/Zrn43/20SD06-2_Zrn43_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD06/Zrn43/20SD06_Zrn43_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD06/Zrn43/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD06/Zrn43/Exports_91500-com")
# push!(offsets,-13.3)
# push!(stretches,1.1)

# push!(newTEs, "ZrnLaserPlots/20SD06/Zrn44/20SD06-2_Zrn44_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD06/Zrn44/20SD06_Zrn44_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD06/Zrn44/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD06/Zrn44/Exports_91500-com")
# push!(offsets,-5.1)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/22SD55E/Zrn14/22SD55E_Zrn14a_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/22SD55E/Zrn14/22SD55E_Zrn14a_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/22SD55E/Zrn14/14aExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/22SD55E/Zrn14/14aExports_91500-com")
# push!(offsets,-0.8)
# push!(stretches,1.0)

# push!(newTEs, "ZrnLaserPlots/22SD55E/Zrn14/22SD55E_Zrn14b_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/22SD55E/Zrn14/22SD55E_Zrn14b_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/22SD55E/Zrn14/14bExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/22SD55E/Zrn14/14bExports_91500-com")
# push!(offsets,-11.8)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/22SD55E/Zrn15/22SD55E_Zrn15_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/22SD55E/Zrn15/22SD55E_Zrn15_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/22SD55E/Zrn15/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/22SD55E/Zrn15/Exports_91500-com")
# push!(offsets,4.7)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/22SD55E/Zrn16/22SD55E_Zrn16a_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/22SD55E/Zrn16/22SD55E_Zrn16a_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/22SD55E/Zrn16/16aExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/22SD55E/Zrn16/16aExports_91500-com")
# push!(offsets,-0.8)
# push!(stretches,12.8)

# push!(newTEs, "ZrnLaserPlots/22SD55E/Zrn16/22SD55E_Zrn16b_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/22SD55E/Zrn16/22SD55E_Zrn16b_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/22SD55E/Zrn16/16bExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/22SD55E/Zrn16/16bExports_91500-com")
# push!(offsets,-0.2)
# push!(stretches,10.2)

# push!(newTEs, "ZrnLaserPlots/22SD55E/Zrn22/22SD55E_Zrn22a_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/22SD55E/Zrn22/22SD55E_Zrn22a_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/22SD55E/Zrn22/22aExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/22SD55E/Zrn22/22aExports_91500-com")
# push!(offsets,6.8)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/22SD55E/Zrn22/22SD55E_Zrn22b_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/22SD55E/Zrn22/22SD55E_Zrn22b_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/22SD55E/Zrn22/22bExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/22SD55E/Zrn22/22bExports_91500-com")
# push!(offsets,0.0)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/22SD55E/Zrn24/22SD55E_Zrn24_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/22SD55E/Zrn24/22SD55E_Zrn24_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/22SD55E/Zrn24/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/22SD55E/Zrn24/Exports_91500-com")
# push!(offsets,-9.5)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/22SD55E/Zrn27/22SD55E_Zrn27_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/22SD55E/Zrn27/22SD55E_Zrn27_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/22SD55E/Zrn27/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/22SD55E/Zrn27/Exports_91500-com")
# push!(offsets,-0.8)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/22SD55E/Zrn30/22SD55E_Zrn30_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/22SD55E/Zrn30/22SD55E_Zrn30_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/22SD55E/Zrn30/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/22SD55E/Zrn30/Exports_91500-com")
# push!(offsets,0.0)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/20SD17A/Zrn4/20SD17A_Zrn4a_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD17A/Zrn4/20SD17A_Zrn4a_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD17A/Zrn4/4aExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD17A/Zrn4/4aExports_91500-com")
# push!(offsets,-3.6)
# push!(stretches,1.3)

# push!(newTEs, "ZrnLaserPlots/20SD17A/Zrn4/20SD17A_Zrn4b_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD17A/Zrn4/20SD17A_Zrn4b_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD17A/Zrn4/4bExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD17A/Zrn4/4bExports_91500-com")
# push!(offsets,11.3)
# push!(stretches,-2.4)

# push!(newTEs, "ZrnLaserPlots/20SD17A/Zrn6/20SD17A_Zrn6_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD17A/Zrn6/20SD17A_Zrn6_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD17A/Zrn6/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD17A/Zrn6/Exports_91500-com")
# push!(offsets,-5.3)
# push!(stretches,1.1)

# push!(newTEs, "ZrnLaserPlots/20SD17A/Zrn14/20SD17A_Zrn14_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD17A/Zrn14/20SD17A_Zrn14_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD17A/Zrn14/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD17A/Zrn14/Exports_91500-com")
# push!(offsets,-11.0)
# push!(stretches,2.0)

# push!(newTEs, "ZrnLaserPlots/20SD17A/Zrn15/20SD17A_Zrn15_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD17A/Zrn15/20SD17A_Zrn15_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD17A/Zrn15/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD17A/Zrn15/Exports_91500-com")
# push!(offsets,-2.9)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/20SD17A/Zrn14/20SD17A_Zrn14_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD17A/Zrn23/20SD17A_Zrn23_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD17A/Zrn23/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD17A/Zrn23/Exports_91500-com")
# push!(offsets,-2.4)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/20SD17A/Zrn25/20SD17A_Zrn25_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD17A/Zrn25/20SD17A_Zrn25_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD17A/Zrn25/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD17A/Zrn25/Exports_91500-com")
# push!(offsets,-3.2)
# push!(stretches,0.0)


# push!(newTEs, "ZrnLaserPlots/20SD17A/Zrn26/20SD17A_Zrn26_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD17A/Zrn26/20SD17A_Zrn26_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD17A/Zrn26/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD17A/Zrn26/Exports_91500-com")
# push!(offsets,4.0)
# push!(stretches,-4.7)

# push!(newTEs, "ZrnLaserPlots/20SD17A/Zrn28/20SD17A_Zrn28a_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD17A/Zrn28/20SD17A_Zrn28a_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD17A/Zrn28/28aExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD17A/Zrn28/28aExports_91500-com")
# push!(offsets,-1.2)
# push!(stretches,7.4)

# push!(newTEs, "ZrnLaserPlots/20SD17A/Zrn28/20SD17A_Zrn28b_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/20SD17A/Zrn28/20SD17A_Zrn28b_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/20SD17A/Zrn28/28bExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/20SD17A/Zrn28/28bExports_91500-com")
# push!(offsets,-6.6)
# push!(stretches,-1.0)

# push!(newTEs, "ZrnLaserPlots/21SD68/Zrn4/21SD68_Zrn4_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/21SD68/Zrn4/21SD68_Zrn4_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/21SD68/Zrn4/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/21SD68/Zrn4/Exports_91500-com")
# push!(offsets,6.6)
# push!(stretches,0.0)


# push!(newTEs, "ZrnLaserPlots/21SD68/Zrn6/21SD68_Zrn6_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/21SD68/Zrn6/21SD68_Zrn6_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/21SD68/Zrn6/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/21SD68/Zrn6/Exports_91500-com")
# push!(offsets,-12.0)
# push!(stretches,-0.7)


# push!(newTEs, "ZrnLaserPlots/21SD68/Zrn7/21SD68_Zrn7a_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/21SD68/Zrn7/21SD68_Zrn7a_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/21SD68/Zrn7/7aExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/21SD68/Zrn7/7aExports_91500-com")
# push!(offsets,-10.6)
# push!(stretches,-3.5)

# push!(newTEs, "ZrnLaserPlots/21SD68/Zrn7/21SD68_Zrn7b_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/21SD68/Zrn7/21SD68_Zrn7b_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/21SD68/Zrn7/7bExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/21SD68/Zrn7/7bExports_91500-com")
# push!(offsets,-4.0)
# push!(stretches,-2.0)

# push!(newTEs, "ZrnLaserPlots/21SD68/Zrn7/21SD68_Zrn7a_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/21SD68/Zrn9/21SD68_Zrn9_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/21SD68/Zrn9/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/21SD68/Zrn9/Exports_91500-com")
# push!(offsets,-2.2)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/21SD68/Zrn17/21SD68_Zrn17_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/21SD68/Zrn17/21SD68_Zrn17_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/21SD68/Zrn17/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/21SD68/Zrn17/Exports_91500-com")
# push!(offsets,-2.2)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/21SD68/Zrn22/21SD68_Zrn22_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/21SD68/Zrn22/21SD68_Zrn22_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/21SD68/Zrn22/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/21SD68/Zrn22/Exports_91500-com")
# push!(offsets,1.5)
# push!(stretches,0.0)

# push!(newTEs, "ZrnLaserPlots/21SD68/Zrn23/21SD68_Zrn23_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/21SD68/Zrn23/21SD68_Zrn23_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/21SD68/Zrn23/Exports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/21SD68/Zrn23/Exports_91500-com")
# push!(offsets,-16.2)
# push!(stretches,-1.5)

# push!(newTEs, "ZrnLaserPlots/21SD68/Zrn27/21SD68_Zrn27a_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/21SD68/Zrn27/21SD68_Zrn27a_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/21SD68/Zrn27/27aExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/21SD68/Zrn27/27aExports_91500-com")
# push!(offsets,-8.6)
# push!(stretches,-2.6)

# push!(newTEs, "ZrnLaserPlots/21SD68/Zrn27/21SD68_Zrn27b_TE.csv")
# push!(newUPbs, "ZrnLaserPlots/21SD68/Zrn27/21SD68_Zrn27b_UPb_91500-com.csv")
# push!(oldExports,"ZrnLaserPlots/21SD68/Zrn27/27bExports/UPbAvgs.csv")
# push!(savedirs,"ZrnLaserPlots/21SD68/Zrn27/27bExports_91500-com")
# push!(offsets,-12.8)
# push!(stretches,0.0)

# for i in 1:lastindex(newTEs)
#     teDF = DataFrame(CSV.File(newTEs[i]))
#     upbDF = DataFrame(CSV.File(newUPbs[i]))
#     selDF = DataFrame(CSV.File(oldExports[i]))
    
#     replot(teDF,upbDF,selDF,normDF,"chondrite",2,800, 1600,savedirs[i];xtranslate = offsets[i], xstretch = stretches[i])
# end
tefile = "ZrnLaserPlots/20SD06/Zrn41/Exports_91500-com/TraceElementData.csv"
upbfile = "ZrnLaserPlots/20SD06/Zrn41/Exports_91500-com/UPbData.csv"
teDF = DataFrame(CSV.File(tefile))
upbDF= DataFrame(CSV.File(upbfile))

rescale_UTh!(teDF, upbDF)

GLMakie.activate!()
fig = Figure(size = FIG_SIZE)

ax = Axis(fig[1,1])

lines!(ax,teDF[!,:Distance],teDF[!,"U238_ppm"])
lines!(ax,upbDF[!,:DistanceMod],upbDF[!,:U_rescaled])
display(GLMakie.Screen(),fig)

fig2 = Figure(size = FIG_SIZE)
ax2 = Axis(fig2[1,1])
lines!(ax2,teDF[!,:Distance],teDF[!,"Th232_ppm"])
lines!(ax2,upbDF[!,:DistanceMod],upbDF[!,:Th_rescaled])
display(GLMakie.Screen(),fig2)