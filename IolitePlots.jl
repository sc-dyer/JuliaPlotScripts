using CairoMakie
using DataFrames
using CSV
using Statistics
# using FileIO
# using Gtk
colourchoice=3
include("PlotDefaults.jl")
DEFAULT_NORM = DataFrame(Element = ["Rb","Sr","Y","La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"],Concentration = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
function filterSampleName(df::DataFrame,sampleName::String)
    return filter([:Sample]=>sam -> occursin(sampleName,sam),df)
end

function toNorm(geochem::DataFrame,sampleName::String;templateNorm::DataFrame = DEFAULT_NORM)
    newNorm = copy(templateNorm)
    convertDf = filter([:Sample]=>sam -> occursin(sampleName,sam),geochem)
    for el in eachrow(newNorm)
        el[:Concentration] = convertDf[1,Regex(el[:Element]*"\$")][1]
    end

    return newNorm
end


function NaNIfNotPositive(x::Number)
    if x <= 0
        return NaN
    end

    return x
end

function isVisible(fig::Figure)
    if length(fig.scene.current_screens) > 0
        return true
    end
    return false
end

function initSpiderFig(normName::String,comparisonData::DataFrame)
    fig = Figure(size = (1200,900))
    elemList = Array{String}([])
    ax = Axis(fig[1,1],
            yticksmirrored=true,
            yscale = log10,
            xminorticksize = 0,
            xminortickwidth = 0)
    ax.ylabel = "sample/"*normName

    for el in eachrow(comparisonData)
        push!(elemList,el[:Element])
    end
    ax.xticks = (1:lastindex(elemList),elemList)

    return fig,ax
end

function initSpiderFigs(normName::String,comparisonData::DataFrame,grid)
    fig = Figure(size = (1000,1000))
    elemList = Array{String}([])
    for el in eachrow(comparisonData)
        push!(elemList,el[:Element])
    end
    axs = []
    for i in 1:grid[1]
        axrow = []
        for j in 1:grid[2]
            ax = Axis(fig[i,j],
                    yticksmirrored=true,
                    yscale = log10,
                    xminorticksize = 0,
                    xminortickwidth = 0,
                    aspect = 1.0)
            ax.ylabel = "sample/"*normName
            ax.xticks = (1:lastindex(elemList),elemList)
            if length(axrow)==0
                axrow = [ax]
            else
                axrow = [axrow ax]
            end
           
        end
        if length(axs) == 0
            axs = axrow
        else
            axs = [axs; axrow]
        end
       
    end
    
    

    return fig,axs
end

#This will assume multiple records with a normal iolite dataframe layout in the sampleData
function plotSpider!(sampleData::DataFrame,comparisonData::DataFrame,ax::Axis,fig::Figure,selectionIndex::Integer,selectionLabel::String)

    # elemError = Array{Float64}([])
    eu_anom_avg = 0.0
    for sample in eachrow(sampleData)
        x = Array{Int64}([])
        elemRatio = Array{Float64}([])
        count = 1
        
        for el in eachrow(comparisonData)
            push!(x,count)
            try
                push!(elemRatio,sample[Regex(el[:Element]*"\\d+.+mean")][1]/el[:Concentration])
            catch
                push!(elemRatio,sample[el[:Element]][1]/el[:Concentration])
            end
            count += 1
        end
        # push!(elemError,em[:Error]/el[:Concentration])
        lines!(ax,x,NaNIfNotPositive.(elemRatio),linewidth = 1, color = Cycled(selectionIndex), label = selectionLabel)
        
        
       
    end
   

    # [delete!(leg) for leg in fig.content if leg isa Legend]
    
    
end

function calcEuAnomaly!(sampleData::DataFrame,normData::DataFrame)

    normSm = filter([:Element]=>norm -> occursin("Sm",norm),normData)[1,:Concentration]
    normGd = filter([:Element]=>norm -> occursin("Gd",norm),normData)[1,:Concentration]
    normEu = filter([:Element]=>norm -> occursin("Eu",norm),normData)[1,:Concentration]
    try
        sampleData[!,:Eu_EuA] = (sampleData[!,r"Eu\d+.+mean"][!,1] / normEu)./sqrt.((sampleData[!,r"Sm\d+.+mean"][!,1] / normSm).*(sampleData[!,r"Gd\d+.+mean"][!,1] / normGd))
    catch
        sampleData[!,:Eu_EuA] = (sampleData[!,r"Eu"][!,1] / normEu)./sqrt.((sampleData[!,r"Sm"][!,1] / normSm).*(sampleData[!,r"Gd"][!,1] / normGd))
    end
end

function normElem!(sampleData::DataFrame,normData::DataFrame,elem::String)
    normElem = filter([:Element]=>norm -> occursin(elem,norm),normData)[1,:Concentration]
    sampleData[!,Symbol(elem*"_norm")] = sampleData[!,Regex(elem*"\\d+.+mean")][!,1]/normElem
end

function sumREEs!(sampleData::DataFrame)

    sampleData[!,:TotalREE] = 
                sampleData[!,r"La\d+.+mean"][!,1] .+ sampleData[!,r"Ce\d+.+mean"][!,1] .+ sampleData[!,r"Pr\d+.+mean"][!,1] .+
                sampleData[!,r"Nd\d+.+mean"][!,1] .+ sampleData[!,r"Sm\d+.+mean"][!,1] .+ sampleData[!,r"Eu\d+.+mean"][!,1] .+
                sampleData[!,r"Gd\d+.+mean"][!,1] .+ sampleData[!,r"Tb\d+.+mean"][!,1] .+ sampleData[!,r"Dy\d+.+mean"][!,1] .+
                sampleData[!,r"Ho\d+.+mean"][!,1] .+ sampleData[!,r"Er\d+.+mean"][!,1] .+ sampleData[!,r"Tm\d+.+mean"][!,1] .+
                sampleData[!,r"Yb\d+.+mean"][!,1] .+ sampleData[!,r"Lu\d+.+mean"][!,1]
end

function plotsamples(df, samplelists, norm, normname)
    fig, axs = initSpiderFigs(normname,norm,[2,2])
    calcEuAnomaly!(df,norm)
    for i in 1:lastindex(samplelists)
        samples = samplelists[i]
        row = Integer(trunc(i/2)+i%2)
        col = (i+1)%2+1
       
        ax = axs[row,col]
        for j in 1:lastindex(samples)
            println(samples[j])
            sampledf = filterSampleName(df,samples[j])
            println("Eu/Eu* = ",mean(sampledf[!,:Eu_EuA]))
            plotSpider!(sampledf,norm,ax,fig,j,samples[j])
            
        end
        axislegend(ax,framevisible = false,merge = true,position = :rb)
        limits!(ax,nothing,nothing,1,3000)
    end

    return fig
end

set_theme!(myTheme)
# samples = ["21SD08-2","21SD08-3","23SD02B-1","21SD09","21SD51D","21SD51F","21SD56A","22SD13C-2","22SD55B3","22SD55D","22SD55E","20SD06","21SD68","23SD03C-1","23SD03C-3","23SD03E",
#     "23SD03H-3","23SD03J-2","23SD20B","23SD20C","23SD20D","23SD20F","23SD20G","2M02061A","2M02061B","M30072-1A","M30072-1C","2M0506-16_2","2M0506-16_3B","M2707-16B_3B","M2707-5B"]
# samplesWGeochem = ["21SD08-2","21SD08-3","23SD02B-1","21SD51D","22SD13C-2","22SD55B3","22SD55D","22SD55E","20SD06","21SD68","23SD03C-1","23SD03C-3","23SD03E",
#     "23SD03H-3","23SD03J-2","23SD20B","23SD20C","23SD20D","23SD20F","23SD20G"]
# geochemSamples = ["21SD08","21SD08","23SD02B","21SD51D","22SD13C","22SD55B","22SD55D","22SD55E","20SD06","21SD68","23SD03C-1","23SD03C-5","23SD03E","23SD03H-Bulk",
#     "23SD03J-2","23SD20B","23SD20C","23SD20D","23SD20F","23SD20G"]

# outcrop1 = ["21SD08-2","21SD08-3","23SD02B-1"]
outcrop2 = ["23SD03C-3","23SD03C-1","23SD03E"]
outcrop3 = ["22SD55E","22SD55B3","22SD55D"]
outcrop4 = ["23SD20B","23SD20C","23SD20F","23SD20G"]
migmatites = ["2M02061A","2M02061B","2M0506-16_2","2M0506-16_3B"]
norm = DataFrame(CSV.File("chondrite.csv"))
ampData = DataFrame(CSV.File("../../LAICPMS/AmphiboleData_UpdatedSi/Amphiboles.csv"))
geochemData = DataFrame(CSV.File("../../Geochem/FluidFlux_BRCb.csv"))

fig = plotsamples(ampData,[outcrop2,outcrop3,outcrop4,migmatites],norm,"chondrite")
# save("SpiderPlots/AmphiboleTEs.svg",fig)

outcrop1 = ["23SD02B-1","23SD02B-5","21SD08"]
outcrop2 = ["23SD03C-5","23SD03C-1","23SD03E"]
outcrop3 = ["22SD55E","22SD55B","22SD55D"]
outcrop4 = ["23SD20B","23SD20C","23SD20F","23SD20G"]
fig = plotsamples(geochemData,[outcrop2,outcrop3,outcrop4,outcrop1],norm,"chondrite")
# save("SpiderPlots/BRC_TEs.svg",fig)
# calcEuAnomaly!(ampData,norm) 
# sumREEs!(ampData)



# sampleIndex = 14
# plotSpider!(filterSampleName(ampData,samples[sampleIndex]),norm,ax,fig,sampleIndex,samples[sampleIndex])
# sampleIndex = 15
# plotSpider!(filterSampleName(ampData,samples[sampleIndex]),norm,ax,fig,sampleIndex,samples[sampleIndex])
# for i in 21:23
#     sampleIndex = i
#     plotSpider!(filterSampleName(ampData,samples[sampleIndex]),norm,ax,fig,sampleIndex,samples[sampleIndex])
# end

# for i in 24:24#lastindex(samples)-1
#     sampleIndex = i
#     plotSpider!(filterSampleName(ampData,samples[sampleIndex]),norm,ax,fig,sampleIndex,samples[sampleIndex])
# end

# for i in 9:11
#     sampleIndex = i
#     plotSpider!(filterSampleName(ampData,samples[sampleIndex]),norm,ax,fig,sampleIndex,samples[sampleIndex])
# end



# fig2 = Figure(size = (1200,900))
# ax2 = Axis(fig2[1,1],yscale = log10,xscale = log10)

# for sam in samples
#     selAmp = filterSampleName(ampData,sam)

#     # scatter!(ax2,filterSampleName(ampData,sam)[!,:TotalREE],filterSampleName(ampData,sam)[!,:Eu_EuA], label = sam)
#     # scatter!(ax2,filterSampleName(ampData,sam)[!,:Rb85_ppm_mean],filterSampleName(ampData,sam)[!,:Eu_EuA], label = sam)
#     # scatter!(ax2,filterSampleName(ampData,sam)[!,:Sr88_ppm_mean],filterSampleName(ampData,sam)[!,:Eu_EuA], label = sam)
#     scatter!(ax2,filterSampleName(ampData,sam)[!,:Y89_ppm_mean],filterSampleName(ampData,sam)[!,:Eu_EuA], label = sam)
# end

# for i in 1:lastindex(samplesWGeochem)
   
#     selAmp = filterSampleName(ampData,samplesWGeochem[i])
#     sampleNorm = toNorm(geochemData,geochemSamples[i])
#     calcEuAnomaly!(selAmp,sampleNorm)
#     normElem!(selAmp,sampleNorm,"Rb")
#     normElem!(selAmp,sampleNorm,"Sr")
#     normElem!(selAmp,sampleNorm,"Y")
#     # scatter!(ax2,filterSampleName(ampData,sam)[!,:TotalREE],filterSampleName(ampData,sam)[!,:Eu_EuA], label = sam)
#     # scatter!(ax2,selAmp[!,:Rb85_ppm_mean],selAmp[!,:Eu_EuA], label = samplesWGeochem[i])
#     # scatter!(ax2,selAmp[!,:Sr88_ppm_mean],selAmp[!,:Eu_EuA], label = samplesWGeochem[i])
#     # scatter!(ax2,selAmp[!,:Y89_ppm_mean],selAmp[!,:Eu_EuA], label = samplesWGeochem[i])
#     scatter!(ax2,selAmp[!,:Rb_norm],selAmp[!,:Eu_EuA], label = samplesWGeochem[i])
#     # scatter!(ax2,selAmp[!,:Sr_norm],selAmp[!,:Eu_EuA], label = samplesWGeochem[i])
#     # scatter!(ax2,selAmp[!,:Y_norm],selAmp[!,:Eu_EuA], label = samplesWGeochem[i])

# end
# fig2[1,2] = Legend(fig2,ax2,framevisible = false,merge = true)
# display(GLMakie.Screen(),fig2)