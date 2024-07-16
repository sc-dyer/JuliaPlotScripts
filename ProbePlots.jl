using CairoMakie
#using AlgebraOfGraphics
using DataFrames
using CSV
using Measurements

function selectLine(probeData::DataFrame,samName::String,lineNum::Integer;reverseProf::Bool = false)
    newDF = filter([:Sample,:Line] => (Sample,Line)->Sample ==samName&&Line==lineNum,probeData)
    newDF[:,:AbsError] = (newDF[!,:PercentError]/100).*newDF[!,:Hf]
    if newDF[1,:Distance] > 0
        distAdjust = newDF[1,:Distance]
        newDF[!,:Distance] = newDF[!,:Distance] .- distAdjust
    end
    if reverseProf
        distVector:: Array{Float64}=[]
        newDF = reverse(newDF)
        lastDist = newDF[1,:Distance]
        newLast = 0
        for row in eachrow(newDF)
            currDist = row[:Distance]
            newDist = lastDist - currDist + newLast
            push!(distVector,newDist)
            newLast = newDist
            lastDist = currDist
        end
        newDF[!,:Distance] = distVector
        #newDF[!,:Distance] = reverse(newDF[!,:Distance])
    end
    return newDF
end

function plotProfile(probeLine::DataFrame, ax::Axis,fig::Figure,fileName::String)
    empty!(ax)
    xMax = round(maximum(probeLine[!,:Distance]),sigdigits=2)
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
    errorbars!(probeLine[!,:Distance],probeLine[!,:Hf], probeLine[!,:AbsError], whiskerwidth = 6)
    scatter!(probeLine[!,:Distance],probeLine[!,:Hf],color = "#FF8CA6")
    save("ProbePlots//"*fileName*".svg",fig)
    #
end
probeData = DataFrame(CSV.File("../../Probe/ZrnHf_Modified.csv"))



fig = Figure(font = "Trebuchet",fontsize = 18)
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

ax.xlabel = "Distance (μm)"
ax.ylabel = "Hf (wt%)"
ax.yticks = 0.75:0.05:1.5
ax.limits = (nothing, nothing, 0.75,1.5)



thisLine = selectLine(probeData,"20SD06-2_Zrn37",1)
plotProfile(thisLine,ax,fig,"20SD06-2_Zrn37")
fig
thisLine = selectLine(probeData,"20SD06-2_Zrn27",1)
plotProfile(thisLine,ax,fig,"20SD06-2_Zrn27")
fig
thisLine = selectLine(probeData,"20SD06-2_Zrn41",1)
plotProfile(thisLine,ax,fig,"20SD06-2_Zrn41_L1")

thisLine = selectLine(probeData,"20SD06-2_Zrn41",2)
plotProfile(thisLine,ax,fig,"20SD06-2_Zrn41_L2")

thisLine = selectLine(probeData,"20SD06-2_Zrn43",1)
plotProfile(thisLine,ax,fig,"20SD06-2_Zrn43")

thisLine = selectLine(probeData,"20SD06-2_Zrn26",1,reverseProf = true)
plotProfile(thisLine,ax,fig,"20SD06-2_Zrn26")

thisLine = selectLine(probeData,"20SD17A_Zrn4",1)
plotProfile(thisLine,ax,fig,"20SD17A_Zrn4_L1")

thisLine = selectLine(probeData,"20SD17A_Zrn4",2)
plotProfile(thisLine,ax,fig,"20SD17A_Zrn4_L2")

thisLine = selectLine(probeData,"20SD17A_Zrn15",1)
plotProfile(thisLine,ax,fig,"20SD17A_Zrn15_L1")

thisLine = selectLine(probeData,"20SD17A_Zrn15",2)
plotProfile(thisLine,ax,fig,"20SD17A_Zrn15_L2")

thisLine = selectLine(probeData,"20SD17A_Zrn25",1)
plotProfile(thisLine,ax,fig,"20SD17A_Zrn25_L1")

thisLine = selectLine(probeData,"20SD17A_Zrn25",2,reverseProf = true)
plotProfile(thisLine,ax,fig,"20SD17A_Zrn25_L2")

thisLine = selectLine(probeData,"20SD17A_Zrn28",1,reverseProf = true)
plotProfile(thisLine,ax,fig,"20SD17A_Zrn28_L1")

thisLine = selectLine(probeData,"20SD17A_Zrn28",2,reverseProf = true)
plotProfile(thisLine,ax,fig,"20SD17A_Zrn28_L2")

thisLine = selectLine(probeData,"21SD68_Zrn9",1,reverseProf = true)
plotProfile(thisLine,ax,fig,"21SD68_Zrn9")

thisLine = selectLine(probeData,"21SD68_Zrn27",1)
plotProfile(thisLine,ax,fig,"21SD68_Zrn27_L1")

thisLine = selectLine(probeData,"21SD68_Zrn27",2)
plotProfile(thisLine,ax,fig,"21SD68_Zrn27_L2")

thisLine = selectLine(probeData,"21SD68_Zrn7",1,reverseProf = true)
plotProfile(thisLine,ax,fig,"21SD68_Zrn7_L1")

thisLine = selectLine(probeData,"21SD68_Zrn7",2,reverseProf = true)
plotProfile(thisLine,ax,fig,"21SD68_Zrn7_L2")

println("Done")
