using GLMakie
using Isoplot
using DataFrames
using CSV
using Statistics
using FileIO
using NativeFileDialog
using JLD2
using SkipNan
using LinearRegression
# using LaTeXStrings
const MAX_SELECT = 20
const FIG_SIZE = (1200,900)
# myColours = ["#FF8CA6","#6D8FFF","#F9C52A","#A9EFA5", "#E6A7FB","#000000","#FF0000","#008000","#0000FF","#4B0082","#FF8C00"]
# myTheme = Theme(
#     fonts = (; regular = "B612"),
#     fontsize = 20,
#     palette = (color =["#FF8CA6","#6D8FFF","#F9C52A","#A9EFA5", "#E6A7FB","#000000","#FF0000","#008000","#0000FF","#4B0082","#FF8C00"], ),
#     Axis = ( 
#         xlabelsize = 28,
#         ylabelsize = 28,
#         xminorticksvisible = true,
#         yminorticksvisible = true,
#         xticksmirrored = true,
#         # yticksmirrored = true,
#         xtickalign = 1,
#         ytickalign = 1,
#         xminortickalign = 1,
#         yminortickalign = 1,
#         xticksize = 10,
#         xtickwidth = 2,
#         yticksize = 10,
#         ytickwidth = 2,
#         xminorticksize = 5,
#         xminortickwidth = 2,
#         yminorticksize = 5,
#         yminortickwidth = 2,
#         xgridvisible = false,
#         ygridvisible = false
#         )
# )
include("PlotDefaults.jl")


mutable struct TransectData
    teData::DataFrame
    UPbData::DataFrame
    UPbAvgs::DataFrame
    teAvgs::DataFrame
    teNorm::DataFrame
    elemList::Array{String}
end

TransectData() = TransectData(DataFrame(),DataFrame(),DataFrame(),DataFrame(),DataFrame(),Array{String}([]))

function isVisible(fig::Figure)
    if length(fig.scene.current_screens) > 0
        return true
    end
    return false
end

function changeElem!(thisData::TransectData,ax::Axis,y::Observable,elem::String)
    y1 = Matrix(thisData.teData[!,Regex(elem*"\\d+")])
    y[] = vec(y1)
    ax.ylabel = elem*" (μg/g)"
    ylims!(ax,0,maximum(y1))
end

function initFig!(lineData::DataFrame, scanSpeed::Int,elem::String)
    
    fig = Figure(size = FIG_SIZE); display(GLMakie.Screen(),fig)
    instanceData = TransectData()
    instanceData.teData = lineData
    elemList = Array{String}([])
    for i=3:lastindex(names(lineData))
        elemName = match(r"[A-z]+",names(instanceData.teData)[i]).match
        push!(elemList,elemName)
    end
    instanceData.elemList = elemList
    
    
    ax = Axis(fig[1,2])
    Makie.deactivate_interaction!(ax, :rectanglezoom)
    Makie.deactivate_interaction!(ax, :scrollzoom)
  
    colsize!(fig.layout,1,Fixed(300))#So dropdown doesnt take too much space
    colsize!(fig.layout,2,Aspect(1,1)) #To enforcce square plot
    rowsize!(fig.layout,1,700)


    instanceData.teData[!,:Distance] = instanceData.teData[!,"Elapsed Time"]*scanSpeed
    x = Observable(Vector(instanceData.teData[!,:Distance]))
    #Have to convert to Matrix then vectorize it
    #Not sure why its not needed for the distance column
    y1 = Matrix(instanceData.teData[!,Regex("$elem\\d+")])
    y = Observable(vec(y1))

    ax.xlabel = "Distance (μm)"
    ax.ylabel = elem*" (μg/g)"
    
    xMax = round(maximum(instanceData.teData[!,:Distance]),sigdigits=2)
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
    xlims!(ax,0,maximum(x[]))
    lines!(ax,x,y,linewidth=2, color = Cycled(1))
    ylims!(ax,0,maximum(y1))

    menu = Menu(fig, options = elemList, default = elem)
    fig[1,1][1,1]=vgrid!(
        Label(fig,"Left Axis:",width=nothing),
        menu,tellheight=false
    )
    on(menu.selection) do s
        changeElem!(instanceData,ax,y,s)
    end

    #This works but Im disabling it for now since its easier to test with it off
    loadUPb = Button(fig[1,1][2,1],label = "Load U-Pb data")
    
    on(loadUPb.clicks) do clicks
        delete!(loadUPb)
        fileName = pick_file(;filterlist="csv")
        instanceData.UPbData = DataFrame(CSV.File(fileName))
        addUPbData!(instanceData,scanSpeed,fig,ax,x)
        
    end

    loadImg = Button(fig[1,1][5,1],label="Load zircon image",tellheight=false)
    selStarts = Observable(Float64[])
    selStops = Observable(Float64[])
    imgFig = Figure()
    # lineSelections = DataSelect(Observable(Float64[]),Observable(Float64[]))
    on(loadImg.clicks) do clicks
        fileName = pick_file(;filterlist="jpg;tif;bmp;png")
        importZrnImage(imgFig,fileName,x[],selStarts,selStops)
    end

    sFig, sAx = initSpiderFig("chondrite")
    isoFig, isoAx = initIsoplotFig(t1 = 800, t2 = 1600)
    initRangeSelect!(ax,fig,instanceData,sAx,sFig,"chondrite.csv",selStarts,selStops,isoplotAx=isoAx, isoplotFig = isoFig)

    # saveState = Button(fig[1,1][6,1],label = "Save state", tellheight=false)
    expData = Button(fig[1,1][6,1],label = "Export data", tellheight = false)

    # on(saveState.clicks) do clicks
    #     fileName = save_dialog("Save the current figures",GtkNullContainer(),(GtkFileFilter("*.jld2")))
    #     jldsave(fileName;fig,sFig,isoFig,imgFig,instanceData)
    # end

    on(expData.clicks) do clicks
        dir = pick_folder()
        if isdir(dir)
            CSV.write(dir*"/TraceElementData.csv",instanceData.teData)
            CSV.write(dir*"/UPbData.csv",instanceData.UPbData)
            CSV.write(dir*"/TraceElementsAvgs.csv",instanceData.teAvgs)
            CSV.write(dir*"/TraceElementsNormalized.csv",instanceData.teNorm)
            CSV.write(dir*"/UPbAvgs.csv",instanceData.UPbAvgs)
            GLMakie.save(dir*"/Transect.png",fig)
            GLMakie.save(dir*"/Spider.png",sFig)
            GLMakie.save(dir*"/Isoplot.png",isoFig)
            GLMakie.save(dir*"/ZrnImage.png",imgFig)
        end
        rescaleUPb
    end
    return fig,ax, x, y
end

function changeAxisVar!(lineData::DataFrame,ax::Axis,y::Observable,varName::String, altData::DataFrame)
    y1 =lineData[!,varName]
    y[] = Vector(y1)
    ax.ylabel = varName
    if occursin("rescaled",varName)
        element = match(r"[^_]",varName).match
        yalt = altData[!,Regex(element)]
        ylims!(ax,0,maximum(yalt[!,1]))
    else
        ylims!(ax,0,mean(skipnan(y1))*2)
    end
end

function adjustXVar(x::Observable,xAdj::Float64)
    x[] = x[] .+ xAdj
end

function addUPbData!(thisData::TransectData,scanSpeed::Int,fig::Figure,teAx::Axis,teX::Observable)


    thisData.UPbData[!,:Distance] = thisData.UPbData[!,"Elapsed Time"]*scanSpeed
    yvar = "Approx_U_PPM"
    ax2 = Axis(fig[1,2],
        yaxisposition = :right)
    hidespines!(ax2,:l,:t,:b)
    hidespines!(teAx,:r)
    hidexdecorations!(ax2)

    
    Makie.deactivate_interaction!(ax2, :rectanglezoom)
    Makie.deactivate_interaction!(ax2, :scrollzoom)
    # hidespines!(teAx,:r)
    xBase = Vector(thisData.UPbData[!,:Distance]) #Used since x gets changed relative to this
    x = Observable(xBase)
    y = Observable(Vector(thisData.UPbData[!,yvar]))

    ax2.ylabel = yvar
    xlims!(ax2,0,maximum(teX[]))

    #Menu set up
    menu = Menu(fig, options = names(thisData.UPbData)[3:end], default = yvar)
    fig[1,1][2,1]=vgrid!(
        Label(fig,"Right axis:",width=nothing),
        menu,tellheight=false
    )
    on(menu.selection) do s
        changeAxisVar!(thisData.UPbData,ax2,y,s,thisData.teData)
    end

    #Slider set up
    xSlider = SliderGrid(fig,
    (range = -20:0.01:20, format = "{:.1f}μm",startvalue = 0,snap = false),
    (range = -20:0.01:20, format = "{:.1f}%",startvalue = 0,snap = false),
    width = 300)

    #This has to be before the plot command
    x = lift(xSlider.sliders[1].value,xSlider.sliders[2].value) do xAdj,xStr
        stretch = 1 + xStr/100
        thisData.UPbData[!,:DistanceMod] = xBase.*stretch .+ xAdj
        xBase.*stretch .+ xAdj 
        # 
    end
    # thisData.UPbData[!,:DistanceMod] = x[]
    fig[1,1][3,1] = vgrid!(
        Label(fig,"x adjustment:",width=nothing),
        xSlider)
    
    lines!(ax2,x,y,linewidth = 2, color = Cycled(2))

    rescaleUPb = Button(fig[1,1][4,1],label = "Rescale U and Th", tellheight=false)
    on(rescaleUPb.clicks) do clicks

        rescale_UTh!(thisData.teData, thisData.UPbData)
        menu = Menu(fig, options = names(thisData.UPbData)[3:end], default = yvar)

        fig[1,1][2,1]=vgrid!(
            Label(fig,"Right axis:",width=nothing),
            menu,tellheight=false
        )
        on(menu.selection) do s
            changeAxisVar!(thisData.UPbData,ax2,y,s,thisData.teData)
        end
    end
    
end

function NaNIfNotPositive(x::Number)
    if x <= 0
        return NaN
    end

    return x
end

function averageLineTE(thisData::TransectData;start::AbstractFloat=0.0, stop::AbstractFloat=-1.0)

    if stop < 0
        stop = thisData.teData[end,:Distance]
    end

    averageElems = DataFrame(Element=[],Concentration=[],Error=[])
    teRow = [start,stop]
    selection = filter([:Distance]=> x -> start<=x<=stop,thisData.teData)
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

function averageUPb(thisData::TransectData;start::AbstractFloat=0.0,stop::AbstractFloat = -1.0)

    if stop < 0
        stop = thisData.UPbData[end,:Distance]
    end

    
    selection = filter([:DistanceMod]=> x -> start<=x<=stop,thisData.UPbData)
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



function initRangeSelect!(ax::Axis,fig::Figure,thisData::TransectData,spiderAx::Axis,spiderFig::Figure,normFile::String,selStarts::Observable,selStops::Observable;isoplotAx::Axis,isoplotFig::Figure)
    
    norm = DataFrame(CSV.File(normFile))
    ax3 = Axis(fig[1,2],
        xaxisposition = :top)
    hidespines!(ax3,:l,:r,:b)
    hidexdecorations!(ax3)
    hideydecorations!(ax3)
    Makie.deactivate_interaction!(ax3, :rectanglezoom)
    Makie.deactivate_interaction!(ax3, :scrollzoom)
    xlims!(ax3,0,maximum(thisData.teData[!,:Distance]))
    ylims!(ax3,0,1)
    selectionsMid = Float64[]

    thisData.UPbAvgs = DataFrame(Selection = Array{String}([]), x1 = Array{Float64}([]),x2=Array{Float64}([]),Analyses = Array{UPbAnalysis{Float64}}([]),
                                Pb207U235=Array{Float64}([]), Pb207U235_se1=Array{Float64}([]), Pb206U238=Array{Float64}([]), Pb206U238_se1=Array{Float64}([]),
                                Pb207Pb206=Array{Float64}([]),Pb207Pb206_se1=Array{Float64}([]),cor_0735_0638=Array{Float64}([]),cor_0706_0638=Array{Float64}([]))

    thisData.teAvgs = DataFrame(Selection = Array{String}([]), x1 = Array{Float64}([]),x2=Array{Float64}([]))
    for elem in thisData.elemList
        thisData.teAvgs[!,Symbol(elem*"_mean")] = Array{Float64}([])
        thisData.teAvgs[!,Symbol(elem*"_1std")] = Array{Float64}([])
    end
   
    thisData.teNorm = DataFrame(Selection = Array{String}([]), x1 = Array{Float64}([]),x2=Array{Float64}([]))
    for elem in norm[!,:Element]
        thisData.teNorm[!,elem] = Array{Float64}([])
    end

    # selectionsStart = Float64[]
    # selectionsStop = Float64[]
    #Implement a third axis on ax with everything hidden but axis labels on the top
    #Then add "S1, S2, S3..." to the axis label as was done with elements on REE plot
    sRange = select_rectangle(ax.scene)
    @show sRange
    on(sRange) do rect
       
        x1 = rect.origin[1]
        x2 = rect.widths[1] + x1
        xmid = mean([x1,x2])

        push!(selectionsMid,xmid)
        
        # push!(selectionsStart,x1)
        # push!(selectionsStop,x2)
        push!(selStarts[],x1)
        push!(selStops[],x2)
        sIndex = length(selectionsMid)
        # selections.start[] = selectionsStart
        # selections.stop[] = selectionsStop
        lowLimitX = Observable([x1,x1])
        upLimitX = Observable([x2,x2])
        yLimits = Observable([0,1])
        
        lines!(ax3,lowLimitX,yLimits,linestyle=:solid,color = :gray, linewidth = 3)
        lines!(ax3,upLimitX,yLimits,linestyle=:solid,color = :gray, linewidth = 3)
        ax3.xticklabelsvisible = true
        ax3.xticks = (selectionsMid,["S"*string(i) for i = 1:sIndex])
        # text!()
        elemAvg, teRow = averageLineTE(thisData,start = x1,stop=x2)
     
        push!(thisData.teAvgs,["S"*string(sIndex);teRow])
        elemRatio = plotSpider!(elemAvg,norm,spiderAx,spiderFig,sIndex)
        push!(thisData.teNorm,["S"*string(sIndex);x1;x2;elemRatio])
        if nrow(thisData.UPbData) > 0
           
            sel = "S"*string(sIndex)
            selAnal, selRaw = averageUPb(thisData,start=x1,stop=x2)
           
            plotUPb!(isoplotFig,isoplotAx,selAnal,sIndex)
            push!(thisData.UPbAvgs,[sel;x1;x2;selAnal;selRaw])
            xmin, xmax, ymin, ymax = datalimits(thisData.UPbAvgs[!,:Analyses])
            limits!(isoplotAx,xmin,xmax,ymin,ymax)
        end
       
    end
        
        

    
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
    # push!(thisData.teNorm,elemRatio)

    lines!(ax,x,NaNIfNotPositive.(elemRatio),linewidth = 4, color = Cycled(selectionIndex), label = "S"*string(selectionIndex))
    
    # lines!(ax,x,NaNIfNotPositive.(elemRatio + elemError),linewidth = 2,linestyle = :dash, color = Cycled(selectionIndex))

    # lines!(ax,x,NaNIfNotPositive.(elemRatio - elemError),linewidth = 2,linestyle = :dash, color = Cycled(selectionIndex))
    ax.xticks = (1:lastindex(elemList),elemList)
    [delete!(leg) for leg in fig.content if leg isa Legend]
    fig[1,2] = Legend(fig,ax,framevisible = false)
    if !isVisible(fig)
        
        display(GLMakie.Screen(),fig)
    end

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

function plotUPb!(fig::Figure,ax::Axis,analysis::UPbAnalysis,selectionIndex::Integer)

    
    plot!(ax,analysis, color = Cycled(selectionIndex),alpha = 0.6)
    
    text!(ax,analysis.μ[1],analysis.μ[2],text="S"*string(selectionIndex),fontsize = 15)

    if !isVisible(fig)
        
        display(GLMakie.Screen(),fig)
    end
end

function lineButtonPressActions!(ax::Axis,linePoint::Observable)
    activateScreen = true    
    mousePos = events(ax).mouseposition
   
    crosshair = lift(mousePos) do pos
        mouseposition(ax)
    end
    Makie.deactivate_interaction!(ax, :rectanglezoom)
    Makie.deactivate_interaction!(ax, :scrollzoom)
    
    chCursor1 = scatter!(ax,crosshair, marker = '+', color = :white, markersize = 50)
    chCursor2 = scatter!(ax,crosshair, marker = '+', color = :black, markersize = 40)
    

    on(events(ax).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.release && activateScreen
             
            #For some reason this doesnt work if I use mouseposition as a point and plot that
            #As I did for the cursor
            linePoint[] = Point2f(Float32(mouseposition(ax)[1]),Float32(mouseposition(ax)[2]))
           
            
            delete!(ax,chCursor1)
            delete!(ax,chCursor2)
            Makie.activate_interaction!(ax, :rectanglezoom)
            Makie.activate_interaction!(ax, :scrollzoom)
            activateScreen = false
            
        end
    end
end

function importZrnImage(fig::Figure,fileName::String,dataX::Array{Float64},selStarts::Observable,selStops::Observable)

    zrnImg = load(assetpath(fileName))

    # fig = Figure()
    imgPlot = image(fig[1,1],zrnImg,axis = (;aspect = DataAspect()))
    hidexdecorations!(imgPlot.axis)
    hideydecorations!(imgPlot.axis)
    hidespines!(imgPlot.axis,:t,:b,:l,:r)
    xlims!(imgPlot.axis,0,size(zrnImg)[1])
    ylims!(imgPlot.axis,0,size(zrnImg)[2])
    selectStart = Button(fig[1,2][1,1],label = "Reference Transect Start",tellheight=false)
    selectEnd = Button(fig[1,2][2,1],label = "Reference Transect End",tellheight=false)
    # pointSelected = false
    # lineStartX = Observable(-50.0)
    # lineStartY = Observable(-50.0)
    
    lineStart = Observable(Point2f(-50.0,-50.0))
    scatter!(imgPlot.axis,lineStart, marker = '+', color = :white, markersize = 50)
    scatter!(imgPlot.axis,lineStart, marker = '+', color = :black, markersize = 40)
    on(selectStart.clicks) do clicks
        lineButtonPressActions!(imgPlot.axis,lineStart)
    end

    # lineEndX = Observable(-50.0)
    # lineEndY = Observable(-50.0)
    lineEnd = Observable(Point2f(-50.0,-50.0))
    scatter!(imgPlot.axis,lineEnd, marker = '+', color = :white, markersize = 50)
    scatter!(imgPlot.axis,lineEnd, marker = '+', color = :black, markersize = 40)
    on(selectEnd.clicks) do clicks
        lineButtonPressActions!(imgPlot.axis,lineEnd)
    end

    
    selectPoints = [Observable(Point2f[]) for i in 1:MAX_SELECT]
    for i in 1:lastindex(selectPoints)
        #Try again without selections as a struct but just seperate observables
        
        selectPoints[i] = lift(lineStart,lineEnd,selStarts, selStops) do p1,p2,starts,stops
            scaleRatio = sqrt((p2[1]-p1[1])^2+(p2[2]-p1[2])^2)/(dataX[end]-dataX[1])
            scaledD = (dataX[end]-dataX[1])*scaleRatio
            dxRatio = scaledD/(p2[1]-p1[1])
            dyRatio = scaledD/(p2[2]-p1[2])
            xstart = 0.0
            xstop = 0.0
            ystart = 0.0
            ystop = 0.0
            
            if i <= lastindex(starts)
               
                xstart = p1[1] + starts[i]*scaleRatio/dxRatio
                xstop = p1[1] + stops[i]*scaleRatio/dxRatio
                ystart = p1[2] + starts[i]*scaleRatio/dyRatio
                ystop = p1[2] + stops[i]*scaleRatio/dyRatio
            end
            [Point2f(xstart,ystart),Point2f(xstop,ystop)]
        end


    end
    for i in 1:lastindex(selectPoints)
        lines!(imgPlot.axis,selectPoints[i],linewidth = 20)

        textX = lift(selectPoints[i]) do sel
           mid = midpoint(sel[1],sel[2])
           mid[1]
        end
        
        textY = lift(selectPoints[i]) do sel
            mid = midpoint(sel[1],sel[2])
            mid[2]
         end
         
        text!(imgPlot.axis,textX,textY,text="S"*string(i), align = (:center,:center))
    end


    display(GLMakie.Screen(),fig)

end

function midpoint(p1::Point2f, p2::Point2f)
    
    x1, x2 = p2[1] > p1[1] ? (p1[1],p2[1]) : (p2[1],p1[1])
    y1, y2 = p2[2] > p1[2] ? (p1[2],p2[2]) : (p2[2],p1[2])

    return Point2f((x2-x1)/2+x1,(y2-y1)/2+y1)

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

    regression = linregress(upbData[!,:Approx_U_PPM],U_TE_match,weights)
    slope = LinearRegression.slope(regression)
    y_intercept = LinearRegression.bias(regression)
    
    upbData[!,:U_rescaled] = upbData[!,:Approx_U_PPM].*slope .+ y_intercept
    println("R^2 U = " *string(r_squared(U_TE_match,upbData[!,:U_rescaled])))
    upbData[!,:Th_rescaled] = upbData[!,:U_rescaled]./upbData[!,"Final U/Th"]
    println("R^2 Th = " *string(r_squared(Th_TE_match,upbData[!,:Th_rescaled])))

end

GLMakie.activate!()
set_theme!(myTheme)
# zrn37 = DataFrame(CSV.File(raw"C:\Users\Sabas\OneDrive - University of Waterloo\Documents\Waterloo\LAICPMS\2023_06_22\20SD06-2_Zrn37.csv"))
# zrn41 = DataFrame(CSV.File(raw"C:\Users\Sabas\OneDrive - University of Waterloo\Documents\Waterloo\LAICPMS\2023_06_22\20SD06-2_Zrn41.csv"))
# zrn10 = DataFrame(CSV.File(raw"C:\Users\Sabas\OneDrive - University of Waterloo\Documents\Waterloo\LAICPMS\2023_06_22\20SD06-2_Zrn10.csv"))
# zrn37UPb = DataFrame(CSV.File(raw"C:\Users\Sabas\OneDrive - University of Waterloo\Documents\Waterloo\LAICPMS\2023_06_06\20SD06_zrn37.csv"))
fileName = pick_file(;filterlist="csv")

zrnTE = DataFrame(CSV.File(fileName))

sIndex = 1
fig,ax,x,conc = initFig!(zrnTE,2,"U")
# ax2,x2,y2 = addUPbData!(zrn37UPb,2,fig,ax,x)
# sFig, sAx = initSpiderFig("chondrite")
# initRangeSelect!(ax,fig,zrn37,sAx,sFig,"Chondrite.csv")

# clFig,clAx = importCLImage(raw"C:\Users\Sabas\OneDrive - University of Waterloo\Documents\Waterloo\SEM\Zircon CL\SVG Files\20SD06\Zrn37.png")
