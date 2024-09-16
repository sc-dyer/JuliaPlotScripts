using GLMakie
#using AlgebraOfGraphics
using DataFrames
using CSV
using Measurements
using FileIO
using NativeFileDialog
using Statistics

const FIG_SIZE = (1200,900)
const MAX_SELECT = 20

include("PlotDefaults.jl")

function initFig!(hfData, sam_name, line_num; reverse_profile = false)
    
    fig = Figure(size = FIG_SIZE); display(GLMakie.Screen(),fig)
   
    probeLine = selectLine(hfData,sam_name,line_num,reverseProf = reverse_profile)
    
    
    ax = Axis(fig[1,2])
    Makie.deactivate_interaction!(ax, :rectanglezoom)
    Makie.deactivate_interaction!(ax, :scrollzoom)
  
    colsize!(fig.layout,1,Fixed(300))#So dropdown doesnt take too much space
    colsize!(fig.layout,2,Aspect(1,1)) #To enforcce square plot
    rowsize!(fig.layout,1,700)

    ax.xlabel = "Distance (μm)"
    ax.ylabel = "Hf (wt%)"
    
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

    x = Observable(probeLine[!,:Distance])
    #Have to convert to Matrix then vectorize it
    #Not sure why its not needed for the distance column
    
    y = Observable(probeLine[!,:Hf])

    y_err = Observable(probeLine[!,:AbsError])

    errorbars!(x,y, y_err, whiskerwidth = 6,color = "black")
    scatter!(x,y)


    ax.limits = (nothing, nothing, 0.75,1.5)

   

    #This works but Im disabling it for now since its easier to test with it off
   
    loadImg = Button(fig[1,1][1,1],label="Load zircon image",tellheight=false)
    selStarts = Observable(Float64[])
    selStops = Observable(Float64[])
    imgFig = Figure()
    # lineSelections = DataSelect(Observable(Float64[]),Observable(Float64[]))
    on(loadImg.clicks) do clicks
        fileName = pick_file(;filterlist="jpg;tif;bmp;png")
        importZrnImage(imgFig,fileName,x[],selStarts,selStops)
    end

    initRangeSelect!(ax,fig,hfData,selStarts,selStops)

    # saveState = Button(fig[1,1][6,1],label = "Save state", tellheight=false)
    expData = Button(fig[1,1][2,1],label = "Export data", tellheight = false)

    # on(saveState.clicks) do clicks
    #     fileName = save_dialog("Save the current figures",GtkNullContainer(),(GtkFileFilter("*.jld2")))
    #     jldsave(fileName;fig,sFig,isoFig,imgFig,instanceData)
    # end

    on(expData.clicks) do clicks
        dir = pick_folder()
        if isdir(dir)
            selections = DataFrame(selection = String[], x1 = Float64[], x2 = Float64[])
            for i in 1:lastindex(selStarts[])
                push!(selections,["S"*string(i),selStarts[][i],selStops[][i]])
            end
            CSV.write(dir*"/Selections.csv",selections)
            GLMakie.save(dir*"/Transect.png",fig)
            GLMakie.save(dir*"/ZrnImage.png",imgFig)
        end
        
    end
    return fig,ax, x, y
end

function initRangeSelect!(ax,fig,thisData,selStarts,selStops)
    
    ax3 = Axis(fig[1,2],
        xaxisposition = :top)
    hidespines!(ax3,:l,:r,:b)
    hidexdecorations!(ax3)
    hideydecorations!(ax3)
    Makie.deactivate_interaction!(ax3, :rectanglezoom)
    Makie.deactivate_interaction!(ax3, :scrollzoom)
    xlims!(ax3,ax.xaxis.attributes.limits[])
    ylims!(ax3,0,1)
    selectionsMid = Float64[]


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

function selectLine(probeData, samName, lineNum;reverseProf = false)
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


GLMakie.activate!()
set_theme!(myTheme)

probeData = DataFrame(CSV.File("../../Probe/ZrnHf_Modified.csv"))




# ax.xlabel = "Distance (μm)"
# ax.ylabel = "Hf (wt%)"
# ax.yticks = 0.75:0.05:1.5
# ax.limits = (nothing, nothing, 0.75,1.5)



# thisLine = selectLine(probeData,"20SD06-2_Zrn37",1)
# plotProfile(thisLine,ax,fig,"20SD06-2_Zrn37")
# fig
# thisLine = selectLine(probeData,"20SD06-2_Zrn27",1)
# plotProfile(thisLine,ax,fig,"20SD06-2_Zrn27")
# fig
# thisLine = selectLine(probeData,"20SD06-2_Zrn41",1)
# plotProfile(thisLine,ax,fig,"20SD06-2_Zrn41_L1")

# thisLine = selectLine(probeData,"20SD06-2_Zrn41",2)
# plotProfile(thisLine,ax,fig,"20SD06-2_Zrn41_L2")

# thisLine = selectLine(probeData,"20SD06-2_Zrn43",1)
# plotProfile(thisLine,ax,fig,"20SD06-2_Zrn43")

# thisLine = selectLine(probeData,"20SD06-2_Zrn26",1,reverseProf = true)
# plotProfile(thisLine,ax,fig,"20SD06-2_Zrn26")

# thisLine = selectLine(probeData,"20SD17A_Zrn4",1)
# plotProfile(thisLine,ax,fig,"20SD17A_Zrn4_L1")

# thisLine = selectLine(probeData,"20SD17A_Zrn4",2)
# plotProfile(thisLine,ax,fig,"20SD17A_Zrn4_L2")

# thisLine = selectLine(probeData,"20SD17A_Zrn15",1)
# plotProfile(thisLine,ax,fig,"20SD17A_Zrn15_L1")

# thisLine = selectLine(probeData,"20SD17A_Zrn15",2)
# plotProfile(thisLine,ax,fig,"20SD17A_Zrn15_L2")

# thisLine = selectLine(probeData,"20SD17A_Zrn25",1)
# plotProfile(thisLine,ax,fig,"20SD17A_Zrn25_L1")

# thisLine = selectLine(probeData,"20SD17A_Zrn25",2,reverseProf = true)
# plotProfile(thisLine,ax,fig,"20SD17A_Zrn25_L2")

# thisLine = selectLine(probeData,"20SD17A_Zrn28",1,reverseProf = true)
# plotProfile(thisLine,ax,fig,"20SD17A_Zrn28_L1")

# thisLine = selectLine(probeData,"20SD17A_Zrn28",2,reverseProf = true)
# plotProfile(thisLine,ax,fig,"20SD17A_Zrn28_L2")

# thisLine = selectLine(probeData,"21SD68_Zrn9",1,reverseProf = true)
# plotProfile(thisLine,ax,fig,"21SD68_Zrn9")

# thisLine = selectLine(probeData,"21SD68_Zrn27",1)
# plotProfile(thisLine,ax,fig,"21SD68_Zrn27_L1")

# thisLine = selectLine(probeData,"21SD68_Zrn27",2)
# plotProfile(thisLine,ax,fig,"21SD68_Zrn27_L2")

# thisLine = selectLine(probeData,"21SD68_Zrn7",1,reverseProf = true)
# plotProfile(thisLine,ax,fig,"21SD68_Zrn7_L1")

# thisLine = selectLine(probeData,"21SD68_Zrn7",2,reverseProf = true)
# plotProfile(thisLine,ax,fig,"21SD68_Zrn7_L2")

# initFig!(probeData,"20SD06-2_Zrn37",1)
# initFig!(probeData,"20SD06-2_Zrn26",1,reverse_profile = true)
# initFig!(probeData,"20SD06-2_Zrn41",1)
initFig!(probeData,"20SD06-2_Zrn43",1)
