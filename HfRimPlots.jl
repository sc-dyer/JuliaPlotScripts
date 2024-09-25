# using GLMakie
using CairoMakie
using DataFrames
using CSV
using Loess
using NativeFileDialog

const FIG_SIZE = (1200,900)
const SPECIAL_SCALE = 40
include("PlotDefaults.jl")

mutable struct RimData
    name::String
    hfdf::DataFrame
    x::Observable
    y::Observable
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

function importfiles(filedir, probedata; special = false)
    selections = DataFrame(CSV.File(joinpath(filedir,"Selections.csv")))

    # hf_df = DataFrame(CSV.File(joinpath(filedir,"TraceElementData.csv")))
    
    
    rim_hf = DataFrame[]
    for row in eachrow(selections)
        
        probedf = selectLine(probedata,row[:Sample],row[:Line],reverseProf=row[:rev])
        
        if special
            push!(rim_hf,select_range(probedf,row[:x1],row[:x2], x0 = row[:x0], special = special))
        else
            push!(rim_hf,select_range(probedf,row[:x1],row[:x2], x0 = row[:x0]))
        end
    end

    return rim_hf
end

function select_range(df,x1,x2;x0 = -1, special = false)
    rimdf = df
    x_start = x1
    if x0 >= 0
        x_start = x0
    end

    if x_start < x2
        rimdf = filter(:Distance => x -> x_start <= x <= x2, df)
    elseif x_start > x2
        rimdf = filter(:Distance => x -> x2 <= x <= x_start, df)

        # dmin = minimum(rimdf[!,:Distance])
        dmax = maximum(rimdf[!,:Distance])


        old_d = rimdf[!,:Distance]
        cumulative = dmax
        d_array = Float64[cumulative]
        for i in 2:lastindex(old_d)

            d = old_d[i] - old_d[i-1]
            cumulative -= d
            push!(d_array,cumulative)
        end

        # reverse!(rimdf)
        rimdf[!,:Distance] = d_array
    end

    #set first to x = 0
    xmin = minimum(rimdf[!,:Distance])
    if x0 >= 0
        if x1 < x2
            xmin = x1
        else
            #distance was reversed so x1 is no longer where the rim start is
            
            xmin = (x0-x1) + x2
        end
    end

    rimdf[!,:Distance] = rimdf[!,:Distance] .- xmin
    
    if special
        arbitrary_scale!(rimdf, x2 = SPECIAL_SCALE)
    else
        arbitrary_scale!(rimdf)
    end
    return rimdf
end

function arbitrary_scale!(df;x1 = 0.0, x2 = 100.0)
    ratio = (x2 - x1)/maximum(df[!,:Distance])
    # x0 = x1
    # if df[1,:Distance] < 0
    #     i = 1
    #     while df[i,:Distance] < 0
    #         i += 1
    #     end
        
    #     increment = (x2 - x1)/(nrow(df)-i)
    #     x0 = 0.0 - (i-1)*increment
        
    # end
    
    x_arb = df[!,:Distance] .* ratio
    # x_arb = [collect(x0:increment:x1);collect(x1+increment:increment:x2)]
   
    df[!,:x_arb] = x_arb

end

function import_all(directory, probefile, is_scaled; special = nothing)
    entries = readdir(directory)
    probedf = DataFrame(CSV.File(probefile))
    zrns = RimData[]
    xmax = 0
    for entry in entries
        if isdir(joinpath(directory,entry))
            name = entry
            hf = DataFrame()
            if name == special 
                hf = importfiles(joinpath(directory,entry),probedf, special = true)
            else
                hf = importfiles(joinpath(directory,entry),probedf)
            end

            for i in 1:lastindex(hf)
                rimName = name*"_R"*string(i)
                y = hf[i][!,:Hf]

                x = hf[i][!,:Distance]
                if is_scaled
                    x = hf[i][!,:x_arb]
                end
                if maximum(x) > xmax
                    xmax = maximum(x)
                end
              
                zrn = RimData(rimName,hf[i],Observable(x),Observable(y))
                push!(zrns,zrn)
            end
            
        end
    end
    
    return zrns, xmax
end


function smooth!(rim,param)
    x = rim.x[]
    y = rim.y[]
   
    index = lastindex(x)
    nindex = index -1
    while x[nindex] == x[index]
        index = nindex
        nindex -= 1
    end
    
    x = convert.(Float64,x[1:index])
    y = y[1:index]
    
    s = param
    # if length(x) < 100
    #     s = 0.2
    # elseif length(x) < 50
    #     s = 0.8
    # end

    model = loess(x,y,span=s)

    new_y = predict(model,rim.x[])

    rim.y[] = new_y


end


function init_fig(directory,probefile;is_scaled = true, special = nothing)

    fig = Figure(size = FIG_SIZE); display(GLMakie.Screen(),fig)

    ax = Axis(fig[1,2])
    lines!(ax, [0,0], [-1000,1000],linestyle = :dash, color = :black, linewidth=3)
    zrns, xmax = import_all(directory, probefile,is_scaled, special = special)


    colsize!(fig.layout,1,Fixed(300))#So dropdown doesnt take too much space
    colsize!(fig.layout,2,Aspect(1,1)) #To enforcce square plot
    rowsize!(fig.layout,1,700)
    ax.xlabel = "x (μm)"
    if is_scaled
        ax.xlabel = "x (arbitrary scale)"
    end
    ax.ylabel = "Hf (wt%)"

    for zrn in zrns
        # scatter!(ax,zrn.x,zrn.y, label = zrn.name)
        scatterlines!(ax,zrn.x,zrn.y,linewidth=3, label = zrn.name)
    end

    # xMax = round(xmax,sigdigits=2)

    # xInterval = 20
    # if xMax <120
    #     xInterval = 10
    #     if xMax <= 60
    #         xInterval = 5
    #         if xMax <= 25
    #             xInterval = 2
    #         end
    #     end
    # end    
    # ax.xticks = 0:xInterval:xMax
    # xlims!(ax,0,xmax)

    maxy = 0
    minx = Inf
    maxx = -Inf
    for rim in zrns

        if maximum(rim.y[])> maxy
            maxy = maximum(rim.y[])
        end

        if minimum(rim.x[]) < minx
            minx = minimum(rim.x[])
        end

        if maximum(rim.x[]) > maxx
            maxx = maximum(rim.x[])
        end

    end

    ylims!(ax,0,1.1*maxy)
    xlims!(ax,minx,maxx)

    smoothbutton = Button(fig[1,1][2,1],label = "Smooth Data")
    
    on(smoothbutton.clicks) do clicks
        for rim in zrns
            smooth!(rim,0.6)
        end
        
    end


    fig[1,3] = Legend(fig,ax)

    expData = Button(fig[1,1][3,1],label = "Export data", tellheight = false)

    # on(saveState.clicks) do clicks
    #     fileName = save_dialog("Save the current figures",GtkNullContainer(),(GtkFileFilter("*.jld2")))
    #     jldsave(fileName;fig,sFig,isoFig,imgFig,instanceData)
    # end

    on(expData.clicks) do clicks
        
        if isdir(directory)
            GLMakie.save(directory*"/HfRims.png",fig)
        end
        
    end
    return fig
end

function save_fig(directory,probefile,savedir;is_scaled = true, special = nothing)

    fig = Figure(size = FIG_SIZE)

    ax = Axis(fig[1,1])
    lines!(ax, [0,0], [-1000,1000],linestyle = :dash, color = :black, linewidth=3)
    zrns, xmax = import_all(directory, probefile,is_scaled, special = special)

    

    ax.xlabel = "x (μm)"
    if is_scaled
        ax.xlabel = "x (arbitrary scale)"
    end
    ax.ylabel = "Hf (wt%)"

    for zrn in zrns
        # scatter!(ax,zrn.x,zrn.y, label = zrn.name)
        scatterlines!(ax,zrn.x,zrn.y,linewidth=3, label = zrn.name)
    end


    maxy = -Inf
    minx = Inf
    maxx = -Inf
    miny = Inf
    for rim in zrns

        if maximum(rim.y[])> maxy
            maxy = maximum(rim.y[])
        end

        if minimum(rim.y[])< miny
            miny = minimum(rim.y[])
        end

        if minimum(rim.x[]) < minx
            minx = minimum(rim.x[])
        end

        if maximum(rim.x[]) > maxx
            maxx = maximum(rim.x[])
        end

    end
    ymax = round(1.05*maxy,digits = 1)
    ymin = round(miny-0.05*maxy,digits = 1)
    @show ymax, ymin
    ylims!(ax,ymin,ymax)
    # xlims!(ax,minx,maxx)
    ax.xticklabelsvisible = false
    
    ax.yticks = ymin:0.1:ymax
    fig[1,2] = Legend(fig,ax)
    save(joinpath(savedir,"HfRimPlots.svg"),fig)
end

# GLMakie.activate!()
set_theme!(myTheme)
# zrns = init_fig("HfPlots/","../../Probe/ZrnHf.csv",is_scaled=true, special = "Zrn26")

save_fig("HfPlots/","../../Probe/ZrnHf.csv","HfPlots/",is_scaled=true, special = "Zrn26")
