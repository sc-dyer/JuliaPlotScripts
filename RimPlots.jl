using GLMakie
using DataFrames
using CSV
using Loess

const FIG_SIZE = (1200,900)

include("PlotDefaults.jl")

mutable struct RimData
    name::String
    tes::DataFrame
    uth::DataFrame
    x::Observable
    y::Observable
end
function importfiles(filedir)
    te_df = DataFrame(CSV.File(joinpath(filedir,"TraceElementData.csv")))
    uth_df = DataFrame(CSV.File(joinpath(filedir,"UPbData.csv")))
    uth_df[!,:Distance] = uth_df[!,:DistanceMod]
    selections = DataFrame(CSV.File(joinpath(filedir,"TraceElementsAvgs.csv")))
    rim_te = DataFrame[]
    rim_uth = DataFrame[]
    for row in eachrow(selections)
        push!(rim_te,select_range(te_df,row[:x1],row[:x2]))
        push!(rim_uth,select_range(uth_df,row[:x1],row[:x2]))
    end

    return rim_te, rim_uth
end

function select_range(df,x1,x2)
    rimdf = df
    if x1 < x2
        rimdf = filter(:Distance => x -> x1 <= x <= x2, df)
    
    elseif x1 > x2
        rimdf = filter(:Distance => x -> x2 <= x <= x1, df)

        d_array = copy(rimdf[!,:Distance])
        reverse!(rimdf)
        rimdf[!,:Distance] = d_array
    end

    #set first to x = 0
    xmin = minimum(rimdf[!,:Distance])
    rimdf[!,:Distance] = rimdf[!,:Distance] .- xmin
    arbitrary_scale!(rimdf)

    return rimdf
end

function arbitrary_scale!(df;x1 = 0.0, x2 = 100.0)
    increment = (x2 - x1)/(nrow(df)-1)

    x_arb = collect(x1:increment:x2)

    df[!,:x_arb] = x_arb

end

function import_all(directory, is_scaled)
    entries = readdir(directory)

    zrns = RimData[]
    xmax = 0
    for entry in entries
        if isdir(joinpath(directory,entry))
            name = entry
            tes, uth = importfiles(joinpath(directory,entry))
            
            for i in 1:lastindex(uth)
                rimName = name*"_R"*string(i)
                y = uth[i][!,:U_rescaled]

                x = uth[i][!,:Distance]
                if is_scaled
                    x = uth[i][!,:x_arb]
                end
                if maximum(x) > xmax
                    xmax = maximum(x)
                end
                zrn = RimData(rimName,tes[i],uth[i],Observable(x),Observable(y))
                push!(zrns,zrn)
            end
            
        end
    end
    
    return zrns, xmax
end



function changeElem!(rim, elem, is_scaled)

    if elem == "U/Th"
        rim.y[] = rim.uth[!,Symbol("Final U/Th")]
        if is_scaled
            rim.x[] = rim.uth[!,:x_arb]
        else
            rim.x[] = rim.uth[!,:Distance]
        end
    
    elseif elem == "U" || elem == "Th"
        rim.y[] = rim.uth[!,Symbol(elem*"_rescaled")]
        if is_scaled
            rim.x[] = rim.uth[!,:x_arb]
        else
            rim.x[] = rim.uth[!,:Distance]
        end
    else
        y = copy(rim.tes[!,Regex(elem*"\\d+")][!,1])
        x = copy(rim.tes[!,:Distance])
        if is_scaled
            x = copy(rim.tes[!,:x_arb])
        end
        if length(x) < length(rim.x[])
            x = [x;fill(x[end],length(rim.x[])-length(x))]
            y = [y;fill(y[end],length(rim.y[])-length(y))]
        end
        rim.x[] = x
        rim.y[] = y
        
    end


end

function changeElem!(rims, elem, ax, is_scaled)
    if elem == "U/Th"
        ax.ylabel = elem
    else
        ax.ylabel = elem*" (μg/g)"
    end
    maxy = 0
    for rim in rims
        changeElem!(rim,elem, is_scaled)

        if maximum(rim.y[])> maxy
            maxy = maximum(rim.y[])
        end
    end

    ylims!(ax,0,maxy)
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
    
    s = 0.1
    if length(x) < 100
        s = 0.2
    elseif length(x) < 50
        s = 0.5
    end

    model = loess(x,y,span=s)

    new_y = predict(model,rim.x[])

    rim.y[] = new_y


end


function init_fig(directory;is_scaled = true)

    fig = Figure(size = FIG_SIZE); display(GLMakie.Screen(),fig)

    ax = Axis(fig[1,2])

    zrns, xmax = import_all(directory, is_scaled)


    colsize!(fig.layout,1,Fixed(300))#So dropdown doesnt take too much space
    colsize!(fig.layout,2,Aspect(1,1)) #To enforcce square plot
    rowsize!(fig.layout,1,700)
    ax.xlabel = "x (μm)"
    if is_scaled
        ax.xlabel = "x (arbitrary scale)"
    end
    ax.ylabel = "U (μg/g)"

    for zrn in zrns
        lines!(ax,zrn.x,zrn.y,linewidth=3, label = zrn.name)
    end

    xMax = round(xmax,sigdigits=2)
    xInterval = 40
    if xmax <240
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
    end
    ax.xticks = 0:xInterval:xMax
    xlims!(ax,0,xmax)

    elemList = Array{String}([])
    headernames = names(zrns[1].tes)
    for i=3:lastindex(headernames)
        elemName = match(r"[A-z]+",headernames[i]).match
        if length(elemName) < 3
            push!(elemList,elemName)
        end
    end
    push!(elemList, "U/Th")
    menu = Menu(fig, options = elemList, default = "U")
    fig[1,1][1,1]=vgrid!(
        Label(fig,"Element:",width=nothing),
        menu,tellheight=false
    )
    on(menu.selection) do s
        changeElem!(zrns,s,ax, is_scaled)
    end

    smoothbutton = Button(fig[1,1][2,1],label = "Smooth Data")
    
    on(smoothbutton.clicks) do clicks
        for rim in zrns
            smooth!(rim,0.8)
        end
        
    end


    fig[1,3] = Legend(fig,ax)
end


GLMakie.activate!()
set_theme!(myTheme)
zrns = init_fig("ZrnLaserPlots/20SD06/TransectScaling/")

