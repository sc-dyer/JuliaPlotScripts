# using GLMakie
using CairoMakie
using DataFrames
using CSV
using Statistics
using Loess

const FIG_SIZE = (1200,900)
const SPECIAL_SCALE = 40
colourchoice = 2
include("PlotDefaults.jl")

mutable struct RimData
    name::String
    tes::DataFrame
    uth::DataFrame
    x::Observable
    y::Observable
end
function importfiles(filedir; special = false)
    te_df = DataFrame(CSV.File(joinpath(filedir,"TraceElementData.csv")))
    uth_df = DataFrame(CSV.File(joinpath(filedir,"UPbData.csv")))
    uth_df[!,:Distance] = uth_df[!,:DistanceMod]
    selections = DataFrame(CSV.File(joinpath(filedir,"TraceElementsAvgs.csv")))
    rim_te = DataFrame[]
    rim_uth = DataFrame[]
    for row in eachrow(selections)
        if special
            push!(rim_te,select_range(te_df,row[:x1],row[:x2], x0 = row[:x0], special = special))
            push!(rim_uth,select_range(uth_df,row[:x1],row[:x2], x0 = row[:x0], special = special))
        else
            push!(rim_te,select_range(te_df,row[:x1],row[:x2],x0 = row[:x0]))
            push!(rim_uth,select_range(uth_df,row[:x1],row[:x2],x0 = row[:x0]))
        end
    end

    return rim_te, rim_uth
end

function select_range(df,x1,x2;x0 = -1,special = false)
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
    if x0 > 0
        if x1 < x2
            xmin = x1
        else
            #distance was reversed so x1 is no longer where the rim start is
            xmin = (x0-x1) + x2
            
        end

    end
    # @show rimdf[!,:Distance]
    rimdf[!,:Distance] = rimdf[!,:Distance] .- xmin
    # @show rimdf[!,:Distance]
    if special
        arbitrary_scale!(rimdf, x2 = SPECIAL_SCALE)
    else
        arbitrary_scale!(rimdf)
    end
    return rimdf
end

function arbitrary_scale!(df;x1 = 0.0, x2 = 100.0)
    # increment = (x2 - x1)/(nrow(df)-1)
    # x0 = x1
    # if df[1,:Distance] < 0
    #     i = 1
    #     while df[i,:Distance] < 0
    #         i += 1
    #     end
        
    #     increment = (x2 - x1)/(nrow(df)-i)
    #     x0 = 0.0 - (i-1)*increment
        
    # end
    
    
    # x_arb = [collect(x0:increment:x1);collect(x1+increment:increment:x2)]
   
    ratio = (x2 - x1)/maximum(df[!,:Distance])
    x_arb = df[!,:Distance] .* ratio
    df[!,:x_arb] = x_arb

end

function import_all(directory, is_scaled; special = nothing)
    entries = readdir(directory)

    zrns = RimData[]
    
    for entry in entries
        if isdir(joinpath(directory,entry))
            name = entry
            tes = DataFrame()
            uth = DataFrame()
            if name == special
                tes, uth = importfiles(joinpath(directory,entry), special = true)
            else
                tes, uth = importfiles(joinpath(directory,entry))
            end

            for i in 1:lastindex(uth)
                rimName = name*"-R"*string(i)
                y = uth[i][!,:U_rescaled]

                x = uth[i][!,:Distance]
                if is_scaled
                    x = uth[i][!,:x_arb]
                end
               
                zrn = RimData(rimName,tes[i],uth[i],Observable(x),Observable(y))
                push!(zrns,zrn)
            end
            
        end
    end
    
    return zrns
end



function changeElem!(rim, elem, is_scaled)

    if elem == "Th/U"
        rim.y[] = 1.0 ./ rim.uth[!,Symbol("Final U/Th")]
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

    elseif contains(elem, "age")
        rim.y[] = rim.uth[!,Symbol(elem)]
        if is_scaled
            rim.x[] = rim.uth[!,:x_arb]
        else
            rim.x[] = rim.uth[!,:Distance]
        end
    elseif contains(elem,"/")
        y = copy(rim.tes[!,elem])
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
    if contains(elem, "age")
        ax.ylabel = elem*" (Ma)"
    elseif contains(elem,"/")
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
        s = 0.5
    elseif length(x) < 50
        s = 0.8
    end

    model = loess(x,y,span=s)

    new_y = predict(model,rim.x[])

    rim.y[] = new_y


end

function average_interval(rim, interval)
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



    neg_intervals = convert(Int64,ceil((-minimum(x))/interval))
    pos_intervals = convert(Int64,ceil(maximum(x)/interval))
    #    @show maximum(x), minimum(x), num_intervals 
    x_mids = Float64[]
    yavgs = Float64[]
    ystds = Float64[]
    start = minimum(x)
    if start < 0
        for i in 1:neg_intervals
            fval = (i-1)*interval + start
            lval  = interval*i + start
            # @show fval, lval
            if lval > 0
                lval = 0
            end
            indices = fval .<= x .< lval
            
            yavg = mean(y[indices])
            ystd = std(y[indices])
            x_mid = mean(x[indices])
            if !isnan(x_mid)
                push!(yavgs,yavg)
                push!(ystds,ystd)
                push!(x_mids,x_mid)
            end
        end
    end

    for i in 1:pos_intervals
        fval = (i-1)*interval
        lval  = interval*i
        # @show fval, lval
        indices = fval .<= x .< lval
        
        yavg = mean(y[indices])
        ystd = std(y[indices])
        x_mid = mean(x[indices])

        push!(yavgs,yavg)
        push!(ystds,ystd)
        push!(x_mids,x_mid)
    end
    # @show x_mids, yavgs
    return x_mids, yavgs, ystds
end

function plot_average_interval!(ax,rims,interval,elem; is_scaled = true, limits = true, scaling = nothing)
    
    
    lines!(ax, [0,0], [-1000,10000],linestyle = :dash, color = :black, linewidth=2)

    ax.xlabel = "x (μm)"
    if is_scaled
        ax.xlabel = "x (arbitrary scale)"
    end
    changeElem!(rims,elem,ax, is_scaled)
    

    for rim in rims
        x, y, y_err = average_interval(rim,interval)
        # @show x, y, y_err
        errorbars!(ax, x,y, y_err, whiskerwidth = 6)
        scatterlines!(ax,x,y,linewidth=2, label = rim.name)
        
    end

    maxy = 0
    minx = Inf
    maxx = -Inf
    for rim in rims

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

    if !isnothing(scaling)
        ylims!(ax,scaling[1],scaling[2])
    end

    if limits
        ylims!(ax,0,maxy)
        xlims!(ax,minx,maxx)
    end
    
end


function init_fig(directory;is_scaled = true, special = nothing)

    fig = Figure(size = FIG_SIZE); display(GLMakie.Screen(),fig)

    ax = Axis(fig[1,2])
    lines!(ax, [0,0], [-1000,10000],linestyle = :dash, color = :black, linewidth=3)
    zrns = import_all(directory, is_scaled,special=special)


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

    
    # xMax = round(xmax,sigdigits=2)
    # xInterval = 40
    # if xmax <240
    #     xInterval = 20
    #     if xMax <120
    #         xInterval = 10
    #         if xMax <= 60
    #             xInterval = 5
    #             if xMax <= 25
    #                 xInterval = 2
    #             end
    #         end
    #     end    
    # end
    # ax.xticks = 0:xInterval:xMax
    # xlims!(ax,0,xmax)
   
    elemList = Array{String}([])
    headernames = names(zrns[1].tes)
    for i=3:lastindex(headernames)
        elemName = match(r"[A-z]+",headernames[i]).match
        if length(elemName) < 3
            push!(elemList,elemName)
        end
    end
    push!(elemList, "Th/U")
    menu = Menu(fig, options = elemList, default = "U")
    fig[1,1][1,1]=vgrid!(
        Label(fig,"Element:",width=nothing),
        menu,tellheight=false
    )
    elem_selection = "U"
    on(menu.selection) do s
        changeElem!(zrns,s,ax, is_scaled)
        elem_selection = s
    end

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

    ylims!(ax,0,maxy)
    xlims!(ax,minx,maxx)
    smoothbutton = Button(fig[1,1][2,1],label = "Smooth Data")
    
    on(smoothbutton.clicks) do clicks
        for rim in zrns
            smooth!(rim,0.8)
        end
        
    end

    avgbutton = Button(fig[1,1][3,1],label = "Average Intervals")
    
    on(avgbutton.clicks) do clicks
        fig2 = Figure(size = FIG_SIZE)
        ax2 = plot_average_interval!(zrns,10,elem_selection, is_scaled=is_scaled)
        fig2[1,1] = ax2
        fig2[1,2] = Legend(fig2,ax2)
        display(GLMakie.Screen(),fig2)
    end

    fig[1,3] = Legend(fig,ax)
end


function save_fig(directory;is_scaled = true, special = nothing, elem = "U",scaling = nothing, averagelines = false)


    zrns = import_all(directory, is_scaled,special=special)
    fig = Figure(size = FIG_SIZE)
    if averagelines
        ax = plot_average_interval!(zrns,10,elem,is_scaled = is_scaled, limits = false, scaling = scaling)
        fig[1,1] = ax
        fig[1,2] = Legend(fig,ax)
        
    else
        

        ax = Axis(fig[1,1])
        lines!(ax, [0,0], [-1000,10000],linestyle = :dash, color = :black, linewidth=3)
        

        
    
        ax.xlabel = "x (μm)"
        if is_scaled
            ax.xlabel = "x (arbitrary scale)"
        end

        if elem == "Th/U"
            ax.ylabel = elem
            
        else
            ax.ylabel = elem*" (μg/g)"
        end
        changeElem!(zrns,elem,ax,is_scaled)
        
        for zrn in zrns
            
            lines!(ax,zrn.x,zrn.y,linewidth=3, label = zrn.name)
        end

        
        # xMax = round(xmax,sigdigits=2)
        # xInterval = 40
        # if xmax <240
        #     xInterval = 20
        #     if xMax <120
        #         xInterval = 10
        #         if xMax <= 60
        #             xInterval = 5
        #             if xMax <= 25
        #                 xInterval = 2
        #             end
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

        ylims!(ax,0,maxy)
        # xlims!(ax,minx,maxx)
    

        if !isnothing(scaling)
            ylims!(ax,scaling[1],scaling[2])
        end
        fig[1,2] = Legend(fig,ax)
    end

    return fig
end

function make_ax!(ax,directory;is_scaled = true, special = nothing, elem = "U",scaling = nothing, averagelines = false)


    zrns = import_all(directory, is_scaled,special=special)
    
    if averagelines
        plot_average_interval!(ax, zrns,10,elem,is_scaled = is_scaled, limits = false, scaling=scaling)
    else
        

        
        lines!(ax, [0,0], [-1000,10000],linestyle = :dash, color = :black, linewidth=2)
        

        
    
        ax.xlabel = "x (μm)"
        if is_scaled
            ax.xlabel = "x (arbitrary scale)"
        end

        if elem == "Th/U" 
            ax.ylabel = elem
        elseif contains(elem,"age")
            ax.ylabel = elem*" (Ma)"
        else
            ax.ylabel = elem*" (μg/g)"
        end
        changeElem!(zrns,elem,ax,is_scaled)
        
        for zrn in zrns
            lines!(ax,zrn.x,zrn.y,linewidth=2, label = zrn.name)
        end

        
        # xMax = round(xmax,sigdigits=2)
        # xInterval = 40
        # if xmax <240
        #     xInterval = 20
        #     if xMax <120
        #         xInterval = 10
        #         if xMax <= 60
        #             xInterval = 5
        #             if xMax <= 25
        #                 xInterval = 2
        #             end
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

        ylims!(ax,0,maxy)
        # xlims!(ax,minx,maxx)
       

        if !isnothing(scaling)
            ylims!(ax,scaling[1],scaling[2])
        end
        
    end
    ax.xticklabelsvisible = false
    # axislegend(ax,framevisible = false)
    
end

function make_ratio_ax!(ax,directory;is_scaled = true, special = nothing, num = "Dy",den="Yb",scaling = nothing)


    zrns = import_all(directory, is_scaled,special=special)
    
   ratio = num*"/"*den
    for zrn in zrns
     
        zrn.tes[!,ratio] = zrn.tes[!,Regex(num*"\\d+")][!,1]./zrn.tes[!,Regex(den*"\\d+")][!,1]
        # @show zrn.tes[!,ratio]
        filter!(Symbol(ratio) => r -> -Inf < r < Inf,zrn.tes)
        # @show zrn.tes[!,ratio]
    end
    
    
    lines!(ax, [0,0], [-1000,10000],linestyle = :dash, color = :black, linewidth=2)
    

    

    ax.xlabel = "x (μm)"
    if is_scaled
        ax.xlabel = "x (arbitrary scale)"
    end

    
    ax.ylabel = ratio
    
    changeElem!(zrns,ratio,ax,is_scaled)
    
    for zrn in zrns
        lines!(ax,zrn.x,zrn.y,linewidth=2, label = zrn.name)
    end

    
    # xMax = round(xmax,sigdigits=2)
    # xInterval = 40
    # if xmax <240
    #     xInterval = 20
    #     if xMax <120
    #         xInterval = 10
    #         if xMax <= 60
    #             xInterval = 5
    #             if xMax <= 25
    #                 xInterval = 2
    #             end
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

    ylims!(ax,0,maxy)
    # xlims!(ax,minx,maxx)
    

    if !isnothing(scaling)
        ylims!(ax,scaling[1],scaling[2])
    end
    
    
    ax.xticklabelsvisible = false
    # axislegend(ax,framevisible = false)
    
end

# GLMakie.activate!()
set_theme!(myTheme)
# zrns = init_fig("ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26")
# zrns = init_fig("ZrnLaserPlots/20SD17A/TransectScaling/",is_scaled=true)
# f = save_fig("ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "U")
# save(joinpath("ZrnLaserPlots/20SD06/TransectScaling2/","URimPlots.svg"),f)

# f = save_fig("ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "Th",scaling = (-10,50))
# save(joinpath("ZrnLaserPlots/20SD06/TransectScaling2/","ThRimPlots.svg"),f)

# f = save_fig("ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "Y",averagelines = true)
# save(joinpath("ZrnLaserPlots/20SD06/TransectScaling2/","YRimPlots.svg"),f)

# f = save_fig("ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "Yb",averagelines = true)
# save(joinpath("ZrnLaserPlots/20SD06/TransectScaling2/","YbRimPlots.svg"),f)

f = Figure()

# ax1 = Axis(f[1,1],aspect=1.5)
# make_ax!(ax1,"ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "U")
# ax1.xlabelvisible = false
# ax2 = Axis(f[1,2],aspect=1.5)
# make_ax!(ax2,"ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "Th",scaling = (-5,50))
# ax2.xlabelvisible = false
ax4 = Axis(f[1,1],aspect=1.5)
make_ax!(ax4,"ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "Y",averagelines = true)
ax4.xlabelvisible = false
ax5 = Axis(f[2,1],aspect=1.5)
make_ax!(ax5,"ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "Yb",averagelines = true)
# ax3 = Axis(f[1,2],aspect=1.5)
# make_ax!(ax3,"ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "Th/U")
# ax3.xlabelvisible = false

# ax1 = Axis(f[1,1],aspect=1.5)
# make_ax!(ax1,"ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "Final Pb206/U238 age",scaling = (900,1400),averagelines = true)
# ax1.xlabelvisible = false
# ax2 = Axis(f[2,1],aspect=1.5)
# make_ax!(ax2,"ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "Final Pb207/U235 age",scaling = (900,1400),averagelines = true)
# ax2.xlabelvisible = false
# ax3 = Axis(f[3,1],aspect=1.5)
# make_ax!(ax3,"ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", elem = "Final Pb207/Pb206 age",scaling = (900,1400),averagelines = true)


# ax6 = Axis(f[3,2],aspect=1.5)
# make_ratio_ax!(ax6,"ZrnLaserPlots/20SD06/TransectScaling2/",is_scaled=true, special = "Zrn26", num="Yb",den="Gd")
# colsize!()
rowsize!(f.layout,1,Fixed(300))
rowsize!(f.layout,2,Fixed(300))
# rowsize!(f.layout,3,Fixed(300))
colsize!(f.layout,1,Fixed(450))
# colsize!(f.layout,2,Fixed(450))

f[1:2,3] = Legend(f,ax1)
text!(ax4, 0.03,0.85,text="A",fontsize=36,space = :relative,offset = (4, 0))
text!(ax5, 0.03,0.85,text="B",fontsize=36,space = :relative,offset = (4, 0))
# text!(ax3, 0.03,0.85,text="C",fontsize=36,space = :relative,offset = (4, 0))
# text!(ax4, 0.03,0.85,text="D",fontsize=36,space = :relative,offset = (4, 0))
# text!(ax5, 0.03,0.85,text="E",fontsize=36,space = :relative,offset = (4, 0))

resize_to_layout!(f)


save(joinpath("ZrnLaserPlots/20SD06/TransectScaling2/","YREE.svg"),f)