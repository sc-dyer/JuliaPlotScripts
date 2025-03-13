using CairoMakie
# using GLMakie
using DataFrames
using CSV
using Peaks
using KernelDensity
using Isoplot
colourchoice = 2
include("PlotDefaults.jl")
DISCORDANCE_LOWER_LIMIT = -3.0
DISCORDANCE_UPPER_LIMIT = 5.0

function filterdiscordance(df::DataFrame)
    concordant = filter(:discordance=> x -> DISCORDANCE_LOWER_LIMIT <= x <=DISCORDANCE_UPPER_LIMIT, df)
    discordant = filter(:discordance=> x -> x < DISCORDANCE_LOWER_LIMIT || x > DISCORDANCE_UPPER_LIMIT, df)
    return concordant, discordant
end

function make_concordia_plot!(ax::Axis,coredf::DataFrame,rimdf::DataFrame,t1,t2)
    concordiacurve!(ax,t1,t2,fontsize = 15)
    ax.xlabel =rich(superscript("207"),"Pb/",superscript("235"),"U")
    ax.ylabel =rich(superscript("206"),"Pb/",superscript("238"),"U")
    rim_concordant, rim_discordant = filterdiscordance(rimdf)
    core_concordant, core_discordant = filterdiscordance(coredf)
    discordant = vcat(rim_discordant,core_discordant)
    # @show UPbAnalysis.(rim_discordant[!,:Pb207U235],rim_discordant[!,:Pb207U235_se1],
    # rim_concordant[!,:Pb206U238],rim_discordant[!,:Pb206U238_se1],rim_discordant[!,:cor_0735_0638])
    # @show UPbAnalysis.(core_discordant[!,:Pb207U235],core_discordant[!,:Pb207U235_se1],
    # core_discordant[!,:Pb206U238],core_discordant[!,:Pb206U238_se1],core_discordant[!,:cor_0735_0638])
    rim_analyses = UPbAnalysis.(rim_concordant[!,:Pb207U235],rim_concordant[!,:Pb207U235_se1],
                                rim_concordant[!,:Pb206U238],rim_concordant[!,:Pb206U238_se1],rim_concordant[!,:cor_0735_0638])
    core_analyses = UPbAnalysis.(core_concordant[!,:Pb207U235],core_concordant[!,:Pb207U235_se1],
                                core_concordant[!,:Pb206U238],core_concordant[!,:Pb206U238_se1],core_concordant[!,:cor_0735_0638])
    discordant_analyses = UPbAnalysis.(discordant[!,:Pb207U235],discordant[!,:Pb207U235_se1],
                            discordant[!,:Pb206U238],discordant[!,:Pb206U238_se1],discordant[!,:cor_0735_0638])
         
    plot!.(ax,discordant_analyses, color = (:grey,0.2))
    plot!.(ax,core_analyses,color = (colourscheme[1], 0.5),label ="Cores")
    plot!.(ax,rim_analyses,color = (colourscheme[2], 0.5),label = "Rims")
    
    
    
    
    
end

function make_concordia_plot!(ax::Axis,pegdf::DataFrame,t1,t2;calcdiscordia = false)
    concordiacurve!(ax,t1,t2,fontsize = 15)
    ax.xlabel =rich(superscript("207"),"Pb/",superscript("235"),"U")
    ax.ylabel =rich(superscript("206"),"Pb/",superscript("238"),"U")
    if !calcdiscordia
        peg_concordant, peg_discordant = filterdiscordance(pegdf)
        discordant =peg_discordant
        # @show UPbAnalysis.(rim_discordant[!,:Pb207U235],rim_discordant[!,:Pb207U235_se1],
        # rim_concordant[!,:Pb206U238],rim_discordant[!,:Pb206U238_se1],rim_discordant[!,:cor_0735_0638])
        # @show UPbAnalysis.(core_discordant[!,:Pb207U235],core_discordant[!,:Pb207U235_se1],
        # core_discordant[!,:Pb206U238],core_discordant[!,:Pb206U238_se1],core_discordant[!,:cor_0735_0638])
        peg_analyses = UPbAnalysis.(peg_concordant[!,:Pb207U235],peg_concordant[!,:Pb207U235_se1],
                                    peg_concordant[!,:Pb206U238],peg_concordant[!,:Pb206U238_se1],peg_concordant[!,:cor_0735_0638])
        discordant_analyses = UPbAnalysis.(discordant[!,:Pb207U235],discordant[!,:Pb207U235_se1],
                                discordant[!,:Pb206U238],discordant[!,:Pb206U238_se1],discordant[!,:cor_0735_0638])
            
        plot!.(ax,discordant_analyses, color = (:grey,0.2))
        plot!.(ax,peg_analyses,color = (colourscheme[1], 0.5)) 
    else
        peg_analyses = UPbAnalysis.(pegdf[!,:Pb207U235],pegdf[!,:Pb207U235_se1],
                                    pegdf[!,:Pb206U238],pegdf[!,:Pb206U238_se1],pegdf[!,:cor_0735_0638])

        upper = upperintercept(peg_analyses)
        lower = lowerintercept(peg_analyses)
        concordialine!(ax,lower,upper)
        plot!.(ax,peg_analyses,color = (colourscheme[1], 0.5))

    end
    # plot!.(ax,rim_analyses,color = (colourscheme[2], 0.5),label = "Rims")
end

function inset_kde!(ax::Axis,coredf::DataFrame,rimdf::DataFrame; plotcore = true, plotrim = true)
    rim_concordant, rim_discordant = filterdiscordance(rimdf)
    core_concordant, core_discordant = filterdiscordance(coredf)
    
    dates::Array{Float64} = rim_concordant[!,:Pb206U238Age]
    if plotcore == true
        dates = [dates;core_concordant[!,:Pb206U238Age]]
        if plotrim == false
            dates = core_concordant[!,:Pb206U238Age]
        end
    end
    density!(ax,dates, color = myColours[3])
    kernal = kde(dates)
    peaks = findmaxima(kernal.density)
    peak_ages = getindex(kernal.x,peaks.indices)
    scatter!(ax,peak_ages,peaks.heights,marker = :vline)
    text!(ax,peak_ages,peaks.heights,text=string.(Int.(round.(peak_ages))),color=:black,fontsize=8,offset = (5,0))
    ylims!(ax,0,1.1*maximum(kernal.density))
end

function inset_kde!(ax::Axis,pegdf::DataFrame;calcdiscordia = false)
    dates::Array{Float64} = []
    if !calcdiscordia
       peg_concordant, peg_discordant = filterdiscordance(pegdf)
       dates = peg_concordant[!,:PbPbAge]
    else
        dates = pegdf[!,:PbPbAge]
    end
    density!(ax,dates, color = colourscheme[3])
    kernal = kde(dates)
    peaks = findmaxima(kernal.density)
    peak_ages = getindex(kernal.x,peaks.indices)
    scatter!(ax,peak_ages,peaks.heights,marker = :vline)
    text!(ax,peak_ages,peaks.heights,text=string.(Int.(round.(peak_ages))),color=:black,fontsize=8,offset = (5,0))
    ylims!(ax,0,1.1*maximum(kernal.density))
end



# GLMakie.activate!()
function zrnpaper_plots()
    set_theme!(myTheme)
    zrn20SD06c = DataFrame(CSV.File("ZrnLaserPlots/UPb_91500-com/20SD06_cores.csv"))
    zrn20SD06r = DataFrame(CSV.File("ZrnLaserPlots/UPb_91500-com/20SD06_rims.csv"))

    zrn20SD17Ac = DataFrame(CSV.File("ZrnLaserPlots/UPb_91500-com/20SD17A_cores.csv"))
    zrn20SD17Ar = DataFrame(CSV.File("ZrnLaserPlots/UPb_91500-com/20SD17A_rims.csv"))

    zrn21SD68c = DataFrame(CSV.File("ZrnLaserPlots/UPb_91500-com//21SD68_cores.csv"))
    zrn21SD68r = DataFrame(CSV.File("ZrnLaserPlots/UPb_91500-com//21SD68_rims.csv"))

    zrn22SD55Ec = DataFrame(CSV.File("ZrnLaserPlots/UPb_91500-com/22SD55E_cores.csv"))
    zrn22SD55Er = DataFrame(CSV.File("ZrnLaserPlots/UPb_91500-com/22SD55E_rims.csv"))

    fig = Figure(;size = (800,800))
    rimdf = [zrn20SD06r zrn20SD17Ar; zrn21SD68r zrn22SD55Er]
    coredf = [zrn20SD06c zrn20SD17Ac; zrn21SD68c zrn22SD55Ec]

    zrndf = [zrn20SD06r, zrn20SD17Ar, zrn21SD68r, zrn22SD55Er, zrn20SD06c, zrn20SD17Ac, zrn21SD68c, zrn22SD55Ec]
    labels = ["20SD06" "20SD17A";"21SD68" "22SD55E"]
    my_analyses =Array{Isoplot.Analysis}([])
    for df in zrndf
        analysis = UPbAnalysis.(df[!,:Pb207U235],df[!,:Pb207U235_se1],df[!,:Pb206U238],df[!,:Pb206U238_se1],df[!,:cor_0735_0638])
        my_analyses= [my_analyses;analysis]
    end

    xmin, xmax, ymin, ymax = datalimits(my_analyses)

    for i in 1:2

        for j in 1:2
            
            isoplot_ax = Axis(fig[i,j])

            make_concordia_plot!(isoplot_ax,coredf[i,j],rimdf[i,j],900,1700)

            limits!(isoplot_ax,xmin,xmax,ymin,ymax)
            if j == 2
                hideydecorations!(isoplot_ax, ticks = false,minorticks=false)
            end
            if i == 1
                hidexdecorations!(isoplot_ax,ticks = false,minorticks=false)
            end

            text!(isoplot_ax, 0.05,0.87, text = labels[i,j],fontsize=24,space = :relative,offset = (4, -2))
            # text!(isoplot_ax,1.75,0.25, text = labels[i,j],fontsize=24)
            axislegend(isoplot_ax, merge = true,labelsize = 12,framevisible = false,halign =0.02,valign=0.85)
            density_ax = Axis(fig[i,j],
                            title = rich(superscript("206"),"Pb/",superscript("238"),"U age KDE"),
                            titlesize = 10,
                            aspect = 1,
                            xticklabelsize = 8,
                            yticklabelsvisible = false,
                            yticklabelsize = 8,
                            xlabelsize = 14,
                            ylabelsize = 14,
                            yticksvisible = false,
                            xminorticksvisible = true,
                            yminorticksvisible = false,
                            xticksmirrored = false,
                            xticksize = 4,
                            xtickwidth = 1.5,
                            yticksize = 6,
                            ytickwidth = 1.5,
                            xminorticksize = 4,
                            xminortickwidth = 1.5,
                            yminorticksize = 3,
                            yminortickwidth = 1.5,
                            width =Relative(0.4), 
                            height = Relative(0.4),
                            halign = 0.9, 
                            valign = 0.15, 
                            backgroundcolor=:white)

            if i == 2 && j == 1
                inset_kde!(density_ax,coredf[i,j],rimdf[i,j], plotrim = false)
            else
                inset_kde!(density_ax,coredf[i,j],rimdf[i,j])
            end
            xlims!(density_ax,800,1600)
            
        end
        
    end

    for i in 1:2
        colsize!(fig.layout,i,Aspect(1,1.0))
    end
    rowgap!(fig.layout,0)

    # display(GLMakie.Screen(),fig)
    # dates = [zrn20SD06c[!,:PbPbAge];zrn20SD06r[!,:PbPbAge]]
    # density!(density_ax,dates, color = colourscheme[3])
    # kernal = kde(dates)
    # peaks = findmaxima(kernal.density)
    # peak_ages = getindex(kernal.x,peaks.indices)
    # scatter!(density_ax,peak_ages,peaks.heights,markercolor = :black, marker = :vline)
    # text!(density_ax,peak_ages,peaks.heights,text=string.(Int.(round.(peak_ages))),color=:black,fontsize=8,offset = (5,0))
    # ylims!(density_ax,0,1.1*maximum(kernal.density))
    # scatter!(density_ax,maxima_helper.())
    CairoMakie.save("ConcordiaPlots/UPb_Plots_v9.svg",fig)
    # display(GLMakie.Screen(),fig)

    # GLMakie.save("ConcordiaPlots/UPb_Plots.png",fig,px_per_unit = 4)
end

function fluidfluxpaper_plots()
    set_theme!(myTheme)
    zrn22SD55D = DataFrame(CSV.File("../../LAICPMS/PegDates/22SD55D.csv"))

    zrn23SD02C = DataFrame(CSV.File("../../LAICPMS/PegDates/23SD02C.csv"))

    zrn23SD03E = DataFrame(CSV.File("../../LAICPMS/PegDates/23SD03E.csv"))

    zrn23SD03K = DataFrame(CSV.File("../../LAICPMS/PegDates/23SD03K.csv"))

    zrn23SD20E = DataFrame(CSV.File("../../LAICPMS/PegDates/23SD20E.csv"))
    fig = Figure(;size = (800,800))
    pegdf = [zrn22SD55D zrn23SD02C; zrn23SD03E zrn23SD03K]
   
    zrndf = [zrn22SD55D, zrn23SD02C, zrn23SD03E, zrn23SD03K,zrn23SD20E]
    labels = ["22SD55D" "23SD02C";"23SD03E" "23SD03K"]
    my_analyses =Array{Isoplot.Analysis}([])
    for df in zrndf
        analysis = UPbAnalysis.(df[!,:Pb207U235],df[!,:Pb207U235_se1],df[!,:Pb206U238],df[!,:Pb206U238_se1],df[!,:cor_0735_0638])
        my_analyses= [my_analyses;analysis]
    end

    xmin, xmax, ymin, ymax = datalimits(my_analyses)

    for i in 1:2

        for j in 1:2
            
            isoplot_ax = Axis(fig[i,j])

            make_concordia_plot!(isoplot_ax,pegdf[i,j],400,2200)

            limits!(isoplot_ax,xmin,xmax,ymin,ymax)
            if j == 2
                hideydecorations!(isoplot_ax, ticks = false,minorticks=false)
            end
            if i == 1 || j == 1
                hidexdecorations!(isoplot_ax,ticks = false,minorticks=false)
            end

            text!(isoplot_ax, 0.05,0.87, text = labels[i,j],fontsize=24,space = :relative,offset = (4, -2))
            # text!(isoplot_ax,1.75,0.25, text = labels[i,j],fontsize=24)
            # axislegend(isoplot_ax, merge = true,labelsize = 12,framevisible = false,halign =0.02,valign=0.85)
            
        end
        
    end
    isoplot_ax = Axis(fig[3,1])

    make_concordia_plot!(isoplot_ax,zrn23SD20E,400,2200,calcdiscordia=true)

    limits!(isoplot_ax,xmin,xmax,ymin,ymax)

    text!(isoplot_ax, 0.05,0.87, text = "23SD20E",fontsize=24,space = :relative,offset = (4, -2))
    # text!(isoplot_ax,1.75,0.25, text = labels[i,j],fontsize=24)
    # axislegend(isoplot_ax, merge = true,labelsize = 12,framevisible = false,halign =0.02,valign=0.85)

    for i in 1:2
        colsize!(fig.layout,i,Aspect(1,1.0))
    end
    rowgap!(fig.layout,0)

    # display(GLMakie.Screen(),fig)
    # dates = [zrn20SD06c[!,:PbPbAge];zrn20SD06r[!,:PbPbAge]]
    # density!(density_ax,dates, color = colourscheme[3])
    # kernal = kde(dates)
    # peaks = findmaxima(kernal.density)
    # peak_ages = getindex(kernal.x,peaks.indices)
    # scatter!(density_ax,peak_ages,peaks.heights,markercolor = :black, marker = :vline)
    # text!(density_ax,peak_ages,peaks.heights,text=string.(Int.(round.(peak_ages))),color=:black,fontsize=8,offset = (5,0))
    # ylims!(density_ax,0,1.1*maximum(kernal.density))
    # scatter!(density_ax,maxima_helper.())
    CairoMakie.save("ConcordiaPlots/Peg_UPb_v1.svg",fig)
    # display(GLMakie.Screen(),fig)

    # GLMakie.save("ConcordiaPlots/UPb_Plots.png",fig,px_per_unit = 4)
end

# fluidfluxpaper_plots()
zrnpaper_plots()