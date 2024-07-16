using CairoMakie
# using GLMakie
using DataFrames
using CSV
using Peaks
using KernelDensity
using Isoplot
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
    plot!.(ax,core_analyses,color = (myColours[1], 0.5),label ="Cores")
    plot!.(ax,rim_analyses,color = (myColours[2], 0.5),label = "Rims")
    
    
    
    
    
end

function inset_kde!(ax::Axis,coredf::DataFrame,rimdf::DataFrame)
    rim_concordant, rim_discordant = filterdiscordance(rimdf)
    core_concordant, core_discordant = filterdiscordance(coredf)
    dates = [core_concordant[!,:PbPbAge];rim_concordant[!,:PbPbAge]]
    density!(ax,dates, color = myColours[3])
    kernal = kde(dates)
    peaks = findmaxima(kernal.density)
    peak_ages = getindex(kernal.x,peaks.indices)
    scatter!(ax,peak_ages,peaks.heights,marker = :vline)
    text!(ax,peak_ages,peaks.heights,text=string.(Int.(round.(peak_ages))),color=:black,fontsize=8,offset = (5,0))
    ylims!(ax,0,1.1*maximum(kernal.density))
end

# GLMakie.activate!()
set_theme!(myTheme)
zrn20SD06c = DataFrame(CSV.File("ZrnLaserPlots/UPb/20SD06_cores.csv"))
zrn20SD06r = DataFrame(CSV.File("ZrnLaserPlots/UPb/20SD06_rims.csv"))

zrn20SD17Ac = DataFrame(CSV.File("ZrnLaserPlots/UPb/20SD17A_cores_ples.csv"))
zrn20SD17Ar = DataFrame(CSV.File("ZrnLaserPlots/UPb/20SD17A_rims_ples.csv"))

zrn21SD68c = DataFrame(CSV.File("ZrnLaserPlots/UPb/21SD68_cores_ples.csv"))
zrn21SD68r = DataFrame(CSV.File("ZrnLaserPlots/UPb/21SD68_rims_ples.csv"))

zrn22SD55Ec = DataFrame(CSV.File("ZrnLaserPlots/UPb/22SD55E_cores_ples.csv"))
zrn22SD55Er = DataFrame(CSV.File("ZrnLaserPlots/UPb/22SD55E_rims_ples.csv"))

fig = Figure(;size = (800,800))
rimdf = [zrn20SD06r zrn20SD17Ar; zrn21SD68r zrn22SD55Er]
coredf = [zrn20SD06c zrn20SD17Ac; zrn21SD68c zrn22SD55Ec]

zrndf = [zrn20SD06r, zrn20SD17Ar, zrn21SD68r, zrn22SD55Er, zrn20SD06c, zrn20SD17Ac, zrn21SD68c, zrn22SD55Ec]
labels = ["20SD06" "20SD17A";"21SD68" "22SD55E"]
my_analyses =Array{Isoplot.Analysis}([])
for df in zrndf
    analysis = UPbAnalysis.(df[!,:Pb207U235],df[!,:Pb207U235_se1],df[!,:Pb206U238],df[!,:Pb206U238_se1],df[!,:cor_0735_0638])
    global my_analyses= [my_analyses;analysis]
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
                        title = rich(superscript("207"),"Pb/",superscript("206"),"Pb age KDE"),
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

        inset_kde!(density_ax,coredf[i,j],rimdf[i,j])
        xlims!(density_ax,800,1600)
        
    end
    
end

for i in 1:2
    colsize!(fig.layout,i,Aspect(1,1.0))
end
rowgap!(fig.layout,0)

# display(GLMakie.Screen(),fig)
# dates = [zrn20SD06c[!,:PbPbAge];zrn20SD06r[!,:PbPbAge]]
# density!(density_ax,dates, color = myColours[3])
# kernal = kde(dates)
# peaks = findmaxima(kernal.density)
# peak_ages = getindex(kernal.x,peaks.indices)
# scatter!(density_ax,peak_ages,peaks.heights,markercolor = :black, marker = :vline)
# text!(density_ax,peak_ages,peaks.heights,text=string.(Int.(round.(peak_ages))),color=:black,fontsize=8,offset = (5,0))
# ylims!(density_ax,0,1.1*maximum(kernal.density))
# scatter!(density_ax,maxima_helper.())
CairoMakie.save("ConcordiaPlots/UPb_Plots_v4.svg",fig)
# display(GLMakie.Screen(),fig)

# GLMakie.save("ConcordiaPlots/UPb_Plots.png",fig,px_per_unit = 4)