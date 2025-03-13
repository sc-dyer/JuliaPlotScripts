using GLMakie
using TernaryDiagrams
using DataFrames
using CSV
using PetroBase
using Migmatites
using JLD2
colourchoice=2
include("../PlotDefaults.jl")
markers = [:circle,:rect,:utriangle,:cross]
# @recipe(TernaryText, x, y, z) do scene
#     Attributes(
#         text = "", 
#         marker = :circle,
#         markersize = 8,
#     )
# end


function solarbrowntern(maindir,sourcename,host,pindex)
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = 1.0)
    
    ternaryaxis!(
        ax; 
        labelx = "Na + Ca",
        labely = "K",
        labelz = "Fe* + Mg + Ti",
        label_fontsize = 24,
        hide_vertex_labels = true,
        tick_fontsize = 0,
        grid_line_width = 0.0,
        arrow_scale = 0.0
    
        # more options available, check out attributes with ?ternaryaxis (same for other plot functions)
        #= Note 
        Depending on the length of the axis labels, they may seem unaligned. 
        Use the kwarg arrow_label_rotation_adjustment to rotate them slightly. 
        For longer labels, use a value closer to 1 (trial and error it).
        =#
    )
    
    # the triangle is drawn from (0,0) to (0.5, sqrt(3)/2) to (1,0).
    xlims!(ax, -0.2, 1.2) # to center the triangle and allow space for the labels
    ylims!(ax, -0.3, 1.1)
    hidedecorations!(ax) # to hide the axis decorations
    
    for i in 1:lastindex(pindex)
        
        data750 = load(joinpath(maindir,sourcename*"_source_"*host*"_host_750/","outputs_x"*string(pindex[i])*".jld2"))
        data800 = load(joinpath(maindir,sourcename*"_source_"*host*"_host_800/","outputs_x"*string(pindex[i])*".jld2"))
       
        data = [data750, data800]
        temps = [750,800]
        
        for j in 1:2

            
            melts = getmelt.(data[j]["hosts"])
            
            melts = melts[mol.(melts) .> 0]
            meltcompos = getcompo.(melts)
            na = mol.(getchemical.(meltcompos,"Na2O")*2)
            k = mol.(getchemical.(meltcompos,"K2O")*2)
            ti = mol.(getchemical.(meltcompos,"TiO2"))
            ca = mol.(getchemical.(meltcompos,"CaO"))
            mg = mol.(getchemical.(meltcompos,"MgO"))
            fe = mol.(getchemical.(meltcompos,"FeO"))

            norm = na .+ ca .+ mg .+ fe .+ ti .+ k
            x = (na .+ ca)./norm
            y = k./norm
            z = (fe .+ ti .+ mg)./norm

            sourceh2o =  round(massfrac(getcompo(data[j]["system_steps"][1]),"H2O")*100,sigdigits = 2)
            grouplabel = string(sourceh2o)*" wt% H2o - "*string(temps[j])*"°C"
            ternaryscatter!(ax,x,y,z,label = grouplabel, color = myColours2[j],markersize = 6,strokewidth = 0.5,marker =markers[i])

            templabels = string.(data[j]["extract_T"]).*"°C"
            # ternarytext!(ax,x,y,z,text = templabels)
        end
    end
    
    # fig[1,2] = Legend(fig,ax)
    # save(joinpath(maindir,sourcename*"solarbrowntern.svg"),fig)
    return fig
end

function migtern(maindir,sourcename,host,pindex)
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = 1.0)
    
    ternaryaxis!(
        ax; 
        labelx = "Qz + Afs + Pl",
        labely = "Amp + Bt",
        labelz = "Melt",
        label_fontsize = 24,
        hide_vertex_labels = true,
        tick_fontsize = 0,
        grid_line_width = 0.0,
        arrow_scale = 0.0
    
        # more options available, check out attributes with ?ternaryaxis (same for other plot functions)
        #= Note 
        Depending on the length of the axis labels, they may seem unaligned. 
        Use the kwarg arrow_label_rotation_adjustment to rotate them slightly. 
        For longer labels, use a value closer to 1 (trial and error it).
        =#
    )
    
    # the triangle is drawn from (0,0) to (0.5, sqrt(3)/2) to (1,0).
    xlims!(ax, -0.2, 1.2) # to center the triangle and allow space for the labels
    ylims!(ax, -0.3, 1.1)
    hidedecorations!(ax) # to hide the axis decorations
    
    for i in 1:lastindex(pindex)
        
        data750 = load(joinpath(maindir,sourcename*"_source_"*host*"_host_750/","outputs_x"*string(pindex[i])*".jld2"))
        data800 = load(joinpath(maindir,sourcename*"_source_"*host*"_host_800/","outputs_x"*string(pindex[i])*".jld2"))
       
        data = [data750, data800]
        temps = [750,800]
        
        for j in 1:2

            
            melts = data[j]["melts_f"]
            
            qz = first.(getphase.(melts,"q"))
            afs = first.(getphase.(melts,"afs"))
            pl = first.(getphase.(melts,"pl"))
            
            melt = first.(getphase.(melts,"melt"))
            
            
            norm = vol.(qz) .+ vol.(afs) .+ vol.(pl) .+ vol.(melt)
            x = ( vol.(afs) .+ vol.(pl))./norm
            y = vol.(qz)./norm
            z = vol.(melt)./norm

            # norm = vol.(qz) .+ vol.(afs) .+ vol.(amp)
            # x = vol.(afs)./norm
            # y = vol.(qz)./norm
            # z = vol.(amp)./norm
            sourceh2o =  round(massfrac(getcompo(data[j]["system_steps"][1]),"H2O")*100,sigdigits = 2)
            grouplabel = string(sourceh2o)*" wt% H2O - "*string(temps[j])*"°C"
           
            
            ternaryscatter!(ax,x,y,z,label = grouplabel, color = myColours2[j],markersize = 10,strokewidth = 0.5,marker =markers[i])
            templabels = string.(round.(data[j]["extract_T"],sigdigits=3)).*"°C"
            
            ternarytext!(ax,x,y,z,text = templabels)
            
            

            

            # 
        end
    end
    
    fig[1,2] = Legend(fig,ax)
    # save(joinpath(maindir,sourcename*"solarbrowntern.svg"),fig)
    return fig
end

function barkertern(filename,savedir)
    fig = Figure();
    ax = Axis(fig[1, 1], aspect = 1.0);
    
    ternaryaxis!(
        ax; 
        labelx = "Ab",
        labely = "Or",
        labelz = "An",
        label_fontsize = 24,
        hide_vertex_labels = true,
        tick_fontsize = 0,
        grid_line_width = 0.0,
        arrow_scale = 0.0
    
        # more options available, check out attributes with ?ternaryaxis (same for other plot functions)
        #= Note 
        Depending on the length of the axis labels, they may seem unaligned. 
        Use the kwarg arrow_label_rotation_adjustment to rotate them slightly. 
        For longer labels, use a value closer to 1 (trial and error it).
        =#
    )
    
    # the triangle is drawn from (0,0) to (0.5, sqrt(3)/2) to (1,0).
    xlims!(ax, -0.2, 1.2) # to center the triangle and allow space for the labels
    ylims!(ax, -0.3, 1.1)
    hidedecorations!(ax) # to hide the axis decorations
    
    gdf = DataFrame(CSV.File(filename))
    
    gdf[!,:norm] = gdf[!,:albite] .+ gdf[!,:anorthite] .+ gdf[!,:orthoclase] 
    gdf[!,:ab] = gdf[!,:albite]./gdf[!,:norm]
    gdf[!,:an] = gdf[!,:anorthite]./gdf[!,:norm]
    gdf[!,:or] = gdf[!,:orthoclase]./gdf[!,:norm]

    ternarylines!(ax,[0,7,0.575],[0,0.25],[0.30,0.175],color = :black,linewidth=2)
    ternarylines!(ax,[0.70,0.575],[0.30,0.25],[0,0.175],color = :black,linewidth=2)
    ternarylines!(ax,[0.57,0.485],[0.25,0.35],[0.175,0.165],color = :black,linewidth=2)
    ternarylines!(ax,[0.48,0.33],[0.35,0.35],[0.165,0.32],color = :black,linewidth=2)
    ternarylines!(ax,[0.60,0.40],[0.20,0.20],[0.20,0.40],color=:black,linewidth=2)
    
    for i in 1:lastindex(SLAGSTAD_CATS)
        cat = SLAGSTAD_CATS[i]
        
        gdfiltered = filter(:Category => c -> c == cat, gdf)
        ternaryscatter!(ax,gdfiltered[!,:ab],gdfiltered[!,:or],gdfiltered[!,:an], label = cat, color = myColours2[i],markersize = 6,marker=:rect,strokewidth = 0.5)
    end
    for i in 1:lastindex(MY_CATS)
        cat = MY_CATS[i]
        
        gdfiltered = filter(:Category => c -> c == cat, gdf)
        ternaryscatter!(ax,gdfiltered[!,:ab],gdfiltered[!,:or],gdfiltered[!,:an], label = cat, color = myColours2[i],markersize = 6,marker=:circle,strokewidth = 0.5)
        ternarytext!(ax,gdfiltered[!,:ab],gdfiltered[!,:or],gdfiltered[!,:an],text = gdfiltered[!,:Sample],fontsize=2)
    end
    fig[1,2] = Legend(fig,ax)
    save(joinpath(savedir,"barkertern_text.svg"),fig)
end
f = migtern("MeltPathPlots/","23SD03B","23SD03B",[2,3,5,10])
# barkertern("../../Geochem/FluidFlux_CIPWnorm2.csv","TernPlots/")

