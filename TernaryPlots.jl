using CairoMakie
using TernaryDiagrams
using DataFrames
using CSV
const MY_CATS = ["Host","Selvage","Pegmatite"]
const SLAGSTAD_CATS = ["Granodiorite mesosome","Stromatic leucosome","Slagstad pegmatite","Diorite mesosome", "Patch leucosome"]
colourchoice=2
include("PlotDefaults.jl")

# @recipe(TernaryText, x, y, z) do scene
#     Attributes(
#         text = "", 
#         marker = :circle,
#         markersize = 8,
#     )
# end

function solarbrowntern(filename,savedir)
    fig = Figure();
    ax = Axis(fig[1, 1], aspect = 1.0);
    
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
    
    gdf = DataFrame(CSV.File(filename))
    
    gdf[!,:norm] = gdf[!,:Na2O]*2 .+ gdf[!,:CaO] .+ gdf[!,:MgO] .+ gdf[!,:FeO] .+ gdf[!,:TiO2].+ gdf[!,:K2O]*2 
    gdf[!,:NC] = (gdf[!,:Na2O]*2 .+ gdf[!,:CaO])./gdf[!,:norm]
    gdf[!,:FMT] = (gdf[!,:MgO] .+ gdf[!,:FeO] .+ gdf[!,:TiO2])./gdf[!,:norm]
    gdf[!,:K] = (gdf[!,:K2O]*2)./gdf[!,:norm]

    
    for i in 1:lastindex(SLAGSTAD_CATS)
        cat = SLAGSTAD_CATS[i]
        
        gdfiltered = filter(:Category => c -> c == cat, gdf)
        ternaryscatter!(ax,gdfiltered[!,:NC],gdfiltered[!,:K],gdfiltered[!,:FMT], label = cat, color = myColours2[i],markersize = 6,marker=:rect,strokewidth = 0.5)
    end

    for i in 1:lastindex(MY_CATS)
        cat = MY_CATS[i]
        
        gdfiltered = filter(:Category => c -> c == cat, gdf)
        ternaryscatter!(ax,gdfiltered[!,:NC],gdfiltered[!,:K],gdfiltered[!,:FMT], label = cat, color = myColours2[i],markersize = 6,marker=:circle,strokewidth = 0.5)
        ternarytext!(ax,gdfiltered[!,:NC],gdfiltered[!,:K],gdfiltered[!,:FMT],text = gdfiltered[!,:Sample],fontsize=2)
    end
    fig[1,2] = Legend(fig,ax)
    save(joinpath(savedir,"solarbrowntern_text.svg"),fig)

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

    ternarylines!(ax,[0.7,0.575],[0,0.25],[0.30,0.175],color = :black,linewidth=2)
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
    save(joinpath(savedir,"barkertern_text2.svg"),fig)
end

function transform_qapf_coords(a,p,q)
    #Assume q,a,p each <1
    #q coord stays the same
    r = 1-q
    a2 = a*r
    p2 = p*r

    return a2,p2,q

end

geta(vec) = vec[1]
getp(vec) = vec[2]
getq(vec) = vec[3
]
function qapftern(filename,savedir)
    fig = Figure();
    ax = Axis(fig[1, 1], aspect = 1.0);
    
    ternaryaxis!(
        ax; 
        labelx = "A",
        labely = "P",
        labelz = "Q",
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
    
    gdf[!,:fsp] = gdf[!,:albite] .+ gdf[!,:anorthite] .+ gdf[!,:orthoclase] 
    gdf[!,:qnorm] = gdf[!,:quartz].+ gdf[!,:fsp]
    gdf[!,:q] = gdf[!,:quartz]./gdf[!,:qnorm]
    gdf[!,:a] = gdf[!,:orthoclase]./gdf[!,:fsp]
    gdf[!,:p] = (gdf[!,:albite] .+ gdf[!,:anorthite])./gdf[!,:fsp]

    coords = transform_qapf_coords.([0.9,0.9],[0.1,0.1],[0,0.6])
    ternarylines!(ax,geta.(coords),getp.(coords),getq.(coords),color = :black,linewidth=1)

    coords= transform_qapf_coords.([0.65,0.65],[0.35,0.35],[0,0.2])
    ternarylines!(ax,geta.(coords),getp.(coords),getq.(coords),color = :black,linewidth=1)

    coords = transform_qapf_coords.([0.35,0.35],[0.65,0.65],[0,0.6])
    ternarylines!(ax,geta.(coords),getp.(coords),getq.(coords),color = :black,linewidth=1)

    coords = transform_qapf_coords.([0.1,0.1],[0.9,0.9],[0,0.6])
    ternarylines!(ax,geta.(coords),getp.(coords),getq.(coords),color = :black,linewidth=1)

    coords = transform_qapf_coords.([1,0],[0,1],[0.05,0.05])
    ternarylines!(ax,geta.(coords),getp.(coords),getq.(coords),color = :black,linewidth=1)

    coords = transform_qapf_coords.([1,0],[0,1],[0.2,0.2])
    ternarylines!(ax,geta.(coords),getp.(coords),getq.(coords),color = :black,linewidth=1)

    coords = transform_qapf_coords.([1,0],[0,1],[0.6,0.6])
    ternarylines!(ax,geta.(coords),getp.(coords),getq.(coords),color = :black,linewidth=1)

    coords = transform_qapf_coords.([1,0],[0,1],[0.9,0.9])
    ternarylines!(ax,geta.(coords),getp.(coords),getq.(coords),color = :black,linewidth=1)

    for i in 1:lastindex(SLAGSTAD_CATS)
        cat = SLAGSTAD_CATS[i]
        
        gdfiltered = filter(:Category => c -> c == cat, gdf)
        coords = transform_qapf_coords.(gdfiltered[!,:a],gdfiltered[!,:p],gdfiltered[!,:q])
        ternaryscatter!(ax,geta.(coords),getp.(coords),getq.(coords), label = cat, color = myColours2[i],markersize = 6,marker=:rect,strokewidth = 0.5)
    end
    for i in 1:lastindex(MY_CATS)
        cat = MY_CATS[i]
        
        gdfiltered = filter(:Category => c -> c == cat, gdf)
        coords = transform_qapf_coords.(gdfiltered[!,:a],gdfiltered[!,:p],gdfiltered[!,:q])
        ternaryscatter!(ax,geta.(coords),getp.(coords),getq.(coords), label = cat, color = myColours2[i],markersize = 6,marker=:circle,strokewidth = 0.5)
        ternarytext!(ax,geta.(coords),getp.(coords),getq.(coords),text = gdfiltered[!,:Sample],fontsize=2)
    end
    fig[1,2] = Legend(fig,ax)
    save(joinpath(savedir,"qapftern_text.svg"),fig)
end
solarbrowntern("../../Geochem/FluidFlux_molpercent.csv","TernPlots/")
barkertern("../../Geochem/FluidFlux_CIPWnorm2.csv","TernPlots/")
qapftern("../../Geochem/FluidFlux_CIPWnorm2.csv","TernPlots/")
