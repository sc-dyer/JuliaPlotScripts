using Migmatites, JPerpleX, CairoMakie, JLD2
colourchoice = 2
COLOR_REF = Dict("Cpx" => 1,
                "Melt" => 2,
                "Amph" =>3,
                "q" => 4,
                "Afs" => 5,
                "Pl" => 6,
                "Bio" => 7,
                "Ilm" => 14,
                "Mica" => 10,
                "Ep" => 13,
                "sph" => 11,
                "H2O" => 12,
                "Gt" => 8,
                "Opx" => 9
                 )
include("../PlotDefaults.jl")
set_theme!(myTheme)

function feldspar_conditions(phase)
    if lowercase(phase.name) == "fsp"
        k2o = getchemical(phase.composition,"K2O")
        na2o = getchemical(phase.composition,"Na2O")
        cao = getchemical(phase.composition,"CaO")
        if mol(k2o*2) >0.1
            return changename(phase,"Afs")
        else
            return changename(phase,"Pl")
        end
    elseif lowercase(phase.name) == "ab"
        return changename(phase,"Pl")
    else
        return phase
    end
end

function plot_sourceh2o_range(dir,source_T, host_T)

    hostlib = init_meemum(joinpath(dir,"Host_CS"))
    host_compo = hostlib.composition
    close_meemum!(hostlib)


    sourcelib = init_meemum(joinpath(dir,"MeltSource"))
    source_compo = getcompo(sourcelib)
    
    
    h2ostart = 1.0
    h2oend = 50.0
    
    source_compo1 = change_list_component(source_compo,h2ostart,"H2O")
    source_compo2 = change_list_component(source_compo,h2oend,"H2O")
    
    comporange = range(source_compo1,source_compo2,10)
    xh2o_range = molfrac.(comporange,"H2O")
    wh2o_range = massfrac.(comporange,"H2O").*100
    println("Closed system: Melt Frac = 2")
    systems = equilibrate_closed_system.((sourcelib,),comporange,(host_compo,),source_T,10000,host_T,9000,2.0,phasefunc = [feldspar_conditions])
    sources = [sys[1] for sys in systems]
    melts = [sys[2] for sys in systems] 
    hosts = [sys[3] for sys in systems]

    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    modebox!(ax,wh2o_range,hosts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Host_CS_2.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    modebox!(ax,wh2o_range,melts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Melt_CS_2.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    modebox!(ax,wh2o_range,sources)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Source_CS_2.svg"),fig)

    println("Closed system: Melt Frac = 10")
    systems = equilibrate_closed_system.((sourcelib,),comporange,(host_compo,),source_T,10000,host_T,9000,10.0,phasefunc = [feldspar_conditions])
    sources = [sys[1] for sys in systems]
    melts = [sys[2] for sys in systems] 
    hosts = [sys[3] for sys in systems]

    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    modebox!(ax,wh2o_range,hosts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Host_CS_10.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    modebox!(ax,wh2o_range,melts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Melt_CS_10.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    modebox!(ax,wh2o_range,sources)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Source_CS_10.svg"),fig)

    close_meemum!(sourcelib)

    println("Open system")
    sources, melts, hosts = equilibrate_open_system(joinpath(dir,"MeltSource"),joinpath(dir,"Host"),source_T,10000,host_T,9000,source_compo1,source_compo2, steps =10, phasefunc = [feldspar_conditions])
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    modebox!(ax,wh2o_range,hosts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Host_OS.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    modebox!(ax,wh2o_range,melts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Melt_OS.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    modebox!(ax,wh2o_range,sources)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Source_OS.svg"),fig)
end


function xh2o_vs_T(dir,xrange,tspan,geotherm, xnodes, tnodes, host_P)
    hostlib = init_meemum(joinpath(dir,"Host_CS"))
    host_compo = hostlib.composition
    close_meemum!(hostlib)


    sourcelib = init_meemum(joinpath(dir,"MeltSource"))
    source_compo = getcompo(sourcelib)
    
    
    h2ostart = xrange[1]
    h2oend = xrange[2]

    tstart = tspan[1]
    tend = tspan[2]
    
    source_compo1 = change_list_component(source_compo,h2ostart,"H2O")
    source_compo2 = change_list_component(source_compo,h2oend,"H2O")
    
    comporange = range(source_compo1,source_compo2,xnodes)
    xh2o_range = molfrac.(comporange,"H2O")
    wh2o_range = massfrac.(comporange,"H2O").*100
    close_meemum!(sourcelib)

    trange = range(tstart,tend,tnodes)
    pstart = host_P
    pend = pstart + (tend-tstart)/(geotherm/1000)
    prange = range(pstart,pend,tnodes)
    sources = []
    melts = []
    hosts = []
    
    for i in 1:lastindex(trange)
        sourcex, meltx, hostx = equilibrate_open_system(joinpath(dir,"MeltSource"),joinpath(dir,"Host"),trange[i],prange[i],tstart,host_P,source_compo1,source_compo2, steps =xnodes, phasefunc = [feldspar_conditions])
        fig = Figure(size = (500,1500))
        ax1 = Axis(fig[1,1],aspect = 1.0)
       
        ax2 = Axis(fig[2,1],aspect=1.0)
        
        ax3 = Axis(fig[3,1],aspect=1.0)
        
        modebox!(ax1,wh2o_range,sourcex)
        fig[1,2] = Legend(fig,ax1)
        modebox!(ax2,wh2o_range,meltx)
        fig[2,2] = Legend(fig,ax2)
        modebox!(ax3,wh2o_range,hostx)
        fig[3,2] = Legend(fig,ax3)
        colsize!(fig.layout,1,Fixed(450))
        save(joinpath(dir,"modeboxes_"*string(trange[i])*".svg"),fig)
        if length(sources) == 0
            sources = sourcex
            melts = meltx
            hosts = hostx
        else
            sources = [sources sourcex]
            melts = [melts meltx]
            hosts = [hosts hostx]
        end
    end
    
    
    # hosts = transpose(hosts)
    # melts = transpose(melts)
    # sources = transpose(sources)
    fig = Figure()
    ax1 = Axis(fig[1,1],aspect=1.0)
    
    meltprops = get_volprop.(hosts,"Melt")
    amphprops = get_volprop.(hosts,"Amph","Melt")
    
    # @show typeof(meltprops)
    # @show meltprops
    c1 = contourf!(ax1,wh2o_range,trange,meltprops)
    Colorbar(fig[1,2],c1)
    # contour!(ax2,wh2o_range,trange,amphprops,labels=true)
    # 
    save(joinpath(dir,"host_melt_contours.svg"),fig)
    fig = Figure()
    ax1 = Axis(fig[1,1],aspect=1.0)
    ax2 = Axis(fig[2,1],aspect=1.0)
    melth2o = massfrac.(getcompo.(melts),"H2O")
   
    hosth2o = massfrac.(getcompo.(hosts),"H2O")
    
    c1 = contourf!(ax1,wh2o_range,trange,melth2o)
    Colorbar(fig[1,2],c1)
    c2 = contourf!(ax2,wh2o_range,trange,hosth2o)
    Colorbar(fig[2,2],c2)
    save(joinpath(dir,"h2ocontent.svg"),fig)
    fig = Figure()
    ax1 = Axis(fig[1,1],aspect=1.0)
    meltprops = get_volprop.(sources,"Melt")


    c1 = contourf!(ax1,wh2o_range,trange,meltprops)
    Colorbar(fig[1,2],c1)
    save(joinpath(dir,"source_melt_contours.svg"),fig)
end

function melt_path_xh2o(dir,xrange,host_T, host_P, source_P, source_t_start, source_t_stop;xnodes = 10,tsteps = 500)
    hostfile = joinpath(dir,"Host")
    sourcefile = joinpath(dir,"MeltSource")
    t_path = [source_t_start,source_t_stop]
    p_path = [source_P,source_P]
    
    sourcelib = init_meemum(joinpath(dir,"MeltSource"))
    source_compo = getcompo(sourcelib)
    
    
    h2ostart = xrange[1]
    h2oend = xrange[2]
    
    source_compo1 = change_list_component(source_compo,h2ostart,"H2O")
    source_compo2 = change_list_component(source_compo,h2oend,"H2O")
    
    comporange = range(source_compo1,source_compo2,xnodes)
    xh2o_range = molfrac.(comporange,"H2O")
    wh2o_range = massfrac.(comporange,"H2O").*100
    close_meemum!(sourcelib)
    t_range = range(source_t_start,source_t_stop,tsteps)
    host_list = []
    t_list =[]
    xaxis_var = []
    count = 1
    for compo in comporange
    
        melts_f, hosts, restite_compos, system_steps, extract_T, extract_P = 
                                    melt_cycle_model_open(hostfile,sourcefile,host_T,host_P,t_path,p_path,numsteps = tsteps,source_compo = compo,phasefunc = [feldspar_conditions])
        
       
        fig = Figure()
        ax1 = Axis(fig[1,1],aspect = 1.0, title = "Melt Source")
        ax1.xlabel = "T (°C)"
        ax1.ylabel = "Modal abundance"
        limits!(ax1,source_t_start,source_t_stop,0.0,1.0)
        modebox!(ax1,t_range,system_steps,colorreference = COLOR_REF,colormap = :Paired_12)
        fig[1,2] = Legend(fig,ax1)
        for t in extract_T
            lines!(ax1,[t,t],[0,1], color=:black, linestyle = :dash)
        end

        if length(hosts) > 0
            ax2 = Axis(fig[2,1], title = "Host",yticklabelcolor = colourscheme[1])
            ax3 = Axis(fig[2,1],yticklabelcolor = colourscheme[2],yaxisposition = :right)
            ax2.xlabel = "Source melt extraction T (°C)"
            ax2.ylabel = "Host melt \nvol proportion"
            ax3.ylabel = "H2O in host \n(weight %)"
            hidespines!(ax3)
            hidexdecorations!(ax3)
            limits!(ax2,source_t_start,source_t_stop,0.0,1.0)
            limits!(ax3,source_t_start,source_t_stop,0.0,20.0)
            extract_index = 1
            while extract_T[extract_index] < host_T && extract_index < lastindex(extract_T)
                extract_index += 1
            end
            push!(t_list,extract_T[extract_index])
            push!(host_list,hosts[extract_index])
            push!(xaxis_var, massfrac(compo,"H2O")*100)
            modalmelt = get_volprop.(hosts,"melt")
            h2oInHost = massfrac.(getcompo.(hosts),"H2O")*100

            scatterlines!(ax2,extract_T,modalmelt,color=colourscheme[1])
            scatterlines!(ax3,extract_T,h2oInHost,color = colourscheme[2])
            # if length(hosts)==1
            #     modebox!(ax2,[extract_T[1],extract_T[1]],[hosts[1],hosts[1]],colorreference = COLOR_REF)
            # else
            #     modebox!(ax2,extract_T,hosts,colorreference = COLOR_REF)
            # end
            
            rowgap!(fig.layout,0)
            colsize!(fig.layout,1,Aspect(1,1.0))
            colsize!(fig.layout,2, Relative(1/5))
            rowsize!(fig.layout, 2, Relative(1/3))
            
        end
        jldname = joinpath(dir,"outputs_x"*string(count)*".jld2")
        jldsave(jldname;melts_f,hosts,restite_compos,system_steps,extract_T,extract_P,host_T,host_P,t_range)
        save(joinpath(dir,"modeboxes_"*string(massfrac(compo,"H2O")*100)*".svg"),fig)
        count +=1
    end

    fig = Figure(size = (500,500))
    ax = Axis(fig[1,1],aspect = 1.0,title = "Host (T = "*string(host_T)*"°C)")
    modebox!(ax,xaxis_var,host_list,colorreference = COLOR_REF,colormap = :Paired_12)
    fig[1,2] = Legend(fig,ax)
    limits!(ax,xrange[1],xrange[2],0,1.0)
    ax.xlabel = "H2O in source (mol %)"
    ax.ylabel = "Modal abundance"


    ax2 = Axis(fig[2,1],aspect=2.0, title = "Melt Source")
    lines!(ax2,xaxis_var,t_list)
    limits!(ax2,xrange[1],xrange[2],700,1000)
    ax2.xlabel = "H2O in source (mol %)"
    ax2.ylabel =  "Source melt extraction T (°C)"
    save(joinpath(dir,"modebox_summary.svg"),fig)

end

function plot_saves(maindir,sourcename,host1,host2,pindex)
    fig = Figure(size = (2000,800))
    for i in 1:lastindex(pindex)
        
        data750host1 = load(joinpath(maindir,sourcename*"_source_"*host1*"_host_750/","outputs_x"*string(pindex[i])*".jld2"))
        data800host1 = load(joinpath(maindir,sourcename*"_source_"*host1*"_host_800/","outputs_x"*string(pindex[i])*".jld2"))
        data750host2 = load(joinpath(maindir,sourcename*"_source_"*host2*"_host_750/","outputs_x"*string(pindex[i])*".jld2"))
        data800host2 = load(joinpath(maindir,sourcename*"_source_"*host2*"_host_800/","outputs_x"*string(pindex[i])*".jld2"))
        system_steps = data750host1["system_steps"]
        extract_T = data750host1["extract_T"]
        wh2o = round(massfrac(getcompo(system_steps[1]),"H2O")*100,sigdigits = 2)
        sourcetitle = "Melt source ("*string(wh2o)*" wt% H2O)"
        ax1 = Axis(fig[1,i],aspect = 1.0, title = sourcetitle)
        ax1.xlabel = "T (°C)"
        if i == 1
            ax1.ylabel = "Modal abundance"
        else
            hideydecorations!(ax1,ticks=false,minorticks=false)
        end
        limits!(ax1,600,1000,0.0,1.0)

        t_range = range(600,1000,500)
        modebox!(ax1,t_range,system_steps,colorreference = COLOR_REF,colormap = :Paired_12)
        
        for t in extract_T
            lines!(ax1,[t,t],[0,1], color=:black, linestyle = :dash)
        end

        ax2a = Axis(fig[2,i], title = "Host (750°C)",yticklabelcolor = colourscheme[1])
        ax2b = Axis(fig[2,i],yticklabelcolor = colourscheme[2],yaxisposition = :right)
        # ax2a.xlabel = "Source melt extraction T (°C)"
        if i == 1
            ax2a.ylabel = "Host melt \nvol proportion"

        else
            hideydecorations!(ax2a,ticks=false,minorticks=false)
        end
        if i == lastindex(pindex)
            ax2b.ylabel = "H2O in host \n(weight %)"
        else
            hideydecorations!(ax2b,ticks=false,minorticks=false)
        end
        limits!(ax2a,600,1000,0.0,1.0)
        limits!(ax2b,600,1000,0.0,20.0)
        ax2b.yticks= [0,10,20]
        hidespines!(ax2b)
        hidexdecorations!(ax2a,ticks=false,minorticks=false)
        hidexdecorations!(ax2b)
        hosts = data750host1["hosts"]
        modalmelt = get_volprop.(hosts,"melt")
        h2oInHost = massfrac.(getcompo.(hosts),"H2O")*100

        scatterlines!(ax2a,extract_T,modalmelt,color=colourscheme[1])
        scatterlines!(ax2b,extract_T,h2oInHost,color = colourscheme[2])

        hosts = data750host2["hosts"]
        modalmelt = get_volprop.(hosts,"melt")
        h2oInHost = massfrac.(getcompo.(hosts),"H2O")*100

        scatterlines!(ax2a,extract_T,modalmelt,color=colourscheme[1],marker=:rect,linestyle=:dash)
        scatterlines!(ax2b,extract_T,h2oInHost,color = colourscheme[2],marker=:rect,linestyle=:dash)

        ax3a = Axis(fig[3,i], title = "Host (800°C)",yticklabelcolor = colourscheme[1])
        ax3b = Axis(fig[3,i],yticklabelcolor = colourscheme[2],yaxisposition = :right)
        ax3a.xlabel = "Source melt extraction T (°C)"
        if i == 1
            ax3a.ylabel = "Host melt \nvol proportion"
        else
            hideydecorations!(ax3a,ticks=false,minorticks=false)
        end
        if i == lastindex(pindex)
            ax3b.ylabel = "H2O in host \n(weight %)"
        else
            hideydecorations!(ax3b,ticks=false,minorticks=false)
        end
        limits!(ax3a,600,1000,0.0,1.0)
        limits!(ax3b,600,1000,0.0,20.0)
        hidespines!(ax3b)
        hidexdecorations!(ax3b)
        hosts = data800host1["hosts"]
        modalmelt = get_volprop.(hosts,"melt")
        h2oInHost = massfrac.(getcompo.(hosts),"H2O")*100

        scatterlines!(ax3a,extract_T,modalmelt,color=colourscheme[1])
        scatterlines!(ax3b,extract_T,h2oInHost,color = colourscheme[2])
        hosts = data800host2["hosts"]
        modalmelt = get_volprop.(hosts,"melt")
        h2oInHost = massfrac.(getcompo.(hosts),"H2O")*100

        scatterlines!(ax3a,extract_T,modalmelt,color=colourscheme[1],marker=:rect,linestyle=:dash)
        scatterlines!(ax3b,extract_T,h2oInHost,color = colourscheme[2],marker=:rect,linestyle=:dash)
        # if length(hosts)==1
        #     modebox!(ax2,[extract_T[1],extract_T[1]],[hosts[1],hosts[1]],colorreference = COLOR_REF)
        # else
        #     modebox!(ax2,extract_T,hosts,colorreference = COLOR_REF)
        # end
        
        rowgap!(fig.layout,0)
        colsize!(fig.layout,i,Aspect(1,1.0))
        rowsize!(fig.layout, 2, Relative(1/6))
        rowsize!(fig.layout, 3, Relative(1/6))
        if i == lastindex(pindex)
            fig[1,i+1] = Legend(fig,ax1)
        end
    end

    save(joinpath(maindir,sourcename*"_source_modeboxes.svg"),fig)
end
# println("Group 1")
# plot_sourceh2o_range("23SD20A_source_23SD20A_host1/",875,800)
# println("Group 2")
# plot_sourceh2o_range("23SD20A_source_23SD20A_host2/",900,825)
# println("Group 3")
# plot_sourceh2o_range("Mesosome_source_23SD20A_host1/",875,800)
# println("Group 4")
# plot_sourceh2o_range("Mesosome_source_23SD20A_host2/",900,825)
# println("Group 5")
# plot_sourceh2o_range("Mesosome_source_Mesosome_host1/",875,800)
# println("Group 6")
# plot_sourceh2o_range("Mesosome_source_Mesosome_host2/",900,825)

# xh2o_vs_T("23SD20A_source_23SD20A_host1/",(1,20),(875,975),100,20,20,8900)
# xh2o_vs_T("Mesosome_source_23SD20A_host1/",(1,20),(875,975),100,20,20,8900)
# xh2o_vs_T("Mesosome_source_Mesosome_host1/",(1,20),(875,975),100,20,20,8900)
# melt_path_xh2o("MeltPathPlots/23SD03B_source_23SD03B_host_750/",[0.001,10.0],750,9000,10000,600,1000)
# melt_path_xh2o("MeltPathPlots/23SD03B_source_23SD20A_host_750/",[0.001,10.0],750,9000,10000,600,1000)
# melt_path_xh2o("MeltPathPlots/23SD20A_source_23SD20A_host_750/",[0.001,10.0],750,9000,10000,600,1000)

# melt_path_xh2o("MeltPathPlots/23SD03B_source_23SD03B_host_800/",[0.001,10.0],800,9000,10000,600,1000)
# melt_path_xh2o("MeltPathPlots/23SD03B_source_23SD20A_host_800/",[0.001,10.0],800,9000,10000,600,1000)
# melt_path_xh2o("MeltPathPlots/23SD20A_source_23SD20A_host_800/",[0.001,10.0],800,9000,10000,600,1000)

plot_saves("MeltPathPlots/","23SD03B","23SD03B","23SD20A",[2,3,5,10])
# plot_saves("MeltPathPlots/","23SD20A","23SD03B",[2,3,5,10])
# plot_saves("MeltPathPlots/","23SD20A","23SD20A",[2,3,5,10])