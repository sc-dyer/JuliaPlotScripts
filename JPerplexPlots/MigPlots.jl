using Migmatites, JPerpleX, CairoMakie
colourchoice = 2
include("../PlotDefaults.jl")

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
    
    systems = equilibrate_closed_system.(fill(sourcelib,10),comporange,fill(host_compo,10),source_T,10000,host_T,9000,2.0)
    sources = [sys[1] for sys in systems]
    melts = [sys[2] for sys in systems] 
    hosts = [sys[3] for sys in systems]

    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    phasemode!(ax,wh2o_range,hosts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Host_CS_2.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    phasemode!(ax,wh2o_range,melts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Melt_CS_2.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    phasemode!(ax,wh2o_range,sources)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Source_CS_2.svg"),fig)


    systems = equilibrate_closed_system.((sourcelib,),comporange,(host_compo,),source_T,10000,host_T,9000,10.0)
    sources = [sys[1] for sys in systems]
    melts = [sys[2] for sys in systems] 
    hosts = [sys[3] for sys in systems]

    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    phasemode!(ax,wh2o_range,hosts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Host_CS_10.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    phasemode!(ax,wh2o_range,melts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Melt_CS_10.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    phasemode!(ax,wh2o_range,sources)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Source_CS_10.svg"),fig)

    close_meemum!(sourcelib)


    sources, melts, hosts = equilibrate_open_system(joinpath(dir,"MeltSource"),joinpath(dir,"Host"),source_T,10000,host_T,9000,source_compo1,source_compo2, steps =10, phasefunc = [feldspar_conditions])
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    phasemode!(ax,wh2o_range,hosts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Host_OS.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    phasemode!(ax,wh2o_range,melts)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Melt_OS.svg"),fig)
    
    fig = Figure(size = (600,450))
    ax = Axis(fig[1,1])
    phasemode!(ax,wh2o_range,sources)
    fig[1,2] = Legend(fig,ax)
    save(joinpath(dir,"Source_OS.svg"),fig)
end

plot_sourceh2o_range("23SD20A_source_23SD20A_host1/",875,800)
plot_sourceh2o_range("23SD20A_source_23SD20A_host2/",900,825)
plot_sourceh2o_range("Mesosome_source_23SD20A_host1/",875,800)
plot_sourceh2o_range("Mesosome_source_23SD20A_host2/",900,825)
plot_sourceh2o_range("Mesosome_source_Mesosome_host1/",875,800)
plot_sourceh2o_range("Mesosome_source_Mesosome_host2/",900,825)