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

sourcelib = init_meemum("23SD20A_melt-test1/MeltSource")
source_compo = getcompo(sourcelib)
close_meemum!(sourcelib)

h2ostart = 1.0
h2oend = 50.0

source_compo1 = change_list_component(source_compo,h2ostart,"H2O")
source_compo2 = change_list_component(source_compo,h2oend,"H2O")

comporange = range(source_compo1,source_compo2,10)
xh2o_range = molfrac.(comporange,"H2O")
wh2o_range = massfrac.(comporange,"H2O").*100

sources, melts, hosts = equilibrate_open_system("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Host",875,10000,800,9000,source_compo1,source_compo2, steps =10, phasefunc = [feldspar_conditions])

fig = Figure(size = (600,450))
ax = Axis(fig[1,1])
phasemode!(ax,wh2o_range,hosts)
fig[1,2] = Legend(fig,ax)
save("23SD20A_melt-test1/Host.svg",fig)

fig = Figure(size = (600,450))
ax = Axis(fig[1,1])
phasemode!(ax,wh2o_range,melts)
fig[1,2] = Legend(fig,ax)
save("23SD20A_melt-test1/Melt.svg",fig)

fig = Figure(size = (600,450))
ax = Axis(fig[1,1])
phasemode!(ax,wh2o_range,sources)
fig[1,2] = Legend(fig,ax)
save("23SD20A_melt-test1/Source.svg",fig)

# sourcelib = init_meemum("23SD20A_melt-test2/MeltSource")
# source_compo = getcompo(sourcelib)
# close_meemum!(sourcelib)

# h2ostart = 1.0
# h2oend = 50.0

# source_compo1 = change_list_component(source_compo,h2ostart,"H2O")
# source_compo2 = change_list_component(source_compo,h2oend,"H2O")
# comporange = range(source_compo1,source_compo2,10)
# xh2o_range = molfrac.(comporange,"H2O")
# wh2o_range = massfrac.(comporange,"H2O").*100
# sources, melts, hosts = equilibrate_open_system("23SD20A_melt-test2/MeltSource","23SD20A_melt-test2/Host",900,10000,825,9000,source_compo1,source_compo2, steps =10, phasefunc = [feldspar_conditions])

# fig = Figure(size = (600,450))
# ax = Axis(fig[1,1])
# phasemode!(ax,wh2o_range,hosts)
# fig[1,2] = Legend(fig,ax)
# save("23SD20A_melt-test2/Host.svg",fig)

# fig = Figure(size = (600,450))
# ax = Axis(fig[1,1])
# phasemode!(ax,wh2o_range,melts)
# fig[1,2] = Legend(fig,ax)
# save("23SD20A_melt-test2/Melt.svg",fig)

# fig = Figure(size = (600,450))
# ax = Axis(fig[1,1])
# phasemode!(ax,wh2o_range,sources)
# fig[1,2] = Legend(fig,ax)
# save("23SD20A_melt-test2/Source.svg",fig)