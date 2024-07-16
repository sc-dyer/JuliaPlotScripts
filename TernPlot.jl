using PlotlyJS

function make_ax(title, tickangle)
    attr(title=title, titlefont_size=20, tickangle=tickangle,
        tickfont_size=15, tickcolor="rgba(0, 0, 0, 0)", ticklen=2,
        showline=true, showgrid=true)
end

function plotTernCompo(plotlabels::Array{String}, compos::Matrix{Float64},elemlabels::Array{String})

    newCompos = compos
    for col in eachcol(newCompos)
        colsum = sum(col)
        for j in 1:lastindex(col)
            col[j] = col[j]/colsum*100
        end
    end
    
    t = scatterternary(mode = "markers", a = newCompos[1,:],b=newCompos[2,:],c=newCompos[3,:],text=elemlabels, marker_color = "red",marker_symbol = "circle",marker_size = 20)
    layout = Layout(ternary=attr(
        sum=100,
        aaxis=make_ax(plotlabels[1], 0),
        baxis=make_ax(plotlabels[2], 45),
        caxis=make_ax(plotlabels[3], -45),
    ))

    Plot(t,layout)
end

Components = ["SiO2","CaO","MgO"]
qz= [1.0,0.0,0.0]
per = [0.0,0.0,1.0]
en = [1.0,0.0,1.0]
di = [2.0,1.0,1.0]
ak = [2.0,2.0,1.0]
wo = [1.0,1.0,0.0]
fo = [1.0,0.0,2.0]
mer = [2.0,3.0,1.0]

minNames = ["Qz","Per","En","Di","Ak","Wo","Fo","Mer"]

compos = [qz per en di ak wo fo mer]

plotTernCompo(Components,compos,minNames)