using GLMakie
using DataFrames
using CSV

include("PlotDefaults.jl")
const FIG_SIZE = (1200,900)
function import_brc(filename)

    df = DataFrame(CSV.File(filename))
    df = permutedims(df,1)
    rename!(df, :Sample => :Element)
    return df
end


function plot_isocon(brcdf,unaltered,altered)
    fig = Figure(size = FIG_SIZE); display(GLMakie.Screen(),fig)
    ax = Axis(fig[1,1])
    lines!(ax,[0,maximum([brcdf[!,unaltered];brcdf[!,altered]])],[0,maximum([brcdf[!,unaltered];brcdf[!,altered]])])
    scatter!(ax,brcdf[!,unaltered],brcdf[!,altered])
    foreach(i -> text!(ax,brcdf[!,unaltered][i],brcdf[!,altered][i],text=brcdf[!,:Element][i],color=:black,fontsize=12,offset = (2,2)),1:lastindex(brcdf[!,:Element]))
    return fig
end
GLMakie.activate!()
set_theme!(myTheme)

mydf = import_brc("BRC.csv")
plot_isocon(mydf,Symbol("23SD20H"),Symbol("23SD20F"))