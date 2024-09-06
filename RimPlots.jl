using GLMakie
using DataFrames
using CSV

mutable struct RimData
    name::String
    tes::Array{DataFrame}
    uth::Array{DataFrame}
end
function importfiles(filedir)
    te_df = DataFrame(CSV.File(joinpath(filedir,"TraceElementData.csv")))
    uth_df = DataFrame(CSV.File(joinpath(filedir,"UPbData.csv")))
    uth_df[!,:Distance] = uth_df[!,:DistanceMod]
    selections = DataFrame(CSV.File(joinpath(filedir,"TraceElementsAvgs.csv")))
    rim_te = DataFrame[]
    rim_uth = DataFrame[]
    for row in eachrow(selections)
        push!(rim_te,select_range(te_df,row[:x1],row[:x2]))
        push!(rim_uth,select_range(uth_df,row[:x1],row[:x2]))
    end

    return rim_te, rim_uth
end

function select_range(df,x1,x2)
    rimdf = df
    if x1 < x2
        rimdf = filter(:Distance => x -> x1 <= x <= x2, df)
    
    elseif x2 > x1
        rimdf = filter(:Distance => x -> x2 <= x <= x1, df)
        d_array = copy(rimdf[!,:Distance])
        reverse!(rimdf)
        rimdf[!,:Distance] = d_array
    end

    #set first to x = 0
    xmin = minimum(rimdf[!,:Distance])
    rimdf[!,:Distance] = rimdf[!,:Distance] .- xmin
    arbitrary_scale!(rimdf)

    return rimdf
end

function arbitrary_scale!(df;x1 = 0.0, x2 = 100.0)
    increment = (x2 - x1)/(nrow(df)-1)

    x_arb = collect(x1:increment:x2)

    df[!,:x_arb] = x_arb

end

function import_all(directory)
    entries = readdir(directory)

    zrns = RimData[]

    for entry in entries
        if isdir(joinpath(directory,entry))
            name = entry
            tes, uth = importfiles(joinpath(directory,entry))
            push!(zrns,RimData(name,tes,uth))
        end
    end
    
    return zrns
end

function plot_rims(rims)


end

function changeElem!(df, y, elem)

    if elem == "U" || elem == "Th"
        y[] = df[!,Symbol(elem*"_rescaled")]
    else
        y[] = df[!,Regex(elem*"\\d+")][!,1]
    end

    return maximum(y)

end

function changeElem!(rims,ys, elem, ax)

    for rim in rims

    end
end

zrns = import_all("ZrnLaserPlots/20SD06/TransectScaling/")

