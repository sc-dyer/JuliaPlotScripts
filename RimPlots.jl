using GLMakie
using DataFrames
using CSV

function importfiles(filedir)
    te_df = DataFrame(CSV.file(filedir*"TraceElementData.csv"))
    uth_df = DataFrame(CSV.file(filedir*"UPbData.csv"))
    uth_df[!,:Distance] = uth_df[!,:DistanceMod]
    selections = DataFrame(CSV.file(filedir*"TraceElementsAvgs.csv"))
    rim_te = DataFrame[]
    rim_uth = DataFrame[]
    for row in eachrow(selections)
        push!(rim_te,select_range(te_df,row[:x1],row[:x2]))
        push!(rim_uth,select_range(uth_df,row[:x1],row[:x2]))
    end

    return rim_te, rim_uth
end

function select_range(df,x1,x2)
    
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

end

function plot_rims()

end