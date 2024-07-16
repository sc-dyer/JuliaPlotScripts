
#using AlgebraOfGraphics
using DataFrames
using CSV
using GLMakie
using Clustering
using Distances
using NearestNeighbors
# myColours = ["#FF8CA6","#6D8FFF","#F9C52A","#A9EFA5", "#E6A7FB","#000000","#FF0000","#008000","#0000FF","#4B0082","#FF8C00"]
include("PlotDefaults.jl")


function clusterXY!(zrnMeas::DataFrame, numClust::Integer)
    xyPos = Matrix(zrnMeas[!,r"(mm)"])
    distMat = pairwise(Euclidean(),transpose(xyPos),dims=2)
    
    # xyClust = hclust(distMat,:single)
    # clustIndex = cutree(xyClust,k=numClust)
    # zrnMeas[!,:Cluster] = clustIndex

    xyClust = kmedoids(distMat,numClust)
    zrnMeas[!,:Cluster] = xyClust.assignments

    # xyClust = kmeans(distMat,numClust)
    # zrnMeas[!,:Cluster] = assignments(xyClust)

end

function plotValues(zrnMeas::DataFrame, xVar::String, yVar::String, ax::Axis,fig::Figure,fileName::String;clustIndex::Int=0)
    empty!(ax)
    ax.xlabel = xVar
    ax.ylabel = yVar
    # zrnMeas = dropmissing(zrnMeas)
    if clustIndex == 0
        clustVal = 1
        iClust = filter(:Cluster=>==(clustVal),zrnMeas)
        while nrow(iClust) > 0
            scatter!(iClust[!,Symbol(xVar)],iClust[!,Symbol(yVar)], color = myColours[clustVal])
            
            clustVal +=1
            iClust = filter(:Cluster=>==(clustVal),zrnMeas)
        end
    else
        iClust = filter(:Cluster=>==(clustIndex),zrnMeas)
        scatter!(iClust[!,Symbol(xVar)],iClust[!,Symbol(yVar)], color = myColours[clustIndex])
    end
    # save("ZrnPlots//"*fileName*".svg",fig)
    
end

function nn_distance(df)
    points = transpose([df[!,"X (mm)"] df[!,"Y (mm)"]])
    kdtree = KDTree(points)
    idxs, dists = knn(kdtree, points, 2)
  
    return [dists[i][1] for i=1:lastindex(dists)]
end

function rangesearch(df,threshold_distance)
    points = transpose([df[!,"X (mm)"] df[!,"Y (mm)"]])
    kdtree = KDTree(points)
    counts = inrangecount(kdtree,points,threshold_distance) .- 1
    return counts
end

GLMakie.activate!()
set_theme!(myTheme)

zrnData =  DataFrame(CSV.File(raw"/home/scdyer/Documents/Waterloo/SEM/Zircon CL/ZrnMeasurements/20SD06.csv"))
# zrnData = dropmissing(zrnData)

fig = Figure()
ax = Axis(fig[1,1])
# hist!(zrnData[!,"a (μm)"], strokewidth = 1, strokecolor = :black, bins = 20)

zrnData[!,:nndist] = nn_distance(zrnData)
zrnData[!,:nncounts] = rangesearch(zrnData,4.0)
zrnData = dropmissing(zrnData)
# scatter!(ax,nndist, zrnData[!,"a-a0"])
# boxplot!(ax,zrnData[!,:nncounts],zrnData[!,"a-a0"])
clusterXY!(zrnData,6)
# plotValues(zrnData, "Core Area (μm^2)", "Rim Area (μm^2)",ax,fig,"20SD06_Areas1",clustIndex = 0)
# plotValues(zrnData, "Co   re Area (μm^2)", "Rim Area (μm^2)",ax,fig,"20SD06_Areas2",clustIndex = 2)
# plotValues(zrnData, "Core Area (μm^2)", "Rim Area (μm^2)",ax,fig,"20SD06_Areas3",clustIndex = 3)
# plotValues(zrnData, "Core Area (μm^2)", "Rim Area (μm^2)",ax,fig,"20SD06_Areas4",clustIndex = 4)
# plotValues(zrnData, "Core Area (μm^2)", "Rim Area (μm^2)",ax,fig,"20SD06_Areas5",clustIndex = 5)
# plotValues(zrnData, "Core Area (μm^2)", "Rim Area (μm^2)",ax,fig,"20SD06_Areas6",clustIndex = 6)

plotValues(zrnData, "a0 (μm)", "a-a0",ax,fig,"20SD06_a1",clustIndex = 4)
# plotValues(zrnData, "a0 (μm)", "a-a0",ax,fig,"20SD06_a2",clustIndex = 2)
# plotValues(zrnData, "a0 (μm)", "a-a0",ax,fig,"20SD06_a3",clustIndex = 3)
# plotValues(zrnData, "a0 (μm)", "a-a0",ax,fig,"20SD06_a4",clustIndex = 4)
# plotValues(zrnData, "a0 (μm)", "a-a0",ax,fig,"20SD06_a5",clustIndex = 5)
# plotValues(zrnData, "a0 (μm)", "a-a0",ax,fig,"20SD06_a6",clustIndex = 6)

# # plotValues(zrnData, "a0 (μm)", "a-a0",ax,fig,"20SD06_Areas")
# plotValues(zrnData, "X (mm)", "Y (mm)",ax,fig,"20SD06_XY")

fig