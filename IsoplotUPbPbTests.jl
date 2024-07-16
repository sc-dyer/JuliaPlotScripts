# using CairoMakie
using DataFrames
using CSV
using Isoplot
zrn20SD06c = DataFrame(CSV.File("ZrnLaserPlots/UPb/20SD06_cores.csv"))
zrn20SD06r = DataFrame(CSV.File("ZrnLaserPlots/UPb/20SD06_rims.csv"))

zrn20SD17Ac = DataFrame(CSV.File("ZrnLaserPlots/UPb/20SD17A_cores.csv"))
zrn20SD17Ar = DataFrame(CSV.File("ZrnLaserPlots/UPb/20SD17A_rims.csv"))

zrn21SD68c = DataFrame(CSV.File("ZrnLaserPlots/UPb/21SD68_cores.csv"))
zrn21SD68r = DataFrame(CSV.File("ZrnLaserPlots/UPb/21SD68_rims.csv"))

zrn22SD55Ec = DataFrame(CSV.File("ZrnLaserPlots/UPb/22SD55E_cores.csv"))
zrn22SD55Er = DataFrame(CSV.File("ZrnLaserPlots/UPb/22SD55E_rims.csv"))

# fig = Figure(;size = (800,800))
rimdf = [zrn20SD06r zrn20SD17Ar; zrn21SD68r zrn22SD55Er]
coredf = [zrn20SD06c zrn20SD17Ac; zrn21SD68c zrn22SD55Ec]

zrndf = [zrn20SD06r, zrn20SD17Ar, zrn21SD68r, zrn22SD55Er, zrn20SD06c, zrn20SD17Ac, zrn21SD68c, zrn22SD55Ec]
labels = ["20SD06" "20SD17A";"21SSD68" "22SD55E"]

my_analyses =Array{Isoplot.Analysis}([])

for i = 1:lastindex(zrndf)
    println(labels[(i+1)÷2])
    if mod(i,2) == 1
        println("Rim")
    else
        println("Core")
    end
    df = zrndf[i]
    analysis = UPbPbAnalysis.(df[!,:Pb206U238],df[!,:Pb206U238_se1],df[!,:Pb207Pb206],df[!,:Pb207Pb206_se1],df[!,:cor_0706_0638])
    @show age76.(analysis,900.0,1900.0)
    @show df[!,:PbPbAge]
end


