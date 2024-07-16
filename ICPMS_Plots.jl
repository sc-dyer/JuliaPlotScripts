#Some simple functions for plotting ICPMS data processed by Iolite
using Plots
using Measurements
using CSV
using DataFrames
using ColorSchemes
using LaTeXStrings

PLOTX = 600
PLOTY = 600
MARK_SIZE = 6
LINE_WIDTH = 3
MIN_SPIDER = 0.01
MAX_SPIDER = 100
MIN_Y = :auto
MAX_Y = :auto
YSCALING = :identity
SPIDERSCALING = :log10
XSCALING = :log10

dscatter() = scatter(size = (PLOTX,PLOTY), xscale = XSCALING, yscale = YSCALING, ylims = (MIN_Y,MAX_Y), 
framestyle = :box, legend_position = :inside, legendfontsize = 12, grid = false, minorticks = true,
tickfontsize = 12, labelfontsize = 16)

CATEGORIES = [["Rim","Core","Mixed"],["Zrn","Ap","Grt","Cpx","Fsp"]]
NULL_CAT = "N/A"
CAT_REGEX = [[r"Rim",r"Core",r"Mixed"],[r"Zrn",r"Ap",r"Grt",r"Cpx",r"Fsp"]]
REES = ["La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"]

MEAS_UNIT = "μg/g"
MARK_SHAPE = [:circle,:diamond,:dtriangle]
MARK_COLOR = [colorant"#FF8CA6",colorant"#6D8FFF",colorant"#F9C52A"]


SAVEDIR = "Plots/"


#Individual ICPMS mass measurement
struct Element
    name::String
    mass::Int64 #Mass of the isotope measured
    concentration::Measurement
end

Element() = Element("None", 0, 0)

struct TopeRatio
    topeNum::Element
    topeDenom::Element
    ratio::Measurement
end

TopeRatio() = TopeRatio(Element(), Element(), 0)
struct Age
    topeNum::Element
    topeDenom::Element
    age::Measurement
end
Age() = Age(Element(), Element(), 0)
#Individual ICPMS analysis of a single sample
struct Sample
    name::String
    elements::Array{Element}
    unit::String #Assumes every element has the same unit of measure in a single sample
    categories::Array{String} #Allows arbitrary categorization of each sample
    xCoord::Float64
    yCoord::Float64
    mineralNum::Integer
    ssRatios::Array{TopeRatio}
    ssAges::Array{Age}
end



function getElem(elems::Array{Element}, el::String)
    for elem in elems
        if elem.name == el
            return elem
        end
    end
end

function getElem(sam::Sample, el::String)
    getElem(sam.elements, el)
end

#Ratio e.g. "Pb206-U238"
function getAge(ages::Array{Age}, ratio::String)
    elemSplit = split(ratio, "-")
    elem1 = match(r"[A-z]+",elemSplit[1]).match
    mass1 = parse(Int64, match(r"\d+",elemSplit[1]).match)
    elem2 = match(r"[A-z]+",elemSplit[2]).match
    mass2 = parse(Int64, match(r"\d+",elemSplit[2]).match)
    for age in ages
        if elem1 == age.topeNum.name && mass1 == age.topeNum.mass && elem2 == age.topeDenom.name && mass2 == age.topeDenom.mass
            return age
        end
    end
end

function getAge(sam::Sample, ratio::String)
    getAge(sam.ssAges, ratio)
end
function conc(elem::Element)
    return elem.concentration
end
function getName(sam::Sample)
    return sam.name
end

function getX(sam::Sample)
    return sam.xCoord
end

function getY(sam::Sample)
    return sam.yCoord
end

value(m::Measurement) = begin
    m.val
end

#Basic function to categorize samples based on names, easily modified based on what you want
function categorize(catString::String)

    cats = String[]

    #Categorize based on core/rim
    for i in 1:lastindex(CAT_REGEX)
        
        containsCat = false
        for j in 1:lastindex(CAT_REGEX[i])
            if occursin(CAT_REGEX[i][j],catString) && !containsCat
                push!(cats,CATEGORIES[i][j])
            end
        end
        if !containsCat
            push!(cats,NULL_CAT)
        end
       

    end

    return cats
    
end

function findMinNum(sampleName::String)
#Assumes format ..._Min##...
    minNumStr = match(r"\d+",match(r"_[A-z]+\d+",sampleName).match).match

    return parse(Int64, minNumStr)
end

function importNorm(fileName::String)
    norms = CSV.File(fileName)
    elems = Element[]

    for row in norms
        name = String(row.Element)
        concentration = row.Concentration
        push!(elems, Element(name,0,concentration))

    end
    return elems
end

function importPartition(fileName::String)
    partitions = CSV.File(fileName)
    slope = Measurement[]
    Yb = Measurement[]
    for row in partitions
        push!(slope,row.logDslope ± row.seSlope)
        push!(Yb, row.logDyb ± row.seYb)
    end

    return slope, Yb
end


function importData(fileName::String, ioliteVersion::Integer)

    icpmsDF = DataFrame(CSV.File(fileName))
    samples = Sample[]
    #This will iterate through each row 3 columns at a time, assuming concentration is listed followed by standard error
    #Assumed headder style for iolite3 files: Element_un(_)it_m# where there may or may not be a second "_" to replace a "/" for units like μg/g (ug_g)
    #This is getting the units
    colName = propertynames(icpmsDF)[5]
    splitColName = split(String(colName),"_")
    unit = ""
    if length(splitColName) > 3 #Check for the case of an extra _
        #Change the ug to μg if present
        if splitColName[2]=="ug"
            unit = "μg/g"
        else
            unit = splitColName[2]*"/"*splitColName[3]
        end
    else
        unit = splitColName[2]
    end
    unit = MEAS_UNIT #Comment out if you dont want to override

    
    for i = 1:lastindex(icpmsDF,1)
        name = String(icpmsDF[i,1])
        cats = String(icpmsDF[i,2])
        if ioliteVersion == 3
            x = icpmsDF[i,3]
            y = icpmsDF[i,4]
            firstDataCol = 5
        elseif ioliteVersion == 4
            xy = icpmsDF[i,3]
            xyRegex = collect(eachmatch(r"\d+",xy))
            x = parse(Int64, xyRegex[1].match)
            y = parse(Int64, xyRegex[2].match)
            firstDataCol = 4
        end
        elems = Element[]
        ssData = false
        ssIndex =  1
        for j = firstDataCol:3:lastindex(icpmsDF,2)
            #Determine the element and the mass based on column name
            colName = String(propertynames(icpmsDF)[j])
            splitColName = split(colName,"_")
            splitColNameSS = split(colName," ")
            if splitColNameSS[1] == "Final"
                ssData = true
                ssIndex = j
                break
            end
            if ioliteVersion == 3
                elemName = splitColName[1]
                mass = parse(Int64, match(r"\d+",last(splitColName)).match)
            elseif ioliteVersion ==4
                elemName = match(r"[A-z]+",splitColName[1]).match
                mass = parse(Int64, match(r"\d+",splitColName[1]).match)
            end
            #retrieve concentration and error (from the next column)
            

            concentration = icpmsDF[i,j]
            error = icpmsDF[i,j+1]
            if !(concentration isa Number)
                
                #println(name,": ",elemName,": ",concentration,"±",error)
                #println(typeof(concentration))
                try
                    concentration = parse(Float64,String(concentration))
                    if !(error isa Number)
                        error = parse(Float64,String(error))
                    end
                    #println(name,": ",elemName,": ",concentration,"±",error)
                catch err
                    concentration = NaN
                    error = NaN
                    #@error "ERROR: " exception=(err, catch_backtrace())
                end
                # println(name,": ",elemName,": ",concentration,"±",error)
                
            end
            
            thisElem = Element(elemName,mass,concentration±error)
            push!(elems,thisElem)
        end
        
        if ssData
            ratios = TopeRatio[]
            ages = Age[]
            for j = ssIndex:2:lastindex(icpmsDF,2)
                colName = String(propertynames(icpmsDF)[j])
                splitColName = split(colName," ")
                if splitColName[1] == "rho"
                    break
                end

                splitElems = split(splitColName[2],"/")
                numerElem =  match(r"[A-z]+",splitElems[1]).match
                numerMass = parse(Int64, match(r"\d+",splitElems[1]).match)
                denomElem =  match(r"[A-z]+",splitElems[2]).match
                denomMass =  parse(Int64, match(r"\d+",splitElems[2]).match)

                val = icpmsDF[i,j] #Tope ratio or age
                error = icpmsDF[i,j+1]

                if !(val isa Number)
                
                    # println(name,": ",elemName,": ",concentration,"±",error)
                    try
                        val = parse(Float64,String(concentration))
                        error = parse(Float64,String(error))
                        
                    catch
                        val = NaN
                        error = NaN
                    end
                    # println(name,": ",elemName,": ",concentration,"±",error)
                    
                end
                if length(splitColName) ==3
                    age = Age(Element(numerElem,numerMass,0), Element(denomElem, denomMass, 0),val±error)
                    push!(ages, age)

                else
                    ratio = TopeRatio(Element(numerElem,numerMass,0), Element(denomElem, denomMass, 0),val±error)
                    push!(ratios, ratio)
                end

            end
            thisSample = Sample(name,elems,unit,categorize(cats),x,y, findMinNum(name), ratios,ages)
        else
            thisSample = Sample(name,elems,unit,categorize(cats),x,y, findMinNum(name), [TopeRatio()],[Age()])
        end


        push!(samples,thisSample)

    end

    return samples    

end


function addXYVal!(x, y, sample::Sample, xElem::String , yElem::String, xElemD::String , yElemD::String,coordAx::Integer = 0)
    if coordAx > 0
        if coordAx == 1
            push!(x,sample.xCoord)
        elseif coordAx == 2
            push!(x,sample.yCoord)
        else
            throw(ErrorException("Error: invalid coordinate axis input"))
        end
        if isempty(yElemD)
            push!(y,getElem(sample,yElem).concentration)
        else
            if yElemD =="Eu*"
                push!(y,getElem(sample,yElem).concentration/sqrt(getElem(sample,"Sm")*getElem(sample,"Gd")))
            else
                push!(y,getElem(sample,yElem).concentration/getElem(sample,yElemD).concentration)
            end
        end
    elseif isempty(xElemD)

        if occursin(r"\d+", xElem)
            push!(x,getAge(sample,xElem).age)
        else
            push!(x,getElem(sample,xElem).concentration)
        end
        if isempty(yElemD)
            #println(y)
            push!(y,getElem(sample,yElem).concentration)
        else
            if yElemD =="Eu*"
                push!(y,getElem(sample,yElem).concentration/sqrt(conc(getElem(sample,"Sm"))*conc(getElem(sample,"Gd"))))
            else
                push!(y,getElem(sample,yElem).concentration/getElem(sample,yElemD).concentration)
            end
        end
    else
        if xElemD =="Eu*"
            push!(x,getElem(sample,xElem).concentration/sqrt(conc(getElem(sample,"Sm"))*conc(getElem(sample,"Gd"))))
        else
            push!(x,getElem(sample,xElem).concentration/getElem(sample,xElemD).concentration)
        end
        if isempty(yElemD)
            
            push!(y,getElem(sample,yElem).concentration)
        else
            if yElemD =="Eu*"
                push!(y,getElem(sample,yElem).concentration/sqrt(conc(getElem(sample,"Sm"))*conc(getElem(sample,"Gd"))))
            else
                push!(y,getElem(sample,yElem).concentration/getElem(sample,yElemD).concentration)
            end
        end
    end

end

"""
Generic function for all xyplotting
"""
function xyPlot(samples::Array{Sample}; xElem::String = "", yElem::String = "", xElemD::String = "",yElemD::String="", coordAx::Integer = 0, catKey::Integer = -1, catSel = 0,saveName::String = "")
    
    pyplot()
    
    
    xLab = xElem*"("*samples[1].unit*")"
    yLab =  yElem*"("*samples[1].unit*")"

    if occursin(r"\d+", xElem)
        xLab = xElem*" age (Ma)"
    end
    #Modify this section if adding more parameters
    if coordAx == 1
        xLab = "X coordinate" 
    elseif coordAx == 2
        xLab = "Y coordinate"
    elseif !isempty(xElemD)
        xLab = xElem*"/"*xElemD
    end

    if !isempty(yElemD)
        yLab = yElem*"/"*yElemD
    end

    if catKey < 0
        x = Measurement[]
        y = Measurement[]
        sampleNames = String[]
        for sample in samples
            
            addXYVal!(x, y, sample, xElem, yElem, xElemD, yElemD,coordAx)
            push!(sampleNames,sample.name)
            
        end
        myPlot = scatter(x,y, xlabel = xLab, ylabel = yLab, hover = sampleNames,size = (PLOTX,PLOTY),markershape = MARK_SHAPE[1],markercolor = MARK_COLOR[1], xscale = XSCALING, yscale = YSCALING)
        
       
    
    #catKey = 0 means to plot a specific mineral grain, or plot based on mineral grain
    elseif catKey ==0
        x = []
        y = []
        sampleNames = []

        catX = Measurement[]
        catY = Measurement[]
        catSam = String[]
        mineralCount = 1
        if catSel > 0
            mineralCount = catSel
        end    
        labels = []
        for sample in samples
            if catSel == 0 || sample.mineralNum == catSel
                if sample.mineralNum > mineralCount

                    push!(labels,mineralCount)
                    push!(x,catX)
                    push!(y,catY)
                    push!(sampleNames, catSam)

                    catX = Measurement[]
                    catY = Measurement[]
                    catSam = String[]
                    mineralCount = sample.mineralNum

                end
               
                addXYVal!(catX, catY, sample, xElem, yElem, xElemD,yElemD, coordAx)
                push!(catSam,sample.name)
            end
            
        end
        push!(labels,mineralCount)
        push!(x,catX)
        push!(y,catY)
        push!(sampleNames, catSam)
       
        
        myPlot = scatter(x[1],y[1],xlabel = xLab, ylabel = yLab, label = labels[1], palette=palette(:roma,maxMinNum(samples)), hover = sampleNames[1],size = (PLOTX,PLOTY),markershape = MARK_SHAPE[1],xscale = XSCALING, yscale = YSCALING)
        if catSel == 0
            for i in 2:lastindex(labels)
                scatter!(myPlot, x[i],y[i], label = labels[i],hover=sampleNames[i],markershape = MARK_SHAPE[1])
            end
        end
        

    elseif catKey > 0
        x = []
        y = []
        sampleNames = []
        for cat in CATEGORIES[catKey]
            if catSel == 0 || cat == CATEGORIES[catKey][catSel]
                catX = Measurement[]
                catY = Measurement[]
                catSam = String[]

                for sample in samples
                    if sample.categories[catKey] == cat #Uses the categorization row assuming they match with CATEGORIES above

                        addXYVal!(catX, catY, sample, xElem, yElem, xElemD, yElemD, coordAx)
                        push!(catSam,sample.name)
                    end
                end
                
                push!(x,catX)
                push!(y,catY)
                push!(sampleNames, catSam)

            end
        end
        
        initialCat = 1
        if catSel > 1
            initialCat = catSel
        end

        #myPlot = DEFAULT_SCATTER
        myPlot = scatter!(dscatter(), x[1],y[1],xlabel = xLab, ylabel = yLab, label = CATEGORIES[catKey][initialCat],hover = sampleNames[1], 
        markershape = MARK_SHAPE[1],markercolor = MARK_COLOR[1], markersize = MARK_SIZE)
        if catSel == 0
            for i in 2:length(CATEGORIES[catKey])
                scatter!(myPlot, x[i],y[i], label = CATEGORIES[catKey][i], 
                hover = sampleNames[i],markershape = MARK_SHAPE[i],markercolor = MARK_COLOR[i], markersize = MARK_SIZE)
            end
        end
        
       
    
    end
    if !isempty(saveName)
        savefig(myPlot,SAVEDIR*saveName*".svg")
        savefig(myPlot,SAVEDIR*saveName*".png")
    end
    return myPlot
end

function addXYZVal!(x, y, z, sample::Sample, zElem::String , zElemD::String)
    push!(x,sample.xCoord)
    push!(y,sample.yCoord)
    if isempty(zElemD)
        push!(z,getElem(sample,zElem).concentration.val)
    else
        push!(z,(getElem(sample,zElem).concentration/getElem(sample,zElemD).concentration).val)
    end
    

end

function xyzPlot(samples::Array{Sample}; zElem::String,zElemD::String="",catKey::Integer = -1, catSel = 0, saveName::String = "")
    plotlyjs()
    #pyplot()
    zAx = zElem*"("*samples[1].unit*")"
    if !isempty(zElemD)
        zAx = zElem*"/"*zElemD
    end
    if catKey < 0
        x = Float64[]
        y = Float64[]
        z = Float64[]
        sampleNames = []
        for sample in samples
            
            addXYZVal!(x,y,z,sample,zElem,zElemD)
            push!(sampleNames,sample.name)
        end
        myPlot = scatter(x,y,z, xlabel = "X coordinate", ylabel = "Y coordinate", zlabel = zAx,hover=sampleNames,size = (PLOTX,PLOTY),markershape = MARK_SHAPE, zscale = YSCALING)
    elseif catKey ==0
        x = []
        y = []
        z = []
        sampleNames = []
        catX = Float64[]
        catY = Float64[]
        catZ = Float64[]
        catSam = String[]
        mineralCount = 1
        if catSel > 0
            mineralCount = catSel
        end    
        labels = []
        for sample in samples
            if catSel == 0 || sample.mineralNum == catSel
                if sample.mineralNum > mineralCount
                    push!(labels,mineralCount)
                    addXYZVal!(catX,catY,catZ,sample,zElem,zElemD)
                    push!(catSam,sample.name)
                    catX = Float64[]
                    catY = Float64[]
                    catZ = Float64[]
                    catSam = String[]
                    mineralCount = sample.mineralNum
                end
                
                
                push!(catX,sample.xCoord)
                push!(catY, sample.yCoord)
                push!(catZ,getElem(sample,zElem).concentration.val)
                push!(catSam,sample.name)
               
            end
            
        end
        push!(labels,mineralCount)
        push!(x,catX)
        push!(y,catY)
        push!(z, catZ)
        push!(sampleNames,catSam)
        
        
        myPlot = scatter(x[1],y[1],z[1],xlabel = "X coordinate", ylabel = "Y coordinate",hover=sampleNames[1], zlabel = zAx, label = labels[1],size = (PLOTX,PLOTY),markershape = MARK_SHAPE, zscale = YSCALING)
        if catSel == 0
            for i in 2:length(labels)
                scatter!(myPlot, x[i],y[i],z[i], label = labels[i],hover=sampleNames[i],markershape = MARK_SHAPE)
            end
        end
        

    elseif catKey > 0
        x = []
        y = []
        z = []
        sampleNames = []
        for cat in CATEGORIES[catKey]
            if catSel == 0 || cat == CATEGORIES[catKey][catSel]
                catX = Float64[]
                catY = Float64[]
                catZ = Float64[]
                catSam = String[]
                for sample in samples
                    if sample.categories[catKey] == cat #Uses the categorization row assuming they match with CATEGORIES above
                        
                        addXYZVal!(catX,catY,catZ,sample,zElem,zElemD)
                        push!(catSam,sample.name)
                    end
                end
                
                push!(x,catX)
                push!(y,catY)
                push!(z,catZ)
                push!(sampleNames,catSam)
            end
        end
        
        initialCat = 1
        if catSel > 1
            initialCat = catSel
        end
        myPlot = scatter(x[1],y[1],z[1],xlabel = "X coordinate", ylabel = "Y coordinate", zlabel = zAx,hover = sampleNames[1],label = CATEGORIES[catKey][initialCat],size = (PLOTX,PLOTY),markershape = MARK_SHAPE)
        if catSel == 0
            for i in 2:length(CATEGORIES[catKey])
                scatter!(myPlot, x[i],y[i],z[i], label = CATEGORIES[catKey][i],hover=sampleNames[i],markershape = MARK_SHAPE)
            end
        end
        
        
    
    end 
    if !isempty(saveName)
        savefig(myPlot,SAVEDIR*saveName*".svg")
        savefig(myPlot,SAVEDIR*saveName*".png")
    end
    return myPlot

end

function spiderPlot(samples::Array{Sample}, elemArray::Array{String}; norm = nothing, normName = nothing, catKey::Integer = -1, catSel::Integer = 0,saveName::String = "")
    plotlyjs()
    x = elemArray
    
        
    sampleNames = String[]
    yAx = "Concentration (",samples[1].unit,")"
    if !isnothing(norm)
        yAx = "sample/"*normName
    end
    myPlot = plot(ylabel = yAx, size = (PLOTX,PLOTY),yscale=SPIDERSCALING,ylims=(MIN_SPIDER,MAX_SPIDER))
    for sample in samples
        y = Measurement[]
        # print("\n",sample.name,": ")
        for el in elemArray
            ratio = getElem(sample,el).concentration #Not really a ratio 
            if !isnothing(norm)
                ratio = getElem(sample,el).concentration/getElem(norm,el).concentration
            end
            # print(el,": ",ratio)
            if ratio.err >= ratio.val #Messes up plotting if this isnt here with the log scale
                push!(y,NaN±NaN)
            else
                push!(y,ratio)
            end
            #print(el,": ",getElem(sample,el).concentration/getElem(norm,el).concentration )
        end
        if catKey == 0 && catSel > 0
            if sample.mineralNum == catSel
                push!(sampleNames,sample.name)
                plot!(myPlot, x, y, label = sample.name,markershape = MARK_SHAPE,linewidth=LINE_WIDTH)
            end
        elseif catKey > 0 && catSel > 0
            if sample.categories[catKey] == CATEGORIES[catKey][catSel]
                push!(sampleNames,sample.name)
                plot!(myPlot, x, y, label = sample.name,markershape = MARK_SHAPE,linewidth=LINE_WIDTH)
            end
        else
            push!(sampleNames,sample.name)
            plot!(myPlot, x, y, label = sample.name,markershape = MARK_SHAPE,linewidth=LINE_WIDTH)
        end  
        
        
    end
    if !isempty(saveName)
        savefig(myPlot,SAVEDIR*saveName*".svg")
        savefig(myPlot,SAVEDIR*saveName*".png")
    end
    return myPlot
    
    
end


#Returns array of sample objects for each item in samples that is in the rim category, these will be a ratio over each item in the core
#category that has the same mineral number
function compareRimCore(samples::Array{Sample})
    
    sampleRatios = Sample[]

    for sam in samples

        if sam.categories[1] == CATEGORIES[1][1]

            for otherSam in samples
                
                if otherSam.categories[1] == CATEGORIES[1][2] && otherSam.mineralNum == sam.mineralNum
                    name = sam.name*"/"*otherSam.name
                    unit = "ratio"
                    categories = sam.categories
                    xCoord = sam.xCoord
                    yCoord = sam.yCoord
                    mineralNum = sam.mineralNum
                    elems = Element[]
                    for i=1:length(sam.elements)
                        elName = sam.elements[i].name
                        mass = sam.elements[i].mass
                        concentration = sam.elements[i].concentration/otherSam.elements[i].concentration
                        elemRatio = Element(elName,mass,concentration)
                        push!(elems,elemRatio)
                    end
                    newSample = Sample(name,elems,unit,categories,xCoord,yCoord,mineralNum)
                    push!(sampleRatios,newSample)
                end

            end
        end

    end

    return sampleRatios

end

function getCategory(samples::Array{Sample},catKey::Integer, catSel::Integer)

    if catKey < 0
        return samples
    end
    newSams = Sample[]
    for sam in samples
        if catKey == 0 
            if sam.mineralNum == catSel
                push!(newSams,sam)
            end
        elseif sam.categories[catKey] == CATEGORIES[catKey][catSel]
            push!(newSams,sam)
        end
    end

    return newSams

end

function average(samples::Array{Sample},elem::String)
    elemArray = getElem.(samples,elem)
    elemConc = conc.(elemArray)
    concAv = sum(elemConc)/length(elemConc)
    return concAv

end

function maxMinNum(samples::Array{Sample})
    maxMinNum = 0
    for sam in samples
        if sam.mineralNum > maxMinNum
            maxMinNum = sam.mineralNum
        end
    end
    return maxMinNum

end
function partitionArray(accMin::Array{Sample},majMin::Array{Sample};catKeyAcc::Integer = -1, catSelAcc::Integer = 0,catKeyMaj::Integer=-1, catSelMaj::Integer = 0, coefFile::String = "" , saveName::String = "")
    #Partition array xy plot as in Taylor et al 2017
    #Will plot D_slope vs D_Yb where D_slope is D_Yb/D_Gd

    #First average the values in majMin of selected category
    plotlyjs()

    majCat = getCategory(majMin,catKeyMaj,catSelMaj)
    majYb =  average(majCat,"Yb")
    majGd = average(majCat,"Gd")
    majSlope = majYb/majGd

    loopCount = 1
    if catSelAcc ==0 && catKeyAcc >0
        loopCount = length(CATEGORIES[catKeyAcc])
    elseif catKeyAcc == 0 && catSelAcc ==0
        loopCount = maxMinNum(accMin)
    end
    x = []
    y = []
    labels = []
    sampleNames = []
    for i=1:loopCount
        newCatSel = catSelAcc
        if catSelAcc == 0
            newCatSel = i
        end
        accCat = getCategory(accMin,catKeyAcc,newCatSel)
        
        accYb = conc.(getElem.(accCat,"Yb"))
        accGd = conc.(getElem.(accCat,"Gd"))
        accSlope = accYb./accGd

        dYb = accYb/majYb
        dSlope = accSlope/majSlope
        push!(x,log10.(dSlope))
        push!(y,log10.(dYb))
        if catKeyAcc == 0
            push!(labels,newCatSel)
        elseif catKeyAcc > 0
            push!(labels,CATEGORIES[catKeyAcc][newCatSel])
        else
            push!(labels, saveName)
        end
        push!(sampleNames,getName.(accCat))
    end

    myPlot = scatter(x[1],y[1],xlabel = "log(D[slope])", ylabel = "log(D[Yb])", label = labels[1], 
            hover = sampleNames[1],size = (PLOTX,PLOTY),markershape = MARK_SHAPE,
            extra_plot_kwargs = KW(
                :include_mathjax => "cdn",
                :yaxis => KW(:automargin => true),
                :xaxis => KW(:domain => "auto")
                ))
            
    if catKeyAcc == 0
        scatter!(myPlot,palette=palette(:roma,maxMinNum(accMin)))
        if catSelAcc == 0
            for i in 2:length(x)
                scatter!(myPlot, x[i],y[i], label = labels[i],hover=sampleNames[i],markershape = MARK_SHAPE)
            end
        end
    
    elseif catSelAcc == 0 && catKeyAcc >0
        for i in 2:length(CATEGORIES[catKeyAcc])
            scatter!(myPlot, x[i],y[i], label = labels[i],hover=sampleNames[i],markershape = MARK_SHAPE)
        end
    end
    
    if !isempty(coefFile)
        pX, pY = importPartition(coefFile)
        scatter!(pX, pY, label = "7 kbar XP data", markershape = MARK_SHAPE)
    end
    if !isempty(saveName)
        savefig(myPlot,SAVEDIR*saveName*".svg")
        savefig(myPlot,SAVEDIR*saveName*".png")
    end
    return myPlot

end
