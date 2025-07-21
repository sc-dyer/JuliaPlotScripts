#myColours = ["#FF8CA6","#6D8FFF","#F9C52A","#A9EFA5", "#E6A7FB","#000000","","#008000","#0000FF","#4B0082","#FF8C00"]

myColours = ["#000000","#E69F00","#56B4E9","#009E73","#0072B2","#D55E00","#CC79A7"]
colourscheme = myColours
myColours2 = ["#E69F00","#56B4E9","#009E73","#0072B2","#D55E00","#CC79A7","#000000"]#"#F0E442"
myColours3 = ["#E69F00","#0072B2","#CC79A7","#56B4E9","#009E73","#D55E00","#000000"]
schemes = [myColours,myColours2,myColours3]
colourscheme = schemes[colourchoice]
myTheme = Theme(
    fonts = (; regular = "B612", bold = "B612-Bold"),
    fontsize = 14,
    palette = (;color =colourscheme,marker = [:circle,:rect,:diamond,:cross, :xcross,:utriangle,:star5,:star8, ],linestyle=[:solid,(:dot,:dense)]),
    Poly = (;cycle = [:color]),
    Scatter = (;cycle = [:color,:marker],markersize = 15),
    BoxPlot = (;cycle = [:color]),
    Lines = (;cycle = [:color,:linestyle]),
    ScatterLines = (;cycle = [:color,:linestyle]),
    TernaryScatter = (;cycle=[:color]),
    Errorbars = (;cycle = [:color]),
    Axis = ( 
        xlabelsize = 16,
        ylabelsize = 16,
        xminorticksvisible = true,
        yminorticksvisible = true,
        xticksmirrored = true,
        yticksmirrored = true,
        xtickalign = 1,
        ytickalign = 1,
        xminortickalign = 1,
        yminortickalign = 1,
        xticksize = 10,
        xtickwidth = 2,
        yticksize = 10,
        ytickwidth = 2,
        xminorticksize = 5,
        xminortickwidth = 2,
        yminorticksize = 5,
        yminortickwidth = 2,
        xgridvisible = false,
        ygridvisible = false
        )
)
