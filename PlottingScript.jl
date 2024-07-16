include("ICPMS_Plots.jl")

Zrn20SD06_te = importData("../../LAICPMS/2022_03_03/220303_Sabastien_Zrc_AP/20SD06Zrn_NoOutliersB.csv",3)
Zrn20SD06_upb = importData("../../LAICPMS/2022_09_02_LASS/20SD06_Zrn_NoOutliers.csv",4)
Zrn21SD68 = importData("../../LAICPMS/2022_09_02_LASS/21SD68_Zrn_NoOutliers.csv",4)
Zrn22SD55B = importData("../../LAICPMS/2022_09_02_LASS/22SD55B_Zrn_NoOutliers.csv",4)


XSCALING = :log10
YSCALING = :identity
plot1 = xyPlot(Zrn20SD06_te;xElem = "Th", xElemD = "U", yElem = "Lu", catKey = 1,saveName = "Zrn20SD06-1_Lu")
plot2 = xyPlot(Zrn20SD06_upb;xElem = "Th", xElemD = "U", yElem = "Lu", catKey = 1,saveName = "Zrn20SD06-2_Lu")
plot3 = xyPlot(Zrn20SD06_te;xElem = "Th", xElemD = "U", yElem = "Hf", catKey = 1,saveName = "Zrn20SD06-1_Hf")
plot4 = xyPlot(Zrn20SD06_upb;xElem = "Th", xElemD = "U", yElem = "Hf", catKey = 1,saveName = "Zrn20SD06-2_Hf")
XSCALING = :identity
plot5 = xyPlot(Zrn20SD06_upb;xElem = "Pb207-Pb206", yElem = "Th",yElemD = "U", catKey = 1,saveName = "Zrn20SD06-2_UPb")

XSCALING = :log10

plot6 = xyPlot(Zrn21SD68;xElem = "Th", xElemD = "U", yElem = "Lu", catKey = 1,saveName = "Zrn20SD68_Lu")
plot7 = xyPlot(Zrn21SD68;xElem = "Th", xElemD = "U", yElem = "Hf", catKey = 1,saveName = "Zrn20SD68_Hf")
XSCALING = :identity
plot8 = xyPlot(Zrn21SD68;xElem = "Pb207-Pb206", yElem = "Th",yElemD = "U", catKey = 1,saveName = "Zrn20SD68_UPb")

XSCALING = :log10
plot9 = xyPlot(Zrn22SD55B;xElem = "Th", xElemD = "U", yElem = "Lu", catKey = 1)
plot10 = xyPlot(Zrn22SD55B;xElem = "Th", xElemD = "U", yElem = "Hf", catKey = 1)
XSCALING = :identity
plot11 = xyPlot(Zrn22SD55B;xElem = "Pb207-Pb206", yElem = "Th",yElemD = "U", catKey = 1)

