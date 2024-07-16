#Simple example plots to illustrate cahnges in gibbs energy

using Plots
using LaTeXStrings

PLOTX = 800
PLOTY = 600

#MgO-SiO2 system binary with XSiO2 on the X axis
#Projecting from H2O

quartz = 1.0
periclase = 0.0
forsterite = 1.0/3.0#Mg2SiO4
enstatite = 1.0/2.0#MgSiO3
talc = 4.0/7.0#Mg3Si4O10(OH)2
serpentine = 2.0/5.0 #Mg3Si2O5(OH)4

X = [periclase,forsterite,enstatite,quartz]
labels = ["Per","Fo","En","Qz"]
G = [-609.43,-2200.94,-3129.72,-923.07]

pyplot()
#scatter(X,G, xlabel = "X", ylabel = "G", ms = 10,mc = "red",legend=false,ylim=[-5000,0],size = [PLOTX,PLOTY],grid = false,xlabelfontsize = 24, ylabelfontsize = 24, ytickfontsize = 0,ytickdirection = :none, xtickdirection = :out,xtickfontsize = 16,reuse = false)

dG = [-100,-40,30,-50]

G = G + 42*dG

scatter(X,G, xlabel = "X", ylabel = "G", ms = 10,mc = "red",legend=false,ylim=[-5000,0],size = [PLOTX,PLOTY],grid = false,xlabelfontsize = 24, ylabelfontsize = 24, ytickfontsize = 0,ytickdirection = :none, xtickdirection = :out,xtickfontsize = 16,reuse = false)
anX = X 
anG = G .+ 80
annotate!(anX,anG,labels,:bottom)


P = range(-1,50)
T = (10000/18.08).*P .+ (-555.82)
#P = (18.08).*x .+ (-555.82-273.15)
plot(T,P,size = [PLOTX,PLOTY],xlabel = "T (K)",ylabel = "P (GPa)" ,xlim=[0,6000],legend = false,linewidth = 4, linecolor = "red",ylim = [0,50],grid = false,xlabelfontsize = 24, ylabelfontsize = 24,framestyle = :box, tickfontsize = 16, tickdirection =:in)
and = [-2588.67e3,92.7,5.153]
sill = [-2585.79e3,95.4,4.986]
ky = [-2592.97e3,83.6,4.414]
P = range(-1,10)
r1 = (sill[3]-and[3])/(sill[2]-and[2]).*P*1000 .+ ((sill[1]-and[1])/(sill[2]-and[2])-273.15)

r2 = (ky[3]-and[3])/(ky[2]-and[2]).*P*1000 .+ ((ky[1]-and[1])/(ky[2]-and[2])-273.15)
r3 = (ky[3]-sill[3])/(ky[2]-sill[2]).*P*1000 .+ ((ky[1]-sill[1])/(ky[2]-sill[2])-273.15)

plot(r1,P,size = [PLOTX,PLOTY],xlabel = "T (C)",ylabel = "P (kbar)" ,xlim=[400,1000],legend = false,linewidth = 4, linecolor = "red",ylim = [0,10],grid = false,xlabelfontsize = 24, ylabelfontsize = 24,framestyle = :box, tickfontsize = 16, tickdirection =:in)
plot!(r2,P, linewidth = 4, linecolor = "blue")
plot!(r3,P,linewidth = 4, linecolor = "green")
