using GLMakie
using ConcaveHull

GLMakie.activate!()

fig = Figure()
ax = Axis(fig[1, 1], aspect = 1.0)

points = [[0.0,-1.0],[0.0,1.0],[10.0,5.0],[10.0,-5.0],[8.0,2.0],[8.0,-2.0],[4.0,2.0],[4.0,-2.0]]
x = [p[1] for p in points]
y = [p[2] for p in points]

hull = concave_hull(points,1)

scatter!(ax,x,y)

vx = [v[1] for v in hull.vertices]
vy = [v[2] for v in hull.vertices]

poly!(ax,vx,vy,alpha = 0.5)
fig