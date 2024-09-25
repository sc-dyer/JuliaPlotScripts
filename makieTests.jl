using GLMakie

fig = Figure(); display(GLMakie.Screen(),fig)
ax = Axis(fig[1,1])

# lines!(ax, [1,2,3],[3,2,1])
# Makie.deactivate_interaction!(ax, :rectanglezoom)
# Makie.deactivate_interaction!(ax, :scrollzoom)

# sRange = select_rectangle(ax.scene)

# on(sRange) do rect
#     x1 = rect.origin[1]
#     x2 = rect.widths[1] + x1
#     @show x1, x2
# end
lines!(ax,[0,0],[-10000,10000])