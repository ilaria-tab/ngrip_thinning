import Pkg
Pkg.add("Plots")

using Plots #or StatsPlots

x = range(0, 10, length=100)
y = sin.(x)
plot(x, y)

