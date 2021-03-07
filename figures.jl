using Pkg
Pkg.activate(".")
Pkg.add("Gadfly")
Pkg.add("Plots")
Pkg.update()

println("Making Figure 1")
include("figures/figureJ1.jl")
figureJ1()

println("Making Figure 2")
include("figures/figureJ2.jl")
figureJ2()
