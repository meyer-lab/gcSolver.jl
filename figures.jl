using Pkg
Pkg.activate(".")

println("Making Figure 1")
include("figures/figureJ1.jl")
figureJ1()

println("Making Figure 2")
include("figures/figureJ2.jl")
figureJ2()

println("Making Figure 3")
include("figures/figureJ3.jl")
figureJ3()

println("Making Figure 4")
include("figures/figureJ4.jl")
figureJ4()

println("Making Figure 5")
include("figures/figureJ5.jl")
figureJ5()

println("Making Figure 6")
include("figures/figureJ6.jl")
figureJ6()
