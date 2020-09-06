using Distributed
addprocs(16; exeflags="--project")

println("Loading packages")
@everywhere using Pkg
@everywhere Pkg.instantiate()
@everywhere using gcSolver

println("Starting fitting")
gcSolver.runFit()
