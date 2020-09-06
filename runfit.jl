using Distributed
addprocs(32; exeflags="--project")

@everywhere using Pkg
@everywhere Pkg.instantiate()
@everywhere using gcSolver

gcSolver.runFit()
