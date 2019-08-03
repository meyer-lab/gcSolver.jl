
include("model.jl")

#using Distributions
#using ForwardDiff
#using BenchmarkTools
#using gcSolver

### Check that runCkine can run

#@benchmark runCkine(tps, rxntfR, false)


#fout = ForwardDiff.jacobian(x -> runCkine(tps, x, false)[1], rxntfR)

## Run two times and check solutions are identical with/without sensitivity, pretreatment, IL2param