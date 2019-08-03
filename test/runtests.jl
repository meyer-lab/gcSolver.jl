println("Testing...")

using Distributions
using ForwardDiff
using BenchmarkTools
using gcSolver

rxntfR = [rand(LogNormal(0.1, 0.25)) for i=1:gcSolver.Nparams]
rxntfR[20] = tanh(rxntfR[20])

include("model.jl")

### Check that runCkine can run
IL2params = [rand(LogNormal(0.1, 0.25)) for i=1:gcSolver.NIL2params]

tps = [0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0]

#@benchmark runCkine(tps, rxntfR, false)

out = runCkine(tps, rxntfR, false)
outTwo = runCkine(tps, rxntfR, false)

@assert out == outTwo


fout = ForwardDiff.jacobian(x -> runCkine(tps, x, false)[1], rxntfR)

## Run two times and check solutions are identical with/without sensitivity, pretreatment, IL2param