println("Testing...")

using Distributions
using ForwardDiff
using BenchmarkTools
using gcSolver

### Check that runCkine can run
params = [rand(LogNormal(0.1, 0.25)) for i=1:gcSolver.Nparams]
params[20] = tanh(params[20])

IL2params = [rand(LogNormal(0.1, 0.25)) for i=1:gcSolver.NIL2params]

tps = [0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0]

#@benchmark runCkine(tps, params, false)

out = runCkine(tps, params, false)
outTwo = runCkine(tps, params, false)

@assert out == outTwo


fout = ForwardDiff.jacobian(x -> runCkine(tps, x, false)[1], params)

#print(out[1])
#print("\n\n\n\n")

#print(fout)

## Run two times and check solutions are identical with/without sensitivity, pretreatment, IL2param