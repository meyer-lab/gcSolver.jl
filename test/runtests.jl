using Test
using gcSolver
using ForwardDiff
using LinearAlgebra
using Random
Random.seed!(123)

rxntfR = exp.(randn(gcSolver.Nparams))
rxntfR[51] = tanh(rxntfR[51])

tps = [0.0, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0]

include("model.jl")
include("gradients.jl")
