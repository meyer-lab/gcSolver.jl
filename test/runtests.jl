using Test
using Profile
using gcSolver
using ForwardDiff
using LinearAlgebra

rxntfR = ones(gcSolver.Nparams) * 0.1

tps = [0.0, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0]

include("model.jl")
include("varprop.jl")
