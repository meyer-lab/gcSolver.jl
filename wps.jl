using gcSolver
using Plots
using OrdinaryDiffEq
using DiffEqDevTools

rxntfR = ones(gcSolver.Nparams) * 0.2

prob = gcSolver.runCkineSetup([500.0], rxntfR)

setups = [Dict(:alg => AutoTsit5(Rodas5())), Dict(:alg => AutoTsit5(Kvaerno5())), Dict(:alg => AutoVern9(Rodas5())), Dict(:alg => AutoVern9(KenCarp5())), Dict(:alg => AutoVern9(Kvaerno5()))]

test_sol = solve(prob, AutoTsit5(Rodas5()), reltol = 1.0e-16, abstol = 1.0e-16)

tols = 1.0 ./ 10.0 .^ LinRange(4, 12, 30)

wp = WorkPrecisionSet(prob, tols, tols, setups; save_everystep = false, appxsol = test_sol, maxiters = Int(1e6))

plot(wp)
