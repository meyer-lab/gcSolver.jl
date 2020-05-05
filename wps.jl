using gcSolver
using Plots
using OrdinaryDiffEq
using DiffEqDevTools

function makeWPS()
    unkVec1 = getUnkVec()

    ILs = [84.0, 0.0, 0.0]

    exprDF = getExpression()
    expr = exprDF.Treg

    rxntfR = fitParams(ILs, unkVec1, expr)

    prob = gcSolver.runCkineSetup([500.0], rxntfR)

    setups = [Dict(:alg => AutoTsit5(Rodas4())), Dict(:alg => AutoTsit5(Rodas4P())), Dict(:alg => AutoTsit5(Rodas5()))]

    test_sol = solve(prob, AutoTsit5(Rodas5()), reltol = 1.0e-16, abstol = 1.0e-16)

    tols = 1.0 ./ 10.0 .^ LinRange(4, 12, 30)

    wp = WorkPrecisionSet(prob, tols, tols, setups; save_everystep = false, appxsol = test_sol, maxiters = Int(1e6))

    plot(wp)
end
