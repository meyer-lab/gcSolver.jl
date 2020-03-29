using gcSolver
using Plots
using OrdinaryDiffEq
using DiffEqDevTools

rxntfR = ones(gcSolver.Nparams) * 0.2

prob = gcSolver.runCkineSetup([500.0], rxntfR)

setups = [Dict(:alg=>AutoTsit5(Rodas4(autodiff=false))),
          Dict(:alg=>AutoTsit5(Rodas4P(autodiff=false))),
          # Dict(:alg=>Rodas4P(autodiff=false)),
          Dict(:alg=>AutoTsit5(Rodas5(autodiff=false))),
          # Dict(:alg=>Rosenbrock23(autodiff=false)), # slow
          # Dict(:alg=>TRBDF2(autodiff=false)), # slow
          # Dict(:alg=>AutoTsit5(RadauIIA5(autodiff=false))),
          # Dict(:alg=>ABDF2(autodiff=false)), # slow
          # Dict(:alg=>QNDF(autodiff=false)), # slow
          # Dict(:alg=>Exprb43(autodiff=false)), # slow
          # Dict(:alg=>Exprb32(autodiff=false)), # slow
]

test_sol = solve(prob, AutoTsit5(Rodas5()), reltol=1.0e-16, abstol=1.0e-16)

tols = 1.0 ./ 10.0 .^ LinRange(4, 12, 30)

wp = WorkPrecisionSet(prob, tols, tols, setups; save_everystep=false, appxsol=test_sol, maxiters=Int(1e6))

plot(wp)
