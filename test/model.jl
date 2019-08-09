using Test
using Distributions
using BenchmarkTools
using OrdinaryDiffEq
using Profile
using gcSolver

rxntfR = [rand(LogNormal(0.1, 0.25)) for i=1:gcSolver.Nparams]
rxntfR[20] = tanh(rxntfR[20])

IL2params = [rand(LogNormal(0.1, 0.25)) for i=1:gcSolver.NIL2params]

tps = [0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0]

# Assert the conservation of species throughout the experiment.
function assertConservation(y)
    # These are all 0-indexed, so add 1
    consList = [[1, 4, 5, 7, 8, 11, 12, 14, 15],  # IL2Rb
                [0, 3, 5, 6, 8],  # IL2Ra
                [9, 10, 12, 13, 15],  # IL15Ra
                [16, 17, 18],  # IL7Ra
                [19, 20, 21],  # IL9R
                [22, 23, 24],  # IL4Ra
                [25, 26, 27],  # IL21Ra
                [2, 6, 7, 8, 13, 14, 15, 18, 21, 24, 27]] #gc

    # Check for conservation of species sum
    for ii in range(1, stop=length(consList))
        diff = sum(y[consList[ii] .+ 1])
        @test isapprox(diff, 0.0, atol=1.0e-12)
    end
end


@testset "Reaction model mass conservation." begin
    dy = zeros(gcSolver.halfL)
    
    gcSolver.dYdT(dy, copy(dy), rxntfR, ones(gcSolver.Nlig))
    
    # Check for conservation of each surface receptor
    assertConservation(dy)
end


@testset "Full model mass conservation." begin
    rr = copy(rxntfR)
    rr[18:end] .= 0.0
    dy = zeros(gcSolver.Nspecies)
    
    gcSolver.fullDeriv(dy, copy(dy), rr, 0.0)
    
    # Check for conservation of each surface receptor
    assertConservation(dy)
    # Check for conservation of each endosomal receptor
    assertConservation(dy[gcSolver.halfL+1:end])
end


@testset "Reproducibility." begin
    @test runCkine(tps, rxntfR) == runCkine(tps, rxntfR)
    @test runCkine(tps, IL2params) == runCkine(tps, IL2params)
    @test runCkineS(tps, rxntfR) == runCkineS(tps, rxntfR)
    @test runCkineS(tps, IL2params) == runCkineS(tps, IL2params)
end


@testset "Steady-state at t=0." begin
    out = gcSolver.solveAutocrine(gcSolver.fullParam(rxntfR)[4])
    IL2out = gcSolver.solveAutocrine(gcSolver.fullParam(IL2params)[4])

    rr = copy(rxntfR)
    IL2rr = copy(IL2params)
    rr[1:6] .= 0.0
    IL2rr[1] = 0.0

    dy = ones(gcSolver.Nspecies)
    IL2dy = ones(gcSolver.Nspecies)

    gcSolver.fullDeriv(dy, out, rr, 0.0)
    gcSolver.fullDeriv(IL2dy, IL2out, IL2rr, 0.0)

    @test all(out .>= 0.0)
    @test all(IL2out .>= 0.0)

    @test isapprox(sum(abs.(dy)), 0.0, atol=1.0e-12)
    @test isapprox(sum(abs.(IL2dy)), 0.0, atol=1.0e-12)
end


@testset "Equilibrium." begin
    out = runCkine([100000.0], rxntfR)
    IL2out = runCkine([100000.0], IL2params)

    dy = ones(gcSolver.Nspecies)
    IL2dy = ones(gcSolver.Nspecies)

    gcSolver.fullDeriv(dy, out[1], rxntfR, 0.0)
    gcSolver.fullDeriv(IL2dy, IL2out[1], IL2params, 0.0)

    @test all(out[1] .>= 0.0)
    @test all(IL2out[1] .>= 0.0)

    @test isapprox(sum(abs.(dy)), 0.0, atol=1.0e-6)
    @test isapprox(sum(abs.(IL2dy)), 0.0, atol=1.0e-6)
end


@testset "Benchmark." begin
    @time gcSolver.fullDeriv(zeros(gcSolver.Nspecies), ones(gcSolver.Nspecies), rxntfR, 0.0)
    
    println("TRBDF2")
    @time runCkine(tps, rxntfR; alg=TRBDF2())
    println("Rosenbrock23")
    @time runCkine(tps, rxntfR; alg=Rosenbrock23())
    println("Kvaerno5")
    @time runCkine(tps, rxntfR; alg=Kvaerno5())
    println("KenCarp4")
    @time runCkine(tps, rxntfR; alg=KenCarp4())
    println("Rodas4P")
    @time runCkine(tps, rxntfR; alg=Rodas4P())
    println("Rodas5")
    @time runCkine(tps, rxntfR; alg=Rodas5())
    println("ImplicitEuler")
    @time runCkine(tps, rxntfR; alg=ImplicitEuler())
    println("Default IL-2")
    @time runCkine(tps, IL2params)
    println("Default runCkineS")
    @time runCkineS(tps, rxntfR)

    for ii in 1:10
        @profile runCkine(tps, rxntfR)
    end
    Profile.print(noisefloor=5.0)
end
