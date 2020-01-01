using Test
using OrdinaryDiffEq
using Profile
using gcSolver

rxntfR = exp.(randn(gcSolver.Nparams))
rxntfR[20] = tanh(rxntfR[20])

IL2params = exp.(randn(gcSolver.NIL2params))

surface = ones(eltype(rxntfR), 21)
endosome = copy(surface)
ILs = zeros(eltype(rxntfR), gcSolver.Nlig)
trafP = zeros(eltype(rxntfR), 13)

tps = [0.0, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0]

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
    
    gcSolver.fullDeriv(dy, copy(dy), (rr, surface, endosome, trafP, ILs), 0.0)
    
    # Check for conservation of each surface receptor
    assertConservation(dy)
    # Check for conservation of each endosomal receptor
    assertConservation(dy[gcSolver.halfL+1:end])
end


@testset "Reproducibility." begin
    output = runCkine(tps, rxntfR)

    @test ndims(output) == 2
    @test output == runCkine(tps, rxntfR)
    @test runCkine(tps, IL2params) == runCkine(tps, IL2params)
end


@testset "Steady-state at t=0." begin
    gcSolver.fullParam!(rxntfR, surface, endosome, trafP, ILs)
    out = gcSolver.solveAutocrine(trafP)
    gcSolver.fullParam!(IL2params, surface, endosome, trafP, ILs)
    IL2out = gcSolver.solveAutocrine(trafP)

    rr = copy(rxntfR)
    IL2rr = copy(IL2params)
    rr[1:6] .= 0.0
    IL2rr[1] = 0.0

    dy = ones(gcSolver.Nspecies)
    IL2dy = ones(gcSolver.Nspecies)

    gcSolver.fullDeriv(dy, out, (rr, surface, endosome, trafP, ILs), 0.0)
    gcSolver.fullDeriv(IL2dy, IL2out, (IL2rr, surface, endosome, trafP, ILs), 0.0)

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

    gcSolver.fullDeriv(dy, out[1, :], (rxntfR, surface, endosome, trafP, ILs), 0.0)
    gcSolver.fullDeriv(IL2dy, IL2out[1, :], (IL2params, surface, endosome, trafP, ILs), 0.0)

    println("runCkineSS")
    @time outSS = runCkineSS(rxntfR)

    @test all(out .>= 0.0)
    @test all(IL2out .>= 0.0)

    @test isapprox(sum(abs.(dy)), 0.0, atol=1.0e-6)
    @test isapprox(sum(abs.(IL2dy)), 0.0, atol=1.0e-6)
end


@testset "Make sure no endosomal species are found when endo=0." begin
    rxntfRR = copy(rxntfR)
    rxntfRR[18:19] .= 0.0  # set endo and activeEndo to 0.0

    yOut = runCkine(tps, rxntfRR)

    @test all(isapprox(sum(abs.(yOut[:, 29:end])), 0.0, atol=1.0e-6))
end


@testset "Test that there is at least 1 non-zero species at T=0." begin
    temp = runCkine(tps, rxntfR)
    @test any(temp[1, :] .> 0.0)
end


@testset "Benchmark." begin
    println("fullDeriv")
    @time gcSolver.fullDeriv(zeros(gcSolver.Nspecies), ones(gcSolver.Nspecies), (rxntfR, surface, endosome, trafP, ILs), 0.0)
    println("fullDeriv IL2")
    @time gcSolver.fullDeriv(zeros(gcSolver.Nspecies), ones(gcSolver.Nspecies), (IL2params, surface, endosome, trafP, ILs), 0.0)

    println("Default runCkine")
    @time runCkine(tps, rxntfR)
    println("Default runCkine IL2")
    @time runCkine(tps, IL2params)

    for ii in 1:10
        @profile runCkine(tps, rxntfR)
    end
    Profile.print(noisefloor=5.0)
end
