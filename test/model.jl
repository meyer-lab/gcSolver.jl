
# Assert the conservation of species throughout the experiment.
function assertConservation(y)
    # These are all 0-indexed, so add 1
    consList = [
        [1, 4, 5, 7, 8, 11, 12, 14, 15],  # IL2Rb
        [0, 3, 5, 6, 8],  # IL2Ra
        [9, 10, 12, 13, 15],  # IL15Ra
        [16, 17, 18],  # IL7Ra
        [2, 6, 7, 8, 13, 14, 15, 18],
    ] #gc

    # Check for conservation of species sum
    for ii in range(1, stop = length(consList))
        diff = sum(y[consList[ii] .+ 1])
        @test isapprox(diff, 0.0, atol = 1.0e-12)
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
    rr[44:end] .= 0.0
    dy = zeros(gcSolver.Nspecies)

    gcSolver.fullDeriv(dy, copy(dy), rr, 0.0)

    # Check for conservation of each surface receptor
    assertConservation(dy)
    # Check for conservation of each endosomal receptor
    assertConservation(dy[(gcSolver.halfL + 1):end])
end


@testset "Reproducibility." begin
    output = runCkine(tps, rxntfR)

    @test ndims(output) == 2
    @test output == runCkine(tps, rxntfR)
end


@testset "Steady-state at t=0." begin
    out = gcSolver.solveAutocrine(rxntfR)

    rr = copy(rxntfR)
    rr[1:(gcSolver.Nlig)] .= 0.0

    dy = ones(gcSolver.Nspecies)

    gcSolver.fullDeriv(dy, out, rr, 0.0)

    @test all(out .>= 0.0)
    @test isapprox(sum(abs.(dy)), 0.0, atol = 1.0e-12)
end


@testset "Equilibrium." begin
    out = runCkine([100000.0], rxntfR)

    dy = ones(gcSolver.Nspecies)

    gcSolver.fullDeriv(dy, out[1, :], rxntfR, 0.0)

    @test all(out .>= 0.0)

    @test isapprox(sum(abs.(dy)), 0.0, atol = 1.0e-6)
end


@testset "Detailed balance using no-trafficking model." begin
    rxntfRR = copy(rxntfR)
    rxntfRR[40:41] .= 0.0  # set endo and activeEndo to 0.0
    out = vec(runCkine([1000000.0], rxntfRR))

    J = ForwardDiff.jacobian((y, x) -> gcSolver.fullDeriv(y, x, rxntfRR, 0.0), ones(gcSolver.Nspecies), out)
    GK = J * diagm(vec(out))

    @test norm(GK - transpose(GK)) < 1.0e-9
end


@testset "Make sure no endosomal species are found when endo=0." begin
    rxntfRR = copy(rxntfR)
    rxntfRR[40:41] .= 0.0  # set endo and activeEndo to 0.0

    yOut = runCkine(tps, rxntfRR)

    @test all(isapprox(sum(abs.(yOut[:, 20:end])), 0.0, atol = 1.0e-6))
end


@testset "Test that there is at least 1 non-zero species at T=0." begin
    temp = runCkine(tps, rxntfR)
    @test any(temp[1, :] .> 0.0)
end


@testset "Benchmark." begin
    println("fullDeriv")
    @time gcSolver.fullDeriv(zeros(gcSolver.Nspecies), ones(gcSolver.Nspecies), rxntfR, 0.0)

    println("runCkine")
    @time runCkine(tps, rxntfR)
end
