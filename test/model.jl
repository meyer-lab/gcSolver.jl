
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

@testset "Model tests." begin
    @testset "Reaction model mass conservation." begin
        dy = zeros(gcSolver.halfL)

        gcSolver.dYdT(dy, copy(dy), rxntfR, ones(3))

        # Check for conservation of each surface receptor
        assertConservation(dy)
    end


    @testset "Full model mass conservation." begin
        rr = copy(rxntfR)
        rr[24:end] .= 0.0
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
        rr = copy(rxntfR)
        rr[1:3] .= 0.0

        out = gcSolver.runCkine([0.0, 1.0], rr)

        @test all(out .>= 0.0)
        @test norm(out[1, :] - out[2, :]) < 1.0e-4
    end


    @testset "Equilibrium." begin
        out = runCkine([100000.0], rxntfR)

        dy = ones(gcSolver.Nspecies)

        gcSolver.fullDeriv(dy, vec(out), rxntfR, 0.0)

        @test all(out .>= 0.0)

        @test norm(dy) < 1.0e-6
    end


    @testset "Detailed balance using no-trafficking model." begin
        prob = gcSolver.runCkineSetup([1000000.0], deepcopy(rxntfR))
        prob.p[20:30] .= 0.0 # Set everything outside of binding reactions to 0

        out = vec(runCkine([1000000.0], prob))

        J = ForwardDiff.jacobian((y, x) -> gcSolver.fullDeriv(y, x, prob.p, 0.0), ones(gcSolver.Nspecies), out)

        # Slice out just the surface species
        GK = J[1:(2gcSolver.halfL), 1:(2gcSolver.halfL)] * diagm(vec(out[1:(2gcSolver.halfL)]))

        @test norm(GK - transpose(GK)) < 1.0e-9
    end


    @testset "Test that there is at least 1 non-zero species at T=0." begin
        temp = runCkine(tps, rxntfR)
        @test any(temp[1, :] .> 0.0)
    end
end
