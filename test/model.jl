using gcSolver

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
        @test isapprox(diff, 0.0, atol=1.0e-12) "$ii species not conserved: $diff"
    end
end


@testset "Full model mass conservation." begin
    rxntfR = copy(rxntfR)
    rxntfR[18:end] .= 0.0
    y0 = ones(gcSolver.Nspecies)
    dy = ones(gcSolver.Nspecies)
    
    gcSolver.fullDeriv(dy, y0, rxntfR, 0.0)
    
    # Check for conservation of each surface receptor
    assertConservation(dy)
    # Check for conservation of each endosomal receptor
    assertConservation(dy[gcSolver.halfL+1:end])
end
