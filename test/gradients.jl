using ForwardDiff


@testset "Profile forward sensitivities." begin
    println("runCkineAS")
    @time runCkineAS(tps, rxntfR, ones(gcSolver.Nspecies))

    @profile runCkineAS(tps, rxntfR, ones(gcSolver.Nspecies))
    Profile.print(noisefloor = 2.0)
end