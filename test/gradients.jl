using ForwardDiff


@testset "Profile forward sensitivities." begin
	runCkineAS(tps, rxntfR, ones(gcSolver.Nspecies), ones(length(tps)))

    println("runCkineAS")
    @time runCkineAS(tps, rxntfR, ones(gcSolver.Nspecies), ones(length(tps)))

    @profile runCkineAS(tps, rxntfR, ones(gcSolver.Nspecies), ones(length(tps)))
    Profile.print(noisefloor = 5.0)
end