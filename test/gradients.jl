using ForwardDiff


@testset "Profile forward sensitivities." begin
    runCkineAS(tps, rxntfR, ones(gcSolver.Nspecies), ones(length(tps)))

    println("runCkineAS")
    @time runCkineAS(tps, rxntfR, ones(gcSolver.Nspecies), ones(length(tps)))
end
