
@testset "Reasonable return from varprop function." begin
    sigma = Matrix{Int}(I, length(rxntfR), length(rxntfR))

    retval = gcSolver.runCkineVarPorp(tps, rxntfR, sigma)

    @time gcSolver.runCkineVarPorp(tps, rxntfR, sigma)

    @test size(retval) == (length(tps), length(tps))
end

