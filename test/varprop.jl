
@testset "Reasonable return from varprop function." begin
    sigma = Matrix{Int}(I, 5, 5)

    retval = runCkineVarProp(tps, rxntfR, sigma)

    @test length(retval) == length(tps)
    @test all(retval .>= 0.0)
end


@testset "Reasonable return from Hessian functions." begin
    retval = gcSolver.runCkineHessian(tps, rxntfR)

    # TODO: Check output
end


@testset "Run fitting." begin
    p = gcSolver.getUnkVec()
    outtG = similar(p)

    outt = gcSolver.resids(p)
    ForwardDiff.gradient!(outtG, gcSolver.resids, p)

    @test isfinite(outt)
    @test all(isfinite.(outtG))

    # minOut = gcSolver.runFit(itern = 10)
    # @test all(isfinite.(minOut))
    # @test length(p) == length(minOut)
end
