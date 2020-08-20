
@testset "Reasonable return from varprop function." begin
    sigma = Matrix{Int}(I, 5, 5)

    retval = runCkineVarProp(tps, rxntfR, sigma)

    @test length(retval) == length(tps)
    @test all(retval .>= 0.0)
end


@testset "Run fitting." begin
    p = gcSolver.getUnkVec()
    outtG = similar(p)

    outt = gcSolver.resids(p)
    ForwardDiff.gradient!(outtG, gcSolver.resids, p)

    @time ForwardDiff.gradient!(outtG, gcSolver.resids, p)

    @test isfinite(outt)
    @test all(isfinite.(outtG))

    minOut = gcSolver.runFit(itern = 4)
    @test all(isfinite.(minOut))
    @test length(p) == length(minOut)
end
