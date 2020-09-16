
@testset "Reasonable return from varprop function." begin
    sigma = Matrix{Int}(I, 5, 5)

    retval = runCkineVarProp(tps, rxntfR, sigma)

    @test length(retval) == length(tps)
    @test all(retval .>= 0.0)
end


@testset "Run fitting." begin
    p = gcSolver.getUnkVec()
    outtG = similar(p)

    outt = gcSolver.resids(p, false)
    ForwardDiff.gradient!(outtG, gcSolver.resids, p, false)

    @time ForwardDiff.gradient!(outtG, gcSolver.resids, p, false)

    @test isfinite(outt)
    @test all(isfinite.(outtG))
end
